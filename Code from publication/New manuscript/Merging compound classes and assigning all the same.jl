include("All functions.jl")

#Import the updated raw data
folder_path = "R:\\PHD2024TH-Q6813\\Research files\\RepoRT updated 08042025\\raw_data"

Final_table = DataFrame(SMILES = String[], InChi = String[], LC_mode = String7[], Column_name = String[], CID = String[], Retention_factors = Float64[], Modifier = Float64[], pH = Float64[], study = Int[])

#These are the fingerprint names for each corresponding column
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")
#Going through each file of the raw data and extracting the necessary information
for i = 1:length(readdir(folder_path))
    #Locating the current file
    studies = readdir(folder_path)[i]
    current_file_raw = joinpath(folder_path, readdir(folder_path)[i])
    println(i)
    #Getting the SMILES, InChI, retention times and CID for each chemical in the current file
    rt_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i])", "rtdata.tsv"],"_"))
    data_compounds = CSV.read(rt_path, DataFrame)
    Compound_SMILES = data_compounds[:,7]
    InChI = data_compounds[:,8]
    RT = data_compounds[:,4]
    CID = data_compounds[:,5]

    #Getting the retention factors and fraction of organic mobile phase at elution if available
    gradient_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i])", "gradient.tsv"],"_"))
    gradient = CSV.read(gradient_path, DataFrame)
    local retention_factors
    Modifier = []
    try 
        run_time = gradient[end,1]
        retention_factors = RT./run_time
        for rtfs in retention_factors
            push!(Modifier, interpolate_B_modifier(rtfs, gradient))
        end

    catch
         retention_factors = repeat(["No gradient info"], length(RT))
         Modifier = repeat(["No gradient info"], length(RT))
    end

    #Getting the LC mode info for the data set
    info_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i])", "info.tsv"],"_"))
    info = CSV.read(info_path, DataFrame)
    #Check if the column info is missing in the meta data
    LC_mode = info[1,3]
    if ismissing(LC_mode) || LC_mode == "NA"

        LC_mode = repeat(["No mode info"], length(RT))

    else

        LC_mode = repeat([LC_mode], length(RT))

    end

    #Getting the column name for each method
    meta_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i])", "metadata.tsv"],"_"))
    meta = CSV.read(meta_path, DataFrame)
    Column_name = meta[1,2]
    pHA = meta[!, Symbol("eluent.A.pH")][1]
    pHB = meta[!, Symbol("eluent.B.pH")][1]
    #Getting all the data in one table to append all the other datasets in this format
    current_table = DataFrame(SMILES = Compound_SMILES, InChi = InChI, LC_mode = LC_mode, Column_name = Column_name, CID = CID, Retention_factors = retention_factors, Modifier = Modifier, pH = pHA, study = studies)
    Final_table = vcat(Final_table, current_table)
end

#Removing some odd entries and entries with missing retention times
Final_table_unique = remove_missing_and_outliers(Final_table)

#Getting indices for RPLC entries
indices_RPLC = findall(row -> occursin("RP", row.LC_mode), eachrow(Final_table_unique))
RPLC_data = Final_table_unique[indices_RPLC,:]

#Removing entries with either a missing SMILES or a missing InChI
RPLC_data = remove_missing_identifiers(RPLC_data)

#Removing chemical entries that are beyond the interquartile range of their set of chemicals with the same identifier
filtered_RPLC_data = remove_outliers_IQR(RPLC_data)

using StatsBase

# Step 1: Assign labels
classes = assign_labels(filtered_RPLC_data.Modifier ./ 100)

# Step 2: Create dictionary to store indices per SMILES
compound_to_indices = Dict{String, Vector{Int}}()
for (i, smiles) in pairs(filtered_RPLC_data.SMILES)
    push!(get!(compound_to_indices, smiles, Int[]), i)
end

# Step 3: Get unique SMILES list
unique_SMILES = unique(filtered_RPLC_data.SMILES)

# Step 4: Initialize dictionary to store majority class per SMILES
compound_majorities = Dict{String, String}()

Random.seed!(42)
# Step 5: For each unique SMILES, find majority class from 'classes'
for (i, smiles) in enumerate(unique_SMILES)
    if i % 100 == 0
        println("Processed $(round(i / length(unique_SMILES) * 100, digits=2))%")
    end

    indices = compound_to_indices[smiles]
    cm = countmap(classes[indices])
    maxval = maximum(values(cm))
    candidates = [k for (k, v) in cm if v == maxval]
    majority_class = length(candidates) == 1 ? candidates[1] : rand(candidates)
    compound_majorities[smiles] = majority_class
end

# Step 6: Create vector of majority classes aligned with all rows
majority_classes = [compound_majorities[smiles] for smiles in filtered_RPLC_data.SMILES]

# Step 7: Copy dataframe and assign new 'Class' column with majority classes
RPLC_data_majority = copy(filtered_RPLC_data)
RPLC_data_majority.Class = majority_classes

RPLC_data_majority.SMILES
RPLC_data_majority.Class

#Obtaining the list of CIDs for the PubChem query for MW and logP
CID_list = DataFrame((CID = unique(collect(skipmissing(RPLC_data_majority.CID)))))

#Get MW and XlogP from PubChem CIDs
CIDs = CID_list[:,1]
index_error = []
# Function to convert strings to integers
for i in eachindex(CIDs)
    try
        if typeof(CIDs[i]) == String || typeof(CIDs[i]) == String15
            CIDs[i] = parse(Int, CIDs[i])
        elseif CIDs[i] == 0
            push!(index_error, i)
        end
    catch
        push!(index_error, i)
    end
end
index_error
CIDs = Int.([CIDs[i] for i in eachindex(CIDs) if i ∉ index_error])
#Getting the MW and XlogP for each chemical from PubCHem
properties = CSV.File(get_for_cids(CIDs; properties="MolecularWeight,XlogP", output="CSV")) |> DataFrame

#Now assigning the unique ones to all the entries
# Create a dictionary mapping for MW and XlogP
MW_dict = Dict(row.CID => row.MolecularWeight for row in eachrow(properties))
XlogP_dict = Dict(row.CID => row.XLogP for row in eachrow(properties))

# Create vectors for MW and XlogP corresponding to the original CIDs vector
MW_vector = [haskey(MW_dict, cid) ? MW_dict[cid] : missing for cid in filtered_RPLC_data[:,5]]
XlogP_vector = [haskey(XlogP_dict, cid) ? XlogP_dict[cid] : missing for cid in filtered_RPLC_data[:,5]]

df_MW_logp = DataFrame(MW = MW_vector, XlogP = XlogP_vector)

#Place the MW and XlogP back with the rest of the data as new columns
RPLC_data_majority = hcat(RPLC_data_majority, df_MW_logp)


RPLC_data_majority = RPLC_data_majority[.!ismissing.(RPLC_data_majority.pH), :]
RPLC_data_majority_pH3 = RPLC_data_majority[in.(RPLC_data_majority.pH, Ref([3.0, 2.6])), :]

sum(RPLC_data_majority_pH3.study .==186)
CSV.write("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_expanded_RPLC_pH3_data.csv", RPLC_data_majority_pH3)

#Load in the fingerprints calculated for the unique SMILES
RepoRT_SMILES_to_FP = CSV.read("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Finerprints\\All RepoRT fingerprints.csv", DataFrame)
missing_FPs = CSV.read("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\missing fingerprints FPS updated data.csv", DataFrame)
missing_FPs_2 = CSV.read("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\missing fingerprints FPS updated data 59.csv", DataFrame)
RepoRT_SMILES_to_FP = append!(RepoRT_SMILES_to_FP, missing_FPs)
RepoRT_SMILES_to_FP = append!(RepoRT_SMILES_to_FP, missing_FPs_2)
# Extract column names
colnames = names(RepoRT_SMILES_to_FP)

# Separate the SMILES column from the fingerprint columns
smiles_col = colnames[1]
fingerprint_cols = colnames[2:end]

# Extract the numerical part of the fingerprint column names and sort them since the current order is: FP0, FP1, FP10, FP100, etc.
sorted_fingerprint_cols = sort(fingerprint_cols, by = x -> parse(Int, match(r"\d+", x).match))

# Combine the SMILES column with the sorted fingerprint columns
sorted_colnames = vcat([smiles_col], sorted_fingerprint_cols)

# Reorder the DataFrame
RepoRT_FP_DF = select(RepoRT_SMILES_to_FP, sorted_colnames)
#allowed_SMILES = RepoRT_FP_DF.x1
#RPLC_data_pH3 = filter(:SMILES => x -> x .∈ Ref(allowed_SMILES), RPLC_data_pH3)
#Create a dictionary of smiles to fingerprints
smiles_to_fingerprints = Dict(RepoRT_FP_DF[i, 1] => RepoRT_FP_DF[i, 2:end] for i in 1:nrow(RepoRT_FP_DF))
#Provide fingerprints for each entry of RepoRT
missing_keys = DataFrame(SMILES = unique(filter(k -> !haskey(smiles_to_fingerprints, k), RPLC_data_majority_pH3.SMILES)))

new_fingerprints = [smiles_to_fingerprints[smiles] for smiles in RPLC_data_majority_pH3.SMILES]


Fingerprints_RepoRT = Int.(zeros(length(RPLC_data_majority_pH3.SMILES), 881))

for i in 1:length(new_fingerprints)
    @show i
    temp = reduce(vcat, new_fingerprints[i])
    Fingerprints_RepoRT[i,:] = temp
end

# Convert the matrix to a DataFrame
df = DataFrame(Fingerprints_RepoRT, sorted_colnames[2:end])

CSV.write("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_expanded_entries_FPS.csv", df)
