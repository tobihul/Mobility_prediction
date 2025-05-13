include("All functions.jl")
#Check the file All funtions.jl for all the specific functions written for this file

#Import the raw data from the RepoRT dataset
folder_path = "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\raw_data"

Final_table = DataFrame(SMILES = String[], InChi = String[], LC_mode = String7[], Column_name = String[], CID = String[], Retention_factors = Float64[], Modifier = Float64[])

#These are the fingerprint names for each corresponding column
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")

#Going through each file of the raw data and extracting the necessary information
for i = 1:374
    #Locating the current file
    current_file_raw = joinpath("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\raw_data", readdir(folder_path)[i+1])
    println(i)
    #Getting the SMILES, InChI, retention times and CID for each chemical in the current file
    rt_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i+1])", "rtdata.tsv"],"_"))
    data_compounds = CSV.read(rt_path, DataFrame)
    Compound_SMILES = data_compounds[:,7]
    InChI = data_compounds[:,8]
    RT = data_compounds[:,4]
    CID = data_compounds[:,5]

    #Getting the retention factors and fraction of organic mobile phase at elution if available
    gradient_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i+1])", "gradient.tsv"],"_"))
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
    info_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i+1])", "info.tsv"],"_"))
    info = CSV.read(info_path, DataFrame)
    #Check if the column info is missing in the meta data
    LC_mode = info[1,3]
    if ismissing(LC_mode) || LC_mode == "NA"

        LC_mode = repeat(["No mode info"], length(RT))

    else

        LC_mode = repeat([LC_mode], length(RT))

    end

    #Getting the column name for each method
    meta_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i+1])", "metadata.tsv"],"_"))
    meta = CSV.read(meta_path, DataFrame)
    Column_name = meta[1,2]
    #Getting all the data in one table to append all the other datasets in this format
    current_table = DataFrame(SMILES = Compound_SMILES, InChi = InChI, LC_mode = LC_mode, Column_name = Column_name, CID = CID, Retention_factors = retention_factors, Modifier = Modifier)
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

#Obtaining the list of columns used for the separation of the chemicals
unique_columns = DataFrame(Column_name = unique(filtered_RPLC_data.Column_name))
CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\List of RPLC columns used.csv", unique_columns)

#Obtaining the list of unique smiles for fingerprint calculation
unique_SMILES = DataFrame(SMILES = string.(unique(filtered_RPLC_data[:,1])))
CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\Unique smiles all chemicals RepoRT.csv", unique_SMILES)

#Obtaining the list of CIDs for the PubChem query for MW and logP
CID_list = DataFrame((CID = unique(collect(skipmissing(filtered_RPLC_data.CID)))))

CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\Unique CID all chemicals RepoRT.csv", CID_list)

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
filtered_RPLC_data = hcat(filtered_RPLC_data, df_MW_logp)
CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\filtered_RPLC_data.csv", filtered_RPLC_data)


####This next part is after calculating all the fingerprints using Fingerprint_generation.jl

#Load in the fingerprints calculated for the unique SMILES
RepoRT_SMILES_to_FP = CSV.read("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Finerprints\\All RepoRT fingerprints.csv", DataFrame)

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

#Create a dictionary of smiles to fingerprints
smiles_to_fingerprints = Dict(RepoRT_FP_DF[i, 1] => RepoRT_FP_DF[i, 2:end] for i in 1:nrow(RepoRT_FP_DF))

#Provide fingerprints for each entry of RepoRT
new_fingerprints = [smiles_to_fingerprints[smiles] for smiles in filtered_RPLC_data[:,1]]

Fingerprints_RepoRT = Int.(zeros(length(RPLC_with_REACH[:,1]), 881))

for i in 1:length(new_fingerprints)
    @show i
    temp = reduce(vcat, new_fingerprints[i])
    Fingerprints_RepoRT[i,:] = temp
end

# Convert the matrix to a DataFrame
df = DataFrame(Fingerprints_RepoRT, sorted_colnames[2:end])

# Save the DataFrame as a CSV file
CSV.write("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT with REACH overlap.csv", df)

##We also need to sort the REACH fingerprints correctly so we can do the same again but for REACH
REACH_data = CSV.read("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Finerprints\\REACH SMILES and FPS.csv", DataFrame)

# Extract column names
colnames = names(REACH_data)

# Separate the SMILES column from the fingerprint columns
properties_cols = colnames[1:3]
fingerprint_cols = colnames[4:end]

# Extract the numerical part of the fingerprint column names and sort them since the current order is: FP0, FP1, FP10, FP100, etc.
sorted_fingerprint_cols = sort(fingerprint_cols, by = x -> parse(Int, match(r"\d+", x).match))

# Combine the SMILES column with the sorted fingerprint columns
sorted_colnames = vcat(properties_cols, sorted_fingerprint_cols)

# Reorder the DataFrame
REACH_ordered_df = select(REACH_data, sorted_colnames)

##Save the ordered DataFrame
CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\REACH fingerprints final.csv", REACH_ordered_df)

RPLC_pH3_pH26_data = CSV.read("C:\\Users\\uqthulle\\Documents\\RPLC_data_pH3and2.6_updated.csv", DataFrame)
ARP_2022_data = CSV.read("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\ARP 2022 paper exp kocs.csv", DataFrame)

ARP_2022_data[:,2]
ARP_2022_data = dropmissing!(ARP_2022_data)
filtered_ARP = filter(row -> occursin("exp", row.Exp_log_koc), ARP_2022_data)
filtered_ARP.Exp_log_koc = parse.(Float64, replace.(filtered_ARP.Exp_log_koc, r"[^\d\.]+" => ""))
indices_unique = findfirst(x-> x == (unique(common_RPLC.SMILES)), RPLC_pH3_pH26_data.SMILES)
indices_unique
RPLC_pH3_pH26_data[indices_unique,:]
indices = findall(x -> x in filtered_ARP.SMILES, RPLC_pH3_pH26_data.SMILES)


common_RPLC = RPLC_pH3_pH26_data[indices, [1,7]]

common_ARP = filtered_ARP[indices,:]

combined = innerjoin(RPLC_pH3_pH26_data, filtered_ARP, on = :SMILES)

first(combined, 5)

correlation = cor(combined.Modifier, combined.Exp_log_koc)

scatter(combined.Modifier, combined.Exp_log_koc, 
xlabel = "% Modifier", ylabel = "exp log Koc"
, dpi = 300, label = "R² = $(correlation)")

using GLM

lm_model = lm(@formula(Exp_log_koc ~ Modifier), combined)

# Generate predicted values
x_vals = sort(combined.Modifier)
y_vals = GLM.predict(lm_model, DataFrame(Modifier = x_vals))

plot!(x_vals, y_vals)