include("All functions.jl")

#Import the updated raw data
folder_path = "R:\\PHD2024TH-Q6813\\Research files\\RepoRT updated 08042025\\raw_data"

Final_table = DataFrame(SMILES = String[], InChi = String[], LC_mode = String7[], Column_name = String[], CID = String[], Retention_factors = Float64[], Modifier = Float64[], pH = Float64[])

#These are the fingerprint names for each corresponding column
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")
#Going through each file of the raw data and extracting the necessary information
for i = 1:length(readdir(folder_path))
    #Locating the current file
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
    current_table = DataFrame(SMILES = Compound_SMILES, InChi = InChI, LC_mode = LC_mode, Column_name = Column_name, CID = CID, Retention_factors = retention_factors, Modifier = Modifier, pH = pHA)
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

unique_compounds = unique(filtered_RPLC_data.SMILES)
classes = assign_labels(filtered_RPLC_data.Modifier./100)

using StatsBase
compound_majorities = Dict{String, String}()

for (i, compound) in enumerate(unique_compounds)
    if i % 100 == 0
        println("Processed $(round(i / length(unique_compounds) * 100, digits = 2))%")
    end

    chemical_entries = findall(x -> x == compound, filtered_RPLC_data.SMILES)
    majority_class = findmax(countmap(classes[chemical_entries]))[2]
    compound_majorities[compound] = majority_class
end
unique_smiles = collect(keys(compound_majorities))

summary_df = DataFrame(
    SMILES = unique_smiles,
    Class = [compound_majorities[smiles] for smiles in unique_smiles],
)

get_first(field, smiles) = first(filtered_RPLC_data[filtered_RPLC_data.SMILES .== smiles, field])
get_mean(field, smiles) = mean(filtered_RPLC_data[filtered_RPLC_data.SMILES .== smiles, field])
function get_mode(field::Symbol, smiles::String)
    values = filtered_RPLC_data[filtered_RPLC_data.SMILES .== smiles, field]
    values = skipmissing(values)
    isempty(values) ? missing : mode(values)
end
summary_df.CID = [get_first(:CID, smiles) for smiles in summary_df.SMILES]
summary_df.pH = [get_mode(:pH, smiles) for smiles in summary_df.SMILES]
summary_df.Modifier = [get_mean(:Modifier, smiles) for smiles in summary_df.SMILES]
summary_df.Retention_factors = [get_mean(:Retention_factors, smiles) for smiles in summary_df.SMILES]


index_error = Int[]
for i in eachindex(summary_df.CID)
    try
        if typeof(summary_df.CID[i]) <: AbstractString
            summary_df.CID[i] = parse(Int, summary_df.CID[i])
        elseif summary_df.CID[i] == 0
            push!(index_error, i)
        end
    catch
        push!(index_error, i)
    end
end

# Remove rows with CID errors
summary_df = summary_df[setdiff(1:end, index_error), :]

CIDs = summary_df.CID
properties = CSV.File(get_for_cids(CIDs; properties="MolecularWeight,XlogP", output="CSV")) |> DataFrame

MW_dict = Dict(row.CID => row.MolecularWeight for row in eachrow(properties))
XlogP_dict = Dict(row.CID => row.XLogP for row in eachrow(properties))

summary_df.MW = [get(MW_dict, cid, missing) for cid in summary_df.CID]
summary_df.XlogP = [get(XlogP_dict, cid, missing) for cid in summary_df.CID]

summary_df

summary_df = summary_df[.!ismissing.(summary_df.pH), :]
summary_df_pH3 = summary_df[in.(summary_df.pH, Ref([3.0, 2.6])), :]

CSV.write("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_pH3_RPLC_data.csv", summary_df_pH3)

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
missing_keys = DataFrame(SMILES = unique(filter(k -> !haskey(smiles_to_fingerprints, k), summary_df_pH3.SMILES)))

new_fingerprints = [smiles_to_fingerprints[smiles] for smiles in summary_df_pH3.SMILES]


Fingerprints_RepoRT = Int.(zeros(length(summary_df_pH3.SMILES), 881))

for i in 1:length(new_fingerprints)
    @show i
    temp = reduce(vcat, new_fingerprints[i])
    Fingerprints_RepoRT[i,:] = temp
end

# Convert the matrix to a DataFrame
df = DataFrame(Fingerprints_RepoRT, sorted_colnames[2:end])

CSV.write("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_entries_FPS.csv", df)
