include("All functions.jl")
folder_path = "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\raw_data"
Final_table = DataFrame(SMILES = String[], InChi = String[], LC_mode = String7[], Column_name = String[], CID = String[], Retention_factors = Float64[], Modifier = Float64[])
for i = 1:374
    #Locating the file
    current_file_raw = joinpath("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\raw_data", readdir(folder_path)[i+1])
    current_file_processed = joinpath("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\processed_data", readdir(processed_folder_path)[1:end-1][i])
    println(i)
    #Getting the SMILES and retention times
    rt_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i+1])", "rtdata.tsv"],"_"))
    data_compounds = CSV.read(rt_path, DataFrame)
    Compound_SMILES = data_compounds[:,7]
    InChI = data_compounds[:,8]
    RT = data_compounds[:,4]
    CID = data_compounds[:,5]
    #Getting the retention factors
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

    #Metadata 
    meta_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i+1])", "metadata.tsv"],"_"))
    meta = CSV.read(meta_path, DataFrame)
    Column_name = meta[1,2]
    #Getting all the data in one table to append all the other datasets in this format
    current_table = DataFrame(SMILES = Compound_SMILES, InChi = InChI, LC_mode = LC_mode, Column_name = Column_name, CID = CID, Retention_factors = retention_factors, Modifier = Modifier)
    Final_table = vcat(Final_table, current_table)
end

Final_table

Final_table_unique = remove_missing_and_outliers(Final_table)


#Getting indices for RPLC and HILIC rows
indices_RPLC = findall(row -> occursin("RP", row.LC_mode), eachrow(Final_table_unique))

RPLC_data = Final_table_unique[indices_RPLC,:]

RPLC_data = remove_missing_identifiers(RPLC_data)

filtered_RPLC_data = remove_outliers_IQR(RPLC_data)

unique_columns = DataFrame(Column_name = unique(filtered_RPLC_data.Column_name))

CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\List of RPLC columns used.csv", unique_columns)

unique_SMILES = DataFrame(SMILES = string.(unique(filtered_RPLC_data[:,1])))

CID_list = DataFrame((CID = unique(collect(skipmissing(filtered_RPLC_data.CID)))))

CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\Unique smiles all chemicals RepoRT.csv", unique_SMILES)
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

properties = CSV.File(get_for_cids(CIDs; properties="MolecularWeight,XlogP", output="CSV")) |> DataFrame
properties

# Create a dictionary mapping for MW and XlogP
MW_dict = Dict(row.CID => row.MolecularWeight for row in eachrow(properties))
XlogP_dict = Dict(row.CID => row.XLogP for row in eachrow(properties))

# Create vectors for MW and XlogP corresponding to the original CIDs vector
MW_vector = [haskey(MW_dict, cid) ? MW_dict[cid] : missing for cid in filtered_RPLC_data[:,4]]
XlogP_vector = [haskey(XlogP_dict, cid) ? XlogP_dict[cid] : missing for cid in filtered_RPLC_data[:,4]]

df_MW_logp = DataFrame(MW = MW_vector, XlogP = XlogP_vector)

filtered_RPLC_data = hcat(filtered_RPLC_data, df_MW_logp)

hist_MW = histogram(filtered_RPLC_data.MW, 
bins = 500, 
linecolor=:transparent,
 xlims = (0,1000),
dpi = 300,
xlabel = ("Molecular Weight (Da)"),
label = false,
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11))

hist_XLOGP = histogram(filtered_RPLC_data.XlogP,
bins = 100,  
linecolor=:transparent,
dpi = 300,
xlabel = ("XlogP"),
label = false, formatter = :plain,
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11))

hist_RF = histogram(filtered_RPLC_data.Retention_factors,  
linecolor=:transparent,
dpi = 300,
xlabel = ("Retention factors"),
label = false, formatter = :plain,
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11))

hist_MOD = histogram(filtered_RPLC_data.Modifier./100,
xlims = (0,1),
linecolor=:transparent,
dpi = 300,
xlabel = ("Φ"),
label = false, formatter = :plain,
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11))

cd("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\Plots")
savefig("Distribution PHI entire dataset.png")
cd("R:\\PHD2024TH-Q6813\\Code\\Regression")

scatter(filtered_RPLC_data.XlogP, filtered_RPLC_data.Modifier)


# Example DataFrame
df = DataFrame(SMILES=["CCO", "CCN", "CCO"], OtherData=[1, 2, 3])

# Unique SMILES and their fingerprints
unique_smiles = ["CCO", "CCN"]
fingerprints = [
    [0, 1, 0, 1, 1, 1],  # Fingerprint for CCO (shortened for illustration)
    [1, 0, 1, 0, 0, 0]   # Fingerprint for CCN (shortened for illustration)
]

# Create a dictionary mapping SMILES to fingerprints
smiles_to_fingerprint = Dict{String, Vector{Int}}()

for (smiles, fingerprint) in zip(unique_smiles, fingerprints)
    smiles_to_fingerprint[smiles] = fingerprint
end

# Add fingerprints to the DataFrame
df[:, :Fingerprint] = [smiles_to_fingerprint[smiles] for smiles in df.SMILES]

# Print the resulting DataFrame
println(df)