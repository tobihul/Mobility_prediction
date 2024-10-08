using PyCall, CSV, DataFrames
# Import RDKit from Python
pd = pyimport("padelpy")
rdk = pyimport("rdkit.Chem")
Fingerprinter = pyimport("rdkit.Chem.EState.Fingerprinter")


smiles_of_interest = CSV.read("R:\\PHD2024TH-Q6813\\Models and other documents\\PubChem_REACH_list.csv", DataFrame).canonicalsmiles
df_smiles = DataFrame(smiles = smiles_of_interest)
@time pubchem_fp = DataFrame(pd.from_smiles(SMILES, fingerprints=true, descriptors=false, timeout = 600, maxruntime = 600))
pubchem_fp = reduce(vcat, pubchem_fp)
missing_Res = hcat(smiles_of_interest, pubchem_fp)

missing_temp = parse.(Int,missing_Res[:,2:end])
missing_fps = hcat(missing_Res[:,1], missing_temp)

PubChem_df
CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\PubChem_fps_part_1.csv", PubChem_df)
MACCS_df
PubChem_df

# Function to extract the numeric part of the column names
function extract_numeric(col_name::String)
    return parse(Int, match(r"\d+", col_name).match)
end
extract_numeric(PubChem_df)
# Sort the column names based on the numeric part
sorted_columns = sort(names(PubChem_df), by=extract_numeric)

# Reorder the DataFrame columns
df = PubChem_df[:, Symbol.(sorted_columns)]

# Combine the results
total_df = hcat(MACCS_df, df)

CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\REACH fingerprints.csv", total_df)


###########For REACH
REACH_data = CSV.read("R:\\PHD2024TH-Q6813\\Models and other documents\\REACH2017_fixedSMILES.csv", DataFrame)

REACH_smiles = sum(ismissing.(REACH_data.SMILES))

unique(REACH_data[:,12])


@time PubChem_df = smiles_to_FPS(REACH_smiles)

# Function to extract the numeric part of the column names
function extract_numeric(col_name::String)
    return parse(Int, match(r"\d+", col_name).match)
end
extract_numeric(PubChem_df)
# Sort the column names based on the numeric part
sorted_columns = sort(names(PubChem_df), by=extract_numeric)

# Reorder the DataFrame columns
df = PubChem_df[:, Symbol.(sorted_columns)]

# Combine the results
total_df = hcat(MACCS_df, df)


