using PyCall, CSV, DataFrames
# Import RDKit from Python
pd = pyimport("padelpy")
rdk = pyimport("rdkit.Chem")
Fingerprinter = pyimport("rdkit.Chem.EState.Fingerprinter")


smiles_of_interest = CSV.read("R:\\PHD2024TH-Q6813\\Models and other documents\\PubChem_REACH_list.csv", DataFrame).canonicalsmiles
smiles_of_interest = CSV.read("C:\\Users\\uqthulle\\Documents\\missing fingerprints updated data.csv", DataFrame)
smiles_of_interest = CSV.read("C:\\Users\\uqthulle\\Documents\\missing fingerprints 59 SMILES.csv", DataFrame)

SMILES = smiles_of_interest.SMILES
SMILES[1]
rows = DataFrame()
errored_smiles = String[]

really_failed =[8,18]
nr_failed = 20
pubchem_fps = pd.from_smiles(errored_smiles[nr_failed], fingerprints = true, descriptors = false)
failed = DataFrame(Dict(Symbol(k) => parse(Int, v) for (k, v) in pubchem_fps))
df_failed = insertcols!(failed, 1, :SMILES => errored_smiles[nr_failed])

append!(rows, df_failed)


batch_size = 10
n = length(SMILES)

for i in 1:batch_size:n
    @show i
    batch_end = min(i + batch_size - 1, n)
    batch = SMILES[i:batch_end]

    try
        # Get vector of Dicts (one per SMILES)
        pubchem_fps = pd.from_smiles(batch, fingerprints = true, descriptors = false)

        for (j, fp_dict) in enumerate(pubchem_fps)
            try
                row = DataFrame(Dict(Symbol(k) => parse(Int, v) for (k, v) in fp_dict))
                insertcols!(row, 1, :SMILES => batch[j])
                append!(rows, row)
            catch e_inner
                @warn "Error converting fingerprint for SMILES: $(batch[j])" exception = e_inner
                push!(errored_smiles, batch[j])
            end
        end

    catch e_outer
        @warn "Entire batch failed from index $i to $batch_end" exception = e_outer
        append!(errored_smiles, batch)
    end
end

PubChem_df = rows[:,2:end]
SMILES = rows[:,1]
rows
#@time pubchem_fp = DataFrame(pd.from_smiles(SMILES, fingerprints=true, descriptors=false, timeout = 600, maxruntime = 600))
#pubchem_fp = reduce(vcat, pubchem_fp)
#missing_Res = hcat(smiles_of_interest[1:10], pubchem_fp)

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
insertcols!(df, 1, :x1 => SMILES)

# Combine the results
total_df = hcat(MACCS_df, df)

CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\REACH fingerprints.csv", total_df)
CSV.write("C:\\Users\\uqthulle\\Documents\\missing fingerprints FPS updated data 59.csv", df)


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


