using PyCall, CSV, DataFrames
# Import RDKit from Python
pd = pyimport("padelpy")
rdk = pyimport("rdkit.Chem")
Fingerprinter = pyimport("rdkit.Chem.EState.Fingerprinter")

smiles_of_interest =  CSV.read("R:\\PHD2024TH-Q6813\\Models and other documents\\S36_UBAPMT_April2022.csv", DataFrame).Single_Constituent_SMILES

# Initialize PyCall with Python dependencies
rdk = pyimport("rdkit")
Chem = rdk.Chem
MACCSkeys = rdk.Chem.rdMolDescriptors.GetMACCSKeysFingerprint
pd = pyimport("padelpy")


n = length(smiles_of_interest)
    
    # Create DataFrames for MACCS and PubChem fingerprints
    MACCS_df = DataFrame([Symbol("MACCS_$i") => Vector{Float64}(undef, n) for i in 1:166]) # MACCS has 166 bits
    PubChem_df = DataFrame()
    
    for i in eachindex(smiles_of_interest)
        @show i
            smiles = smiles_of_interest[i]
            mol = Chem.MolFromSmiles(smiles)
            
            # Generate MACCS fingerprints
            maccs_fp = MACCSkeys(mol).GetOnBits()
            maccs_fp_array = [Int(i in maccs_fp) for i in 0:165]
            MACCS_df[i, :] = maccs_fp_array
            
            # Generate PubChem fingerprints using PaDEL-Descriptor
            pubchem_fp = DataFrame(pd.from_smiles(smiles, fingerprints=true, descriptors=false))
            PubChem_df =vcat(PubChem_df, pubchem_fp)
    end
    
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


