include("All functions.jl")
using YAML
#Import the raw data from the RepoRT dataset
folder_path = "R:\\PHD2024TH-Q6813\\Research files\\RepoRT updated 08042025\\raw_data"
expected_columns =  ["flowrate"
"id"
"length"
"name"
"particle.size"
"temperature"
"usp.code"]
Column_table = DataFrame()
Mode_table = Vector{String}()

# Define the expected column names

for i = 1:length(readdir(folder_path))
    @show i
    current_file_raw = joinpath("R:\\PHD2024TH-Q6813\\Research files\\RepoRT updated 08042025\\raw_data", readdir(folder_path)[i])
    
    info_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i])", "info.tsv"],"_"))
    info = CSV.read(info_path, DataFrame)
    LC_mode = info[1,3]
    if ismissing(LC_mode) || LC_mode == "NA"
        LC_mode = "No mode info"
    end
    
    meta_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i])", "metadata.yaml"],"_"))

    try 
        Column_data = YAML.load_file(meta_path)["column"]

    catch
        Column_data = Dict("column" => missing)
    end
    Df_column = DataFrame(Column_data)
    
    
    if length(Df_column[1,:]) < 7 
        for col in expected_columns
            if !(col in names(Df_column))
                Df_column[!, Symbol(col)] = fill(missing, nrow(Df_column))
            end
        end
    end
    
    # Select the expected columns and ensure order
    Df_column = Df_column[:, expected_columns]
    Column_table = vcat(Column_table, Df_column)
    push!(Mode_table, LC_mode)
end


indices_RPLC = findall(x-> x == "RP", Mode_table)

RPLC_column_table = Column_table[indices_RPLC,:]

CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\All_column_data.csv", RPLC_column_table)
CSV.write("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\All_column_data_updated", RPLC_column_table)

###ph table

RPLC_pH3_pH26_data = CSV.read("C:\\Users\\uqthulle\\Documents\\RPLC_data_pH3and2.6_updated.csv", DataFrame)

unique_pH = sort(unique(RPLC_pH3_pH26_data.pH))

counts = []
for pH in unique_pH
    push!(counts, (sum(RPLC_pH3_pH26_data.pH .== pH)))
end
counts
percentages = round.(counts./sum(counts) *100, digits = 1)

pH_counts = DataFrame(pH = unique_pH, entries = counts, Percentage = percentages)

smiles_list = string.(RPLC_pH3_pH26_data[:,1])

open("molecules.smi", "w") do file
    for smiles in smiles_list
        write(file, smiles * "\n")
    end
end

