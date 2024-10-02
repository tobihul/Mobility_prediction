include("All functions.jl")
using YAML
#Import the raw data from the RepoRT dataset
folder_path = "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\raw_data"

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

for i = 1:374
    @show i
    current_file_raw = joinpath("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\raw_data", readdir(folder_path)[i+1])
    
    info_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i+1])", "info.tsv"],"_"))
    info = CSV.read(info_path, DataFrame)
    LC_mode = info[1,3]
    if ismissing(LC_mode) || LC_mode == "NA"
        LC_mode = "No mode info"
    end
    
    meta_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i+1])", "metadata.yaml"],"_"))

    try 
        Column_data = YAML.load_file(meta_path)["column"]

    catch
        Column_data = Dict("column" => missing)
    end
    Df_column = DataFrame(Column_data)
    
    pH = 
    
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


##############################
Final_table = DataFrame(SMILES = String[], InChi = String[], LC_mode = String7[], Column_name = String[], pH_A = Int[], pH_B = Int[], CID = String[], Retention_factors = Float64[], Modifier = Float64[])

#Getting the pH
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
    pH_A = meta[!, "eluent.A.pH"]
    if ismissing(pH_A) || pH_A == "NA"

        pH_A = repeat(["No mode info"], length(Compound_SMILES))

    else

        pH_A = repeat(pH_A, length(Compound_SMILES))

    end
    
    pH_B = meta[!, "eluent.A.pH"]
    if ismissing(pH_B) || pH_B == "NA"

        pH_B = repeat(["No mode info"], length(Compound_SMILES))

    else

        pH_B = repeat(pH_B, length(Compound_SMILES))

    end
    #Getting all the data in one table to append all the other datasets in this format
    current_table = DataFrame(SMILES = Compound_SMILES, InChi = InChI, LC_mode = LC_mode, Column_name = Column_name, pH_A = pH_A, pH_B = pH_B, CID = CID, Retention_factors = retention_factors, Modifier = Modifier)
    Final_table = vcat(Final_table, current_table)
end
#Going through each file of the raw data and extracting the necessary information


#Removing some odd entries and entries with missing retention times
Final_table_unique = remove_missing_and_outliers(Final_table)

#Getting indices for RPLC entries
indices_RPLC = findall(row -> occursin("RP", row.LC_mode), eachrow(Final_table_unique))
RPLC_data = Final_table_unique[indices_RPLC,:]

#Removing entries with either a missing SMILES or a missing InChI
RPLC_data = remove_missing_identifiers(RPLC_data)

#Removing chemical entries that are beyond the interquartile range of their set of chemicals with the same identifier
filtered_RPLC_data = remove_outliers_IQR(RPLC_data)



histogram(filtered_RPLC_data.pH_A, bins = 50, formatter = :plain,
dpi = 300, xlabel = "Mobile phase pH", ylabel = "Frequency", legend = false, xlims = (1,14),
xticks = (1:12),xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), grid = :y, ylims = (0,140000))


filtered_RPLC_data.pH_A == filtered_RPLC_data.pH_B

sum(filtered_RPLC_data.pH_B .== 3.0)

sum(filtered_RPLC_data.pH_B .== 2.0)

unique_pH = sort(unique(filtered_RPLC_data.pH_B))

counts = []
for pH in unique_pH
    push!(counts, (sum(filtered_RPLC_data.pH_B .== pH)))
end
counts
percentages = round.(counts./sum(counts) *100, digits = 1)

pH_counts = DataFrame(pH = unique_pH, entries = counts, Percentage = percentages)

CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\pH table.csv", pH_counts)

smiles_list = string.(filtered_RPLC_data[:,1])

open("molecules.smi", "w") do file
    for smiles in smiles_list
        write(file, smiles * "\n")
    end
end