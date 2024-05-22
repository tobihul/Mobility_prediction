using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots, MultivariateStats
folder_path = "C:\\Users\\uqthulle\\Documents\\RepoRT-master\\raw_data"
Final_table = DataFrame(canonical_smiles = String[],  LC_mode = String7[], pH_A = Float64[], pH_B = Float64)

for i = 2:376
    #Locating the file
    current_file = joinpath("C:\\Users\\uqthulle\\Documents\\RepoRT-master\\raw_data", readdir(folder_path)[i])
    println(i)
    #Getting the compounds name, canonical SMILES and retention time 
    rt_path = joinpath(current_file, join(["$(readdir(folder_path)[2:end][i-1])", "rtdata.tsv"],"_"))
    data_compounds = CSV.read(rt_path, DataFrame)
    compound_canonical_SMILES = data_compounds[:,7]
    

    metadata_path = joinpath(current_file, join(["$(readdir(folder_path)[2:end][i-1])", "metadata.tsv"],"_"))
    meta_data = CSV.read(metadata_path, DataFrame)
    #Check if the column info is missing in the meta data
    pH_A = ones(length(data_compounds[:,1])) .* meta_data[!, Symbol("eluent.A.pH")][1]
    try 
        pH_B = ones(length(data_compounds[:,1])) .* meta_data[!, Symbol("eluent.B.pH")][1]
        
    catch
        pH_B = ones(length(data_compounds[:,1])) .* 3
    end

    #Getting the LC mode info for the data set
    info_path = joinpath(current_file, join(["$(readdir(folder_path)[2:end][i-1])", "info.tsv"],"_"))
    info = CSV.read(info_path, DataFrame)
    #Check if the column info is missing in the meta data
    LC_mode = info[1,3]
    if ismissing(LC_mode) || LC_mode == "NA"

        LC_mode = repeat(["No mode info"], length(data_compounds[:,1]))

    else

        LC_mode = repeat([LC_mode], length(data_compounds[:,1]))

    end
 
    #Getting all the data in one table to append all the other datasets in this format
    current_table = DataFrame(canonical_smiles = compound_canonical_SMILES, LC_mode = LC_mode, pH_A = pH_A, pH_B = pH_B)
    Final_table = vcat(Final_table, current_table)
end

a = histogram(Final_table.pH_A, grid = false, xticks = (1:1:14), title = "Mobile phase A pH n = $(length(Final_table.pH_A))", 
xlabel = "pH", yticks = false, bottom_margin = 5Plots.mm)
b = histogram(Final_table.pH_B, grid = false, xticks = (1:1:14), title = "Mobile phase A pH n = $(length(Final_table.pH_A))", 
xlabel = "pH", yticks = false, bottom_margin = 5Plots.mm)

plot(a,b, size = (800,400), dpi = 300)

savefig("C:\\Users\\uqthulle\\Documents\\Plots\\Ph mobile phases.png")