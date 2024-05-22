using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots, MultivariateStats

folder_path = "C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data"
Final_table = DataFrame(Inchikey = String[], LC_mode = String7[], MACCS = String[], Pubchem_fps = String[], Retention_factors = Float64[])

for i = 1:374
    #Locating the file
    current_file = joinpath("C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data", readdir(folder_path)[i])
    println(i)
    #Getting the SMILES and retention times
    rt_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "rtdata_canonical_success.tsv"],"_"))
    data_compounds = CSV.read(rt_path, DataFrame)
    Compound_Inchi = data_compounds[:,6]
    RT = data_compounds[:,4]

     #Getting the retention factors
     metadata_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "metadata.tsv"],"_"))
     metadata = CSV.read(metadata_path, DataFrame)
     t_0 = metadata[:,9][1]
     if t_0 != 0
         retention_factors = (RT.-t_0)./t_0
     else
         retention_factors = String.(repeat(["No t0 info"], length(RT)))
     end


    #Getting the MACCS
    MACCS_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "fingerprints_maccs_canonical_success.tsv"],"_"))
    MACCS_info = CSV.read(MACCS_path, DataFrame)
    df_merged = innerjoin(data_compounds, MACCS_info, on = :id)
    MACCS = df_merged[:,end]
    #Getting the PubChemFingerprints
    Pubchem_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "fingerprints_pubchem_canonical_success.tsv"],"_"))
    Pubchem_info = CSV.read(Pubchem_path, DataFrame)
    df_merged_2 = innerjoin(data_compounds, Pubchem_info, on = :id)
    Pubchem_fps = df_merged_2[:,end]

    #Getting the LC mode info for the data set
    info_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "info.tsv"],"_"))
    info = CSV.read(info_path, DataFrame)
    #Check if the column info is missing in the meta data
    LC_mode = info[1,3]
    if ismissing(LC_mode) || LC_mode == "NA"

        LC_mode = repeat(["No mode info"], length(MACCS))

    else

        LC_mode = repeat([LC_mode], length(MACCS))

    end

    
  
    #Getting all the data in one table to append all the other datasets in this format
    current_table = DataFrame(Inchikey = Compound_Inchi, LC_mode = LC_mode, MACCS = MACCS, Pubchem_fps = Pubchem_fps, Retention_factors = retention_factors)
    Final_table = vcat(Final_table, current_table)
end

Final_table

# Separate entries with RP and HILIC modes
rp_Inchikey = Final_table[Final_table.LC_mode .== "RP", :].Inchikey
hilic_Inchikey = Final_table[Final_table.LC_mode .== "HILIC", :].Inchikey

# Find the intersection of CIDss for both modes
common_Inchikey = intersect(rp_Inchikey, hilic_Inchikey)

#Removing the ones that have been separated by both as we already have them
condition = in.(Final_table.Inchikey, Ref(common_Inchikey))
Final_table_unique = Final_table[.!condition, :]




Final_table_unique
values = Final_table_unique.Retention_factors
indexes = Int[]
for i = 1:length(values)
    if typeof(values[i]) == Float64
        push!(indexes, i)
    end
end
indexes
Final_table_unique = Final_table_unique[indexes,:]
println(Final_table_unique[:,5])

MACCS = Int.(zeros(length(Final_table_unique[:,2]),166))
for i = 1:length(Final_table_unique[:,3])
    string_bits = Final_table_unique[i,3]
    vectorized_bits = parse.(Int, split(string_bits, ","))
    for v in vectorized_bits
        MACCS[i,v] = 1
    end
end
MACCS
Pubchem_fps = Int.(zeros(length(Final_table_unique[:,4]),881))
for i = 1:length(Final_table_unique[:,4])
    string_bits = Final_table_unique[i,4]
    vectorized_bits = parse.(Int, split(string_bits, ","))
    for v in vectorized_bits
        Pubchem_fps[i,v] = 1
    end
end
Pubchem_fps

indices_RPLC = findall(row -> occursin("RP", row.LC_mode), eachrow(Final_table_unique))
indices_HILIC = findall(row -> occursin("HILIC", row.LC_mode), eachrow(Final_table_unique))

X = [MACCS Pubchem_fps]

#Random Forest Regression
X = X[indices_HILIC,:]
y = Final_table[indices_HILIC,5]

using ScikitLearn
using ScikitLearn.CrossValidation: cross_val_score
using ScikitLearn: @sk_import
@sk_import ensemble: RandomForestRegressor
@sk_import model_selection: train_test_split
using ScikitLearn.GridSearch: GridSearchCV

n_estimators = collect(100:100:600)
        min_samples_leaf = collect(2:2:8)
        max_features = ["sqrt", "log2"]
        param_grid = Dict("n_estimators" => n_estimators, "min_samples_leaf" => min_samples_leaf, "max_features" => max_features)

# Create a random forest regressor
rf_regressor = RandomForestRegressor()

# Perform grid search with cross-validation
grid_search = GridSearchCV(rf_regressor, param_grid, cv = 3)

# Fit the grid search to the data
fit!(grid_search, X_train, y_train)

# Get the best parameters and best score
best_params = grid_search.best_params_
best_score = grid_search.best_score_

println(y)

