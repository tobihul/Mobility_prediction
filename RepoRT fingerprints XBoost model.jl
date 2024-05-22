using CSV, Statistics, DataFrames, StatsPlots, LinearAlgebra

folder_path = "C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data"
Final_table = DataFrame(Inchikey = String[], LC_mode = String7[], MACCS = String[], Pubchem_fps = String[], XlogP = Float64[], Retention_factors = Float64[])
MACCS_keys = readlines( "C:\\Users\\uqthulle\\Documents\\MACCS keys.txt")
PubChem_keys = readlines( "C:\\Users\\uqthulle\\Documents\\PubChem keys.txt")
All_keys =  [MACCS_keys ;PubChem_keys]
All_keys_y = [MACCS_keys ;PubChem_keys; "y"]
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
     gradient_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "gradient.tsv"],"_"))
     gradient = CSV.read(gradient_path, DataFrame)
     local retention_factors
     try 
         run_time = gradient[end,1]
         retention_factors = RT./run_time
     catch
         retention_factors = repeat(["No gradient info"], length(RT))
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

    #Getting MW and XLogP
    descriptor_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "descriptors_canonical_success.tsv"],"_"))
    data_descriptors = CSV.read(descriptor_path, DataFrame)
    MW = data_descriptors.MW
    XLogP = data_descriptors.XLogP
  
    #Getting all the data in one table to append all the other datasets in this format
    current_table = DataFrame(Inchikey = Compound_Inchi, LC_mode = LC_mode, MACCS = MACCS, Pubchem_fps = Pubchem_fps, XlogP = XLogP, Retention_factors = retention_factors)
    Final_table = vcat(Final_table, current_table)
end

# Separate entries with RP and HILIC modes
rp_Inchikey = Final_table[Final_table.LC_mode .== "RP", :].Inchikey
hilic_Inchikey = Final_table[Final_table.LC_mode .== "HILIC", :].Inchikey

# Find the intersection of CIDss for both modes
common_Inchikey = intersect(rp_Inchikey, hilic_Inchikey)

#Removing the ones that have been separated by both as we already have them
condition = in.(Final_table.Inchikey, Ref(common_Inchikey))
Final_table_unique = Final_table[.!condition, :]

valuess = Final_table_unique.Retention_factors
indexes = Int[]
for i = 1:length(valuess)
    if !ismissing(valuess[i]) && typeof(valuess[i]) == Float64
        push!(indexes, i)
    end
end
indexes
Final_table_unique = Final_table_unique[indexes,:]

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




#############################################################
#Data preparation
using Random, MLJ
Random.seed!(42)
X = [MACCS Pubchem_fps]

#Selecting either HILIC or RPLC data
y = Final_table_unique[indices_RPLC,6]
X = X[indices_RPLC,:]

#Shuffling the data
shuffle_indices = shuffle(collect(1:length(y)))
X = X[shuffle_indices,:]
y= y[shuffle_indices]

#Keeping a test set for after parameter optimization
train,test = partition(eachindex(y), 0.9, rng=42)
X_train,X_test = X[train,:],X[test,:]
y_train, y_test = y[train], y[test]

#Preparing the model and model parameter grid search
using XGBoost, OrderedCollections, MLDataUtils
df_results = DataFrame(best_iter = Int[],eval_metric = String[],esr=Int[],max_leaves= Int[],max_depth = Int[],nu=Float32[],score_train=Float32[],score_valid=Float32[],score_test = Float32[])
num_round = 10000
eval_metric = "rmse"
esr_m = [5]
max_depth_m=[3,5,10, 18, 80]
nu_m = [0.3]
max_leaves_m = [5,10,100,1000]
set_data = "RPLC Kfold"
df_results = CSV.read("C:\\Users\\uqthulle\\Documents\\$set_data.csv", DataFrame)
k = 5
kf = MLDataUtils.kfolds(y_train,k)

@time for esr in esr_m
    for max_depth in max_depth_m
        for nu in nu_m
            for max_leaves in max_leaves_m
                    scores_valid = zeros(k)
                    scores_train = zeros(k)
                    scores_test = zeros(k)
                    for fold = collect(1:k)
                        train_fold = kf.train_indices[fold]
                        valid_fold = kf.val_indices[fold]

                        dtrain = DMatrix((X_train[train_fold,:], y_train[train_fold]))
                        dvalid = DMatrix((X_train[valid_fold,:], y_train[valid_fold]))


                        bst = xgboost(dtrain, num_round = num_round,max_leaves=max_leaves, eval_metric = eval_metric, 
                        watchlist = OrderedDict(["train" => dtrain, "eval" => dvalid]), early_stopping_rounds = esr, max_depth=max_depth, η=nu)

                        y_hat_train = XGBoost.predict(bst, X_train[train_fold,:])
                        y_hat_valid = XGBoost.predict(bst, X_train[valid_fold,:])
                        y_hat_test = XGBoost.predict(bst, X_test)

                        scores_train[fold] = cor(y_train[train_fold], y_hat_train)^2
                        scores_valid[fold] = cor(y_train[valid_fold], y_hat_valid)^2
                        scores_test[fold] = cor(y_test, y_hat_test)^2



                        
                    end
                    cv_score_train = mean(scores_train)
                    cv_score_valid = mean(scores_valid)
                    cv_score_test = mean(scores_test)
                    push!(df_results,[bst.best_iteration,eval_metric,esr,max_leaves,max_depth,nu,cv_score_train,cv_score_valid,cv_score_test])
            
            end
        end
    end
end
CSV.write("C:\\Users\\uqthulle\\Documents\\$set_data.csv ", df_results)

vscodedisplay(df_results)

###########
#Best model
bm_index = argmax(df_results.score_valid)

train,valid, test = partition(eachindex(y), 0.7,0.2, rng=42)

dtrain = DMatrix((X[train,:], y[train]))
dvalid = DMatrix((X[valid,:], y[valid]))

bst = xgboost(dtrain, num_round = 10000,max_leaves=df_results.max_leaves[bm_index], eval_metric = eval_metric, watchlist = OrderedDict(["train" => dtrain, "eval" => dvalid]), early_stopping_rounds = df_results.esr[bm_index], max_depth=df_results.max_depth[bm_index], η=df_results.nu[bm_index])

y_hat_train = XGBoost.predict(bst, X[train,:])
y_hat_valid = XGBoost.predict(bst, X[valid,:])
y_hat_test = XGBoost.predict(bst, X[test,:])
score_train = cor(y[train], y_hat_train)^2
score_valid = cor(y[valid], y_hat_valid)^2
score_test = cor(y[test], y_hat_test)^2

scatter(y[train], y_hat_train, xlabel = "Test", ylabel = "Predicted", label = "train = $(round(score_train, digits = 2))", dpi = 300, ylims = (0,1))
scatter!(y[valid], y_hat_valid, label = "valid = $(round(score_test, digits = 2))",
title = "XGBoost model fingerprints to retention factors RPLC")
scatter!(y[test], y_hat_test, label = "test = $(round(score_test, digits = 2))",
title = "XGBoost model fingerprints to retention factors RPLC")


savefig("C:\\Users\\uqthulle\\Documents\\Plots\\RF model XGBoost best parameters.png")


    train, valid, test = partition(eachindex(y), 0.7, 0.2, rng=42)
    dtrain = DMatrix((X[train,:], y[train]))
    dvalid = DMatrix((X[valid,:], y[valid]))
    dm = DMatrix(X, feature_names=All_keys)
    XGBoost.setfeaturenames!(dm, All_keys)
    bst = xgboost(dtrain, num_round = 10000,max_leaves=df_results.max_leaves[bm_index], eval_metric = eval_metric, watchlist = OrderedDict(["train" => dtrain, "eval" => dvalid]), early_stopping_rounds = df_results.esr[bm_index], max_depth=df_results.max_depth[bm_index], η=df_results.nu[bm_index])

    y_hat_train = XGBoost.predict(bst, X[train,:])
    y_hat_valid = XGBoost.predict(bst, X[valid,:])
    y_hat_test = XGBoost.predict(bst, X[test,:])
    score_train = cor(y[train], y_hat_train)^2
    score_valid = cor(y[valid], y_hat_valid)^2
    score_test = cor(y[test], y_hat_test)^2


# can accept tabular data, will keep feature names
df = DataFrame(Float64.([X[train,:] y[train]]), All_keys_y)
bst = xgboost((df[!, All_keys], df.y), num_round = 10000,max_leaves=df_results.max_leaves[bm_index], eval_metric = eval_metric, watchlist = OrderedDict(["train" => dtrain, "eval" => dvalid]), early_stopping_rounds = df_results.esr[bm_index], max_depth=df_results.max_depth[bm_index], η=df_results.nu[bm_index])

imp = importance(bst)

# Initialize empty arrays to store keys and values
keys_arr = String[]
values_arr = Vector{Float32}[]

# Iterate over key-value pairs of the OrderedDict
for (key, value) in imp
    push!(keys_arr, key)
    push!(values_arr, value)
end
values_flat = vcat(values_arr...)
bar_plot = bar(keys_arr[1:15], values_flat[1:15], xtickfont = font(7), 
xrotation=40, dpi = 300, title = "RPLC important variables", bottom_margin = 5Plots.mm)

XGBoost.importancereport(bst)
# Show the plot
show(bar_plot)

sortperm



findfirst(x-> x=="QHAAQH", All_keys)
MACCS
compounds_of_interest = []
for i = eachindex(X[:,54]) 
    if X[i,54] == 1
        push!(compounds_of_interest,i)
    end
end
compounds_of_interest

hello_there = Final_table_unique[indices_HILIC,:][compounds_of_interest,:]

vscodedisplay(hello_there)

println(hello_there[2,1])