using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots, MultivariateStats
using GLM, TypedTables, LinearAlgebra, ScikitLearn, Random, MLJ, MLDataUtils
using ScikitLearn: @sk_import
@sk_import ensemble: RandomForestRegressor
@sk_import model_selection: StratifiedKFold
folder_path = "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\processed_data"
Final_table = DataFrame(Inchikey = String[], LC_mode = String7[], MACCS = String[], Pubchem_fps = String[],MW = Float64[], XlogP = Float64[], Retention_factors = Float64[], Modifier = Float64[] )
MACCS_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\MACCS keys.txt")
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")
All_keys =  [MACCS_keys ;PubChem_keys]
All_keys_y = [MACCS_keys ;PubChem_keys]
function interpolate_B_modifier(time::Float64, gradient::DataFrame)
    times = gradient[:,1]./run_time
    idx = searchsortedlast(times, time)
    if idx == 0
        return gradient[1,3]
    elseif idx == length(times)
        return gradient[end,3]
    elseif time > 1
        return gradient[end,3]
    else
        t1, t2 = times[idx], times[idx+1]
        B1, B2 = gradient[:,3][idx], gradient[:,3][idx+1]
        return B1 + (B2 - B1) * (time - t1) / (t2 - t1)
    end
end
function calculate_leverage(X,x)
    leverage = zeros(length(x[:,1]))
    b = pinv(X'*X)
    for i = 1:length(x[:,1])
        leverage[i] = (x[i,:])' * b * x[i,:] 
    end
    return leverage
end
function remove_overlapping(Results_table)

    # Separate entries with RP and HILIC modes
    rp_Inchikey = Results_table[Results_table.LC_mode .== "RP", :].Inchikey
    hilic_Inchikey = Results_table[Results_table.LC_mode .== "HILIC", :].Inchikey

    # Find the intersection of CIDss for both modes
    common_Inchikey = intersect(rp_Inchikey, hilic_Inchikey)

    #Removing the ones that have been separated by both as we already have them
    condition = in.(Results_table.Inchikey, Ref(common_Inchikey))
    Final_table_unique = Results_table[.!condition, :]

    return Final_table_unique
end
function remove_missing_and_outliers(Results_unique)
    #This removes rows that have missing retention factors due to no gradient data or other 
    valuess = Results_unique.Retention_factors
    indexes = Int[]
    for i = 1:length(valuess)
        if !ismissing(valuess[i]) && typeof(valuess[i]) == Float64
            push!(indexes, i)
        end
    end
    indexes
    Results_unique = Results_unique[indexes,:]

    #This removes outliers with high retention factor but low %B
    rtfs = Results_unique[:,end-1]
    pmbph = Results_unique[:,end]
    indexes = Int[]
    for i = 1:length(rtfs)
        if rtfs[i] > 0.65 && pmbph[i] < 12.5
            continue
        else
            push!(indexes, i)
        end
    end
    indexes
    Results_unique = Results_unique[indexes,:]
    return Results_unique
end
function format_fingerprints(Final_table_unique)
    #First we get MACCS
    MACCS = Int.(zeros(length(Final_table_unique[:,2]),166))
    for i = 1:length(Final_table_unique[:,3])
        string_bits = Final_table_unique[i,3]
        vectorized_bits = parse.(Int, split(string_bits, ","))
        for v in vectorized_bits
            MACCS[i,v] = 1
        end
    end
    #Now PubChem FPs
    Pubchem_fps = Int.(zeros(length(Final_table_unique[:,4]),881))
    for i = 1:length(Final_table_unique[:,4])
        string_bits = Final_table_unique[i,4]
        vectorized_bits = parse.(Int, split(string_bits, ","))
        for v in vectorized_bits
            Pubchem_fps[i,v] = 1
        end
    end
    return MACCS, Pubchem_fps
end
function custom_binning(y, num_bins)
    min_val = minimum(y)
    max_val = maximum(y)
    bin_width = (max_val - min_val) / num_bins
    bins = Array{String}(undef, length(y))
    
    for i in 1:length(y)
        bin_index = ceil(Int, (y[i] - min_val) / bin_width)
        bins[i] = string( min(bin_index, num_bins))
    end
    
    return bins
end
function Stratifiedcv(X, groups, n_folds)
    # Initialize StratifiedKFold
    skf = StratifiedKFold(n_splits=n_folds)

    # Initialize variables to store indices for each fold
    train_indices = Vector{Vector{Int}}(undef, n_folds)
    test_indices = Vector{Vector{Int}}(undef, n_folds)

    # Iterate over the splits and store indices for each fold
    for (i, (train_idx, test_idx)) in enumerate(skf.split(X, groups))
        train_indices[i] = train_idx .+ 1  # Adjust for 1-based indexing
        test_indices[i] = test_idx .+ 1  # Adjust for 1-based indexing
    end
    return train_indices, test_indices
end
function train_test_split_no_leakage_regressor(filtered_RPLC_data, split)
    fingerprints = PubChem_fps

    unique_compounds = unique(filtered_RPLC_data[:,1])

    index_first_occurrence = Vector{Int}(undef, length(unique_compounds))

    train, test = train_test_split(unique_compounds, test_size = 0.1, random_state = 42)
    
    # Initialize variables
    train_indices = []
    test_indices = []

    # Find indices for train and test data
    @time for i in eachindex(unique_compounds)
        if (i % 100) == 0
            println("$(round(i/length(unique_compounds)*100, digits = 2))%")
        end
        occurrences = findall(x -> x == unique_compounds[i], filtered_RPLC_data[:,1])
        if unique_compounds[i] in train
            append!(train_indices, occurrences)
        else
            append!(test_indices, occurrences)
        end
    end

    r_f = filtered_RPLC_data[:,end-1]
    # Extract train and test data
    X_train = fingerprints[train_indices, :]
    X_test = fingerprints[test_indices, :]
    y_train = r_f[train_indices]
    y_test = r_f[test_indices]

    return X_train, X_test, y_train, y_test, train_indices, test_indices, y
end
function remove_outliers_IQR(RPLC_data)
    unique_compounds = unique(RPLC_data[:,1])

    non_outliers = []
    for i in eachindex(unique_compounds)
        if i % 100 == 0
            println("$(round(i/length(unique_compounds)*100, digits = 2))%")
        end

        occurrences = findall(x-> x == unique_compounds[i],RPLC_data[:,1])
        group = RPLC_data[occurrences,end]./100
        Q1 = quantile(group, 0.25)
        Q3 = quantile(group, 0.75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        indices_non_outliers = occurrences[(group .>= lower_bound) .& (group .<= upper_bound)]
        append!(non_outliers, indices_non_outliers)
    end

    filtered_RPLC_data = RPLC_data[non_outliers,:]

    return filtered_RPLC_data
end
for i = 1:374
    #Locating the file
    current_file = joinpath("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\processed_data", readdir(folder_path)[i])
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
     Modifier = []
     try 
        run_time = gradient[end,1]
        retention_factors = RT./run_time
        for rtfs in retention_factors
            push!(Modifier, interpolate_B_modifier(rtfs, gradient))
        end
     catch
         retention_factors = repeat(["Not gradient info"], length(RT))
         Modifier = repeat(["Not gradient info"], length(RT))
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
    current_table = DataFrame(Inchikey = Compound_Inchi, LC_mode = LC_mode, MACCS = MACCS, Pubchem_fps = Pubchem_fps,MW = MW, XlogP = XLogP, Retention_factors = retention_factors, Modifier = Modifier)
    Final_table = vcat(Final_table, current_table)
end

Final_table
#Removing compounds that occur both in RPLC and HILIC 
Final_table_unique = remove_overlapping(Final_table)

#Removing rows with missing retention factor and outliers
Final_table_unique = remove_missing_and_outliers(Final_table_unique)

#Formatting MACCS and PubChem Fingerprints for modelling
MACCS, PubChem_fps = format_fingerprints(Final_table_unique)

#Getting indices for RPLC and HILIC rows
indices_RPLC = findall(row -> occursin("RP", row.LC_mode), eachrow(Final_table_unique))

RPLC_data = Final_table_unique[indices_RPLC,:]

#############################################################################################################
#Removing outliers if IQR > 1.5 * threshold
@time filtered_RPLC_data = remove_outliers_IQR(RPLC_data)
filtered_RPLC_data
MACCS, PubChem_fps = format_fingerprints(filtered_RPLC_data)
MACCS
PubChem_fps

#Train/test split without data leakage

X_train, X_test, y_train, y_test, train_indices, test_indices, y = train_test_split_no_leakage_regressor(filtered_RPLC_data, 0.1)
###########################################################################################################
#Random forest classification
Random.seed!(42)

X = [MACCS PubChem_fps]

X = (X[indices_RPLC,:])

y = Float64.(Final_table_unique[:,end][indices_RPLC]./100)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size= 0.1, random_state = 42)

#######################
#Manual Hyperparameter optimisation for ScikitLearn RF
n_estimators_t = 20

criterion = "squared_error"

max_depth_t = [200]

min_samples_split_t = [2]

min_samples_leaf_t = [1]

random_state = 42

f_name = "Random forest RPLC Stratification with strings"

df_results_3 = DataFrame(n_estimators = Int[], min_samples_split = Int[], min_samples_leaf = Int[], train_score = Float64[], test_score = Float64[])

df_results_3 = CSV.read("C:\\Users\\uqthulle\\Documents\\$f_name.csv ", DataFrame)



#StratifiedKFold
bins = 10
k = 3
classes = custom_binning(y[model], bins)

histogram(classes, xticks = (1:bins),
dpi = 300,
label = false,
title = "Binning of retention factors for Stratification")

train, test = Stratifiedcv(X[model,:], classes, k)

for n_estimators in n_estimators_t
    @show n_estimators
        for min_samples_split in min_samples_split_t
            @show min_samples_split
            for min_samples_leaf in min_samples_leaf_t
                @show min_samples_leaf
                scores_train = zeros(k)
                scores_test = zeros(k)
                for fold = collect(1:k)
                    #train_fold = kf.train_indices[fold]
                    #val_fold = kf.val_indices[fold]

                    train_fold = train[fold]
                    val_fold = test[fold]

                    rf_regressor = RandomForestRegressor(n_estimators = n_estimators, criterion = criterion, 
                                                        min_samples_split = min_samples_split, min_samples_leaf = min_samples_leaf, random_state = random_state, n_jobs = -1)

                    ScikitLearn.fit!(rf_regressor, X[model, :][train_fold,:], y[model][train_fold])

                    scores_train[fold] = ScikitLearn.score(rf_regressor, X[model, :][train_fold,:], y[model][train_fold])
                    scores_test[fold] = ScikitLearn.score(rf_regressor, X[model, :][val_fold,:], y[model][val_fold])
                end
                score_train_cv = mean(scores_train)
                score_test_cv = mean(scores_test)
                @show score_train_cv
                @show score_test_cv
                push!(df_results_3, [n_estimators, min_samples_split, min_samples_leaf, score_train_cv, score_test_cv])
            end
    end
end

CSV.write("C:\\Users\\uqthulle\\Documents\\$f_name.csv ", df_results_3)

df_results_3

scatter(df_results_3[8:end,1],df_results_3[8:end,end-1], ylims = (0,1),
xlabel = "number of trees",
ylabel = "CV 3 score",
label ="test",
title = "RPLC RF Sratified CV 3 scores",
dpi = 300)
scatter!(df_results_3[8:end,1],df_results_3[8:end,end],
label = "train")

#Getting the model for the test data using the best hyperparameters
ind_best = argmax(df_results_3.test_score)

#Using the hyperparameter data
rf_regressor = RandomForestRegressor(n_estimators = df_results_3.n_estimators[ind_best], criterion = "squared_error", max_depth = df_results_3.max_depth[ind_best], 
                                     min_samples_split = df_results_3.min_samples_split[ind_best], min_samples_leaf = 1, random_state = 42, n_jobs = -1)

#Manually
rf_regressor = RandomForestRegressor(n_estimators = 20, criterion = "squared_error", 
                                     min_samples_split = 2, min_samples_leaf = 1, random_state = 42, n_jobs = -1)

ScikitLearn.fit!(rf_regressor, X_train, y_train)

y_hat_train = ScikitLearn.predict(rf_regressor,X_train)
y_hat_test = ScikitLearn.predict(rf_regressor,X_test)

score_train = ScikitLearn.score(rf_regressor, X_train, y_train)
score_test = ScikitLearn.score(rf_regressor, X_test, y_test)

scatter(y_train, y_hat_train, label = "train = $(round(score_train, digits= 2))", dpi = 300)
scatter!(y_test,y_hat_test, label = "test = $(round(score_test, digits= 2))")

importance = rf_regressor.feature_importances_

sorted_importance = sortperm(importance, rev = true)

labels = All_keys[sorted_importance]

bar(labels[1:15],sort(importance, rev = true)[1:15],
xrotation=40, 
dpi = 300,
title = "RPLC important variables", 
bottom_margin = 8Plots.mm,
legend = false)

classes_pred = []
for i in eachindex(y_hat_test)
    if y_hat_test[i] >= 0.6
        push!(classes_pred, "Non-mobile")
    elseif p_b[i] <= 0.2
        push!(classes_pred, "Very mobile")
    else 
        push!(classes_pred, "Mobile")
    end
end

confusion_matrix(classes_test, classes_pred)
#leverage

leverage = calculate_leverage(X[model,:], X[test_final,:])

histogram(leverage, xlims = (0,0.35), label = false, dpi = 300,
title = "Leverage calculation Train vs Test set")

percentile_95 = quantile(leverage, 0.95)

vline!([percentile_95], label="95th Percentile", color="red", linestyle=:dash)

#Applicability domain
indices_inside_AD = findall(x-> x <percentile_95, leverage)
indices_outside_AD = findall(x-> x >=percentile_95, leverage)

scatter(y[model], y_hat_train, label = "train = $(round(score_train, digits= 2))", xlims = (0,1), ylims = (0,1), dpi = 300)
scatter!(y[test_final][indices_inside_AD],y_hat_test[indices_inside_AD], label = "test = $(round(score_test, digits= 2))")
scatter!(y[test_final][indices_outside_AD],y_hat_test[indices_outside_AD], label = "Outside AD test", shape = :star)




cd("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\Plots")
savefig("TPR FDR all data.png")
cd("R:\\PHD2024TH-Q6813\\Code\\Regression")


rmse(y_hat_test,y[test_final])