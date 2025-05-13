using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots
using LinearAlgebra, ScikitLearn, Random, MLJ, PyCall, Conda
using ScikitLearn: @sk_import

joblib = pyimport("joblib")
np = pyimport("numpy")
sklearn = pyimport("sklearn")
@sk_import ensemble: RandomForestClassifier
@sk_import model_selection: StratifiedKFold
@sk_import model_selection: train_test_split


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
    rtfs = Results_unique.Retention_factors
    pmbph = Results_unique.Modifier
    indexes = Int[]
    for i = 1:length(rtfs)
        if rtfs[i] > 0.7 && pmbph[i] < 10  || rtfs[i] > 1
           # continue
        else
            push!(indexes, i)
        end
    end
    indexes
    Results_unique = Results_unique[indexes,:]
    return Results_unique
end
function TPR_FDR(c_matrix)
    #c_matrix= ConfusionMatrices.matrix(c_matrix)
    TP_M = c_matrix[1,1]
    FN_M = c_matrix[2,1] + c_matrix[3,1]
    FP_M = c_matrix[1,2] + c_matrix[1,3]
    TPR_M = TP_M/(TP_M + FN_M)
    FDR_M = FP_M/(TP_M + FP_M)
    F1_M = (2*TP_M)/(2*TP_M + FP_M + FN_M)

    TP_NM = c_matrix[2,2]
    FN_NM = c_matrix[1,2] + c_matrix[3,2]
    FP_NM = c_matrix[2,1] + c_matrix[2,3]
    TPR_NM = TP_NM/(TP_NM + FN_NM)
    FDR_NM = FP_NM/(TP_NM + FP_NM)
    F1_NM = (2*TP_NM)/(2*TP_NM + FP_NM + FN_NM)

    TP_VM = c_matrix[3,3]
    FN_VM = c_matrix[1,3] + c_matrix[2,3]
    FP_VM = c_matrix[3,1] + c_matrix[3,2]
    TPR_VM = TP_VM/(TP_VM + FN_VM)
    FDR_VM = FP_VM/(TP_VM + FP_VM)
    F1_VM = (2*TP_VM)/(2*TP_VM + FP_VM + FN_VM)
    
    Classes = ["Mobile", "Non-Mobile", "Very mobile"]
    TPRs = [TPR_M, TPR_NM, TPR_VM]
    FDRs = [FDR_M, FDR_NM, FDR_VM]
    F1_scores = [F1_M, F1_NM, F1_VM]
    df = DataFrame(Class = Classes, TPR = TPRs, FDR = FDRs, F1_scores = F1_scores)

    return df
end
function train_test_split_no_leakage_classifier(filtered_RPLC_data, split, fingerprints)
   

    unique_compounds = unique(filtered_RPLC_data.InChi)

    p_b = filtered_RPLC_data.Modifier./100

    y = assign_labels(p_b)

    first_labels = strat_labels(unique_compounds,filtered_RPLC_data.InChi, y)
    
    train, test = train_test_split(unique_compounds, test_size = split, random_state = 42, stratify = first_labels, shuffle = true)
    
    # Initialize variables
    train_indices = []
    test_indices = []

    # Find indices for train and test data
    @time for i in eachindex(unique_compounds)
        if (i % 100) == 0
            println("$(round(i/length(unique_compounds)*100, digits = 2))%")
        end
        occurrences = findall(x -> x == unique_compounds[i], filtered_RPLC_data.InChi)
        if unique_compounds[i] in train
            append!(train_indices, occurrences)
        else
            append!(test_indices, occurrences)
        end
    end
    Random.seed!(42)
     # Shuffle the indices
     shuffle!(train_indices)
     shuffle!(test_indices)
    # Extract train and test data
    X_train = fingerprints[train_indices, :]
    X_test = fingerprints[test_indices, :]
    y_train = y[train_indices]
    y_test = y[test_indices]

    return X_train, X_test, y_train, y_test, train_indices, test_indices, y
end
function assign_labels(Φ)
    y = String[]
    for i in eachindex(Φ)
        if Φ[i] >= 0.6
            push!(y, "Non-mobile")
        elseif Φ[i] <= 0.2
            push!(y, "Very mobile")
        else 
            push!(y, "Mobile")
        end
    end
    return y
end
function strat_labels(unique_vector, original_data, labels)

    index_first_occurrence = Vector{Int}(undef, length(unique_vector))

    for i in eachindex(unique_vector)

        if (i % 100) == 0
            println("$(round(i/length(unique_vector)*100, digits = 2))%")
        end
        index_first_occurrence[i] = findfirst(x-> x == unique_vector[i],original_data)

    end

    strat_labels = labels[index_first_occurrence]
    return strat_labels
end
function remove_outliers_IQR(RPLC_data)
    unique_compounds = unique(RPLC_data.SMILES)

    non_outliers = []
    for i in eachindex(unique_compounds)
        if i % 100 == 0
            println("$(round(i/length(unique_compounds)*100, digits = 2))%")
        end

        occurrences = findall(x-> x == unique_compounds[i],RPLC_data.SMILES)
        group = RPLC_data.Modifier[occurrences]./100
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
function remove_missing_identifiers(RPLC_data)      
    InChI = RPLC_data.InChi
    indexes = Int[]
    for i = 1:length(InChI)
        if !ismissing(InChI[i])
            push!(indexes, i)
        end
    end
    indexes
    RPLC_data = RPLC_data[indexes,:]

    SMILES = RPLC_data.SMILES
    indexes = Int[]
    for i = 1:length(SMILES)
        if !ismissing(SMILES[i])
            push!(indexes, i)
        end
    end
    indexes
    RPLC_data = RPLC_data[indexes,:]
    
    return RPLC_data
end
function append_occurrences_folds(unique_train::Vector{String}, fold::Int)
    train_set = Set(unique_train[train_fold[fold]])
    # Initialize variables
    train_indices_fold = Int[]
    test_indices_fold = Int[]

    # Find indices for train and test data
    @time for i in eachindex(unique_train)
        if (i % 10) == 0
            println("$(round(i/length(unique_train)*100, digits = 2))%")
        end
        occurrences = findall(x -> x == unique_train[i], filtered_RPLC_data.InChi[train_indices])
        if unique_train[i] in train_set
            append!(test_indices_fold, occurrences)
        else
            append!(train_indices_fold, occurrences)
        end
    end

    shuffle!(train_indices_fold)
    shuffle!(test_indices_fold)
    return train_indices_fold, test_indices_fold

end

export interpolate_B_modifier, calculate_leverage, remove_missing_and_outliers, TPR_FDR, 
train_test_split_no_leakage_classifier, assign_labels, strat_labels, remove_outliers_IQR, remove_missing_identifiers,
append_occurrences_folds