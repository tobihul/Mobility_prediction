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
    c_matrix= ConfusionMatrices.matrix(c_matrix)
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
function train_test_split_no_leakage_classifier(data::DataFrame, split::Float64, fingerprints::Matrix)
    # Step 1: Identify unique compounds and their majority class
    unique_compounds = unique(data.SMILES)
    y = data.Class

    # Step 2: Get majority labels (already done outside)
    first_labels = strat_labels(unique_compounds, data.SMILES, y)

    # Step 3: Stratified train/test split
    train, test = train_test_split(unique_compounds, test_size=split, stratify=first_labels, random_state = 42, shuffle=true)

    # Step 4: Map SMILES to indices
    compound_to_indices = Dict{String, Vector{Int}}()
    for (i, smiles) in pairs(data.SMILES)
        push!(get!(compound_to_indices, smiles, Int[]), i)
    end

    # Step 5: Collect train/test indices
    train_indices = vcat([compound_to_indices[smiles] for smiles in train]...)
    test_indices = vcat([compound_to_indices[smiles] for smiles in test]...)

    # Step 6: Shuffle and extract data
    Random.seed!(42)
    shuffle!(train_indices)
    shuffle!(test_indices)

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
    # Map each unique value in original_data to its first index
    first_occurrence = Dict{eltype(original_data), Int}()
    for (i, val) in pairs(original_data)
        if !haskey(first_occurrence, val)
            first_occurrence[val] = i
        end
    end

    # Use the dictionary to quickly get first indices for unique_vector
    index_first_occurrence = Vector{Int}(undef, length(unique_vector))
    for i in eachindex(unique_vector)
        if (i % 100) == 0
            println("$(round(i/length(unique_vector)*100, digits = 2))%")
        end
        index_first_occurrence[i] = first_occurrence[unique_vector[i]]
    end

    return labels[index_first_occurrence]
end
function remove_outliers_IQR(RPLC_data)
    unique_compounds = unique(RPLC_data.SMILES)
    smiles_col = RPLC_data.SMILES

    non_outliers = []

    for (i, compound) in enumerate(unique_compounds)
        if i % 100 == 0
            println("$(round(i / length(unique_compounds) * 100, digits = 2))%")
        end

        occurrences = findall(x -> x == compound, smiles_col)
        group = RPLC_data.Modifier[occurrences] ./ 100

        Q1 = quantile(group, 0.25)
        Q3 = quantile(group, 0.75)
        IQR = Q3 - Q1

        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR

        indices_non_outliers = occurrences[(group .>= lower_bound) .& (group .<= upper_bound)]
        append!(non_outliers, indices_non_outliers)
    end

    return RPLC_data[non_outliers, :]
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
function append_occurrences_folds(unique_train::Vector{String}, fold::Int, train_indices::Vector{Int})
    train_set = Set(unique_train[train_fold[fold]])

    # Map from SMILES to global indices in merged_data (subset of train_indices)
    SMILES_to_indices = Dict{String, Vector{Int}}()
    for i in train_indices
        SMILES = merged_data.SMILES[i]
        push!(get!(SMILES_to_indices, SMILES, Int[]), i)
    end

    # Map from global index to local index in train_indices
    global_to_local = Dict{Int, Int}()
    for (local_idx, global_idx) in enumerate(train_indices)
        global_to_local[global_idx] = local_idx
    end

    train_indices_fold_local = Int[]
    test_indices_fold_local = Int[]

    for (i, compound) in enumerate(unique_train)
        if (i % 10) == 0
            println("$(round(i / length(unique_train) * 100, digits=2))%")
        end

        indices_global = get(SMILES_to_indices, compound, Int[])

        # Convert to local indices and skip if not found (shouldn't happen if consistent)
        indices_local = [global_to_local[idx] for idx in indices_global if haskey(global_to_local, idx)]

        if compound in train_set
            append!(test_indices_fold_local, indices_local)
        else
            append!(train_indices_fold_local, indices_local)
        end
    end

    shuffle!(train_indices_fold_local)
    shuffle!(test_indices_fold_local)

    return train_indices_fold_local, test_indices_fold_local
end
export interpolate_B_modifier, calculate_leverage, remove_missing_and_outliers, TPR_FDR, 
train_test_split_no_leakage_classifier, assign_labels, strat_labels, remove_outliers_IQR, remove_missing_identifiers,
append_occurrences_folds