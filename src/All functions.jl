using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots
using LinearAlgebra, ScikitLearn, Random, MLJ, PyCall, Conda
using ScikitLearn: @sk_import

@sk_import ensemble: RandomForestClassifier
@sk_import model_selection: StratifiedKFold
@sk_import model_selection: train_test_split

# You can access `rf_cl` and `Precompiled_smiles` from the main module
const rf_cl = Mobility_prediction.rf_cl
const Precompiled_smiles = Mobility_prediction.Precompiled_smiles
const pd = Mobility_prediction.pd
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
    rtfs = Results_unique[:,end-1]
    pmbph = Results_unique[:,end]
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
function smiles_to_mobility(SMILES::String)
    #test
    precompiled_smiles_df = Precompiled_smiles[]

    if SMILES in precompiled_smiles_df.SMILES

        index_smiles = findfirst(x-> x == SMILES, precompiled_smiles_df.SMILES)

        predicted_class = precompiled_smiles_df[index_smiles,2]

        predicted_probability = precompiled_smiles_df[index_smiles,3]
    else

        pubchem_fp = try

            DataFrame(pd[].from_smiles(SMILES, fingerprints=true, descriptors = false, timeout = 600, maxruntime = 600))

        catch

            throw("There was either an error calculating the fingerprints or the SMILES is incorrect")
        
        end

        # Extract column names from the DataFrame
        colnames = names(pubchem_fp)

        # Extract the numerical part of the fingerprint column names and sort them since the current order is: FP0, FP1, FP10, FP100, etc.
        sorted_fingerprint_cols = sort(colnames, by = x -> parse(Int, match(r"\d+", x).match))

        # Reorder the DataFrame
        ordered_pubchem_fp = select(pubchem_fp, sorted_fingerprint_cols)

        fingerprints = parse.(Int, Matrix(ordered_pubchem_fp))

        #Predict the class
        predicted_class = ScikitLearn.predict(rf_cl[],fingerprints)

        predicted_probability = Int(round.(maximum(ScikitLearn.predict_proba(rf_cl[],fingerprints))*100, digits = 0))
    end

    return predicted_class, predicted_probability

end
function smiles_to_mobility(path::String,SMILES::String)

    precompiled_smiles_df = Precompiled_smiles[]

    if SMILES in precompiled_smiles_df.SMILES

            index_smiles = findfirst(x-> x == SMILES, precompiled_smiles_df.SMILES)

            predicted_class = precompiled_smiles_df[index_smiles,2]

            predicted_probability = precompiled_smiles_df[index_smiles,3]
    else

        pubchem_fp = try

            DataFrame(pd[].from_smiles(SMILES, fingerprints=true, descriptors = false, timeout = 600, maxruntime = 600))

        catch

            throw("There was either an error calculating the fingerprints or the SMILES is incorrect")
        
        end

        # Extract column names from the DataFrame
        colnames = names(pubchem_fp)

        # Extract the numerical part of the fingerprint column names and sort them since the current order is: FP0, FP1, FP10, FP100, etc.
        sorted_fingerprint_cols = sort(colnames, by = x -> parse(Int, match(r"\d+", x).match))

        # Reorder the DataFrame
        ordered_pubchem_fp = select(pubchem_fp, sorted_fingerprint_cols)

        fingerprints = parse.(Int, Matrix(ordered_pubchem_fp))

        #Predict the class
        predicted_class = ScikitLearn.predict(rf_cl[],fingerprints)

        predicted_probability = Int(round.(maximum(ScikitLearn.predict_proba(rf_cl[],fingerprints))*100, digits = 0))
    end

    df_results = DataFrame(SMILES = SMILES, Predicted_class = predicted_class, Probability = predicted_probability)

    CSV.write(joinpath(path,"Predicted mobility of $(SMILES).csv"), df_results)

    return predicted_class, predicted_probability

end
function smiles_to_mobility(path::String,SMILES::Vector{String})

    precompiled_smiles_df = Precompiled_smiles[]

    #Checking if the fps are already precomputed
    precompiled = [string in precompiled_smiles_df.SMILES for string in SMILES]

    #Splitting into precomputed and not computed yet
    yes_precompiled = [string for (string, is_found) in zip(SMILES, precompiled) if is_found]
    not_precompiled = [string for (string, is_found) in zip(SMILES, precompiled) if !is_found]

    println("$(length(yes_precompiled)) already pre-computed")
    println("$(length(not_precompiled)) need to be manually calculated")
    
    #If there was a match with precomputed smiles, get the mobility from there
    if !isempty(yes_precompiled)
        #Find the indices of the precomputed ones
        precomp = [findall(x -> x == string, precompiled_smiles_df.SMILES) for string in yes_precompiled]

        indices_precomp = vcat(precomp...)

        predicted_class_precomp = precompiled_smiles_df[indices_precomp,2]

        predicted_probability_precomp = precompiled_smiles_df[indices_precomp,3]

        df_comp =  DataFrame(SMILES = yes_precompiled, Predicted_mobility = predicted_class_precomp, Probability = predicted_probability_precomp)
    end

    #Check if there are entries to be manually calculated
    if !isempty(not_precompiled)

        #If there are less than 10 entries run them all at once
        if length(not_precompiled) < 10

            fingerprints = DataFrame()
            Smiles_list = String[]
            try
                    if length(not_precompiled) == 1

                        pubchem_fp = DataFrame.(pd[].from_smiles(not_precompiled, fingerprints=true, descriptors = false, timeout = 600, maxruntime = 600))

                    else 
                        pubchem_fp = DataFrame.(pd[].from_smiles(not_precompiled, fingerprints=true, descriptors = false, timeout = 600, maxruntime = 600))

                    end
                fingerprints = reduce(vcat, pubchem_fp)
                #If they all succeed the list of smiles simply becomes all of them
                Smiles_list = not_precompiled

            catch


                #In case the batch went wrong, do them all individually
                for j = collect(1:length(not_precompiled))
                    
                    try
                        pubchem_fp = DataFrame(pd[].from_smiles(not_precompiled[j], fingerprints=true, descriptors = false, timeout = 600, maxruntime = 600))
                        append!(fingerprints, pubchem_fp)
                        #Keep the smiles that were succesfully calculated
                        push!(Smiles_list , not_precompiled[j])
                    catch
                        #If an entry errors
                        println("$(not_precompiled[j]) has failed")
                    end
                end  
                
            
            end

        else 
            
            #This is when there are more than 10 entries that need to be manually calculated
            #Calculating fingerprints in batches of 10 is around 8 times faster than one-by-one
            
            iterations = collect(1:10:length(not_precompiled)) #Split into gorups of 10
            fingerprints = DataFrame()
            Smiles_list = String[]
            for i in iterations

                    if i != iterations[end]
                        try

                            pubchem_fp = DataFrame.(pd.from_smiles(not_precompiled[i:i+9], fingerprints=true, descriptors = false, timeout = 600, maxruntime = 600))
                            temp_df = reduce(vcat, pubchem_fp)
                            append!(fingerprints, temp_df)
                            append!(Smiles_list, not_precompiled[i:i+9])
                            
                    
                        catch 
                            for j = collect(i:i+9)
                                
                                try
                                    pubchem_fp = DataFrame(pd[].from_smiles(not_precompiled[j], fingerprints=true, descriptors = false, timeout = 600, maxruntime = 600))
                                    append!(fingerprints, pubchem_fp)
                                    push!(Smiles_list, not_precompiled[j])
                                
                                catch
                                    println("$(not_precompiled[j]) has failed")
                                    
                                end
                            end  
                        end
                    else #When the loop reaches the last batch, to avoid having indexing errors if the batch has < 10 entries
                        try

                            pubchem_fp = DataFrame.(pd[].from_smiles(not_precompiled[i:end], fingerprints=true, descriptors = false, timeout = 600, maxruntime = 600))
                            temp_df = reduce(vcat, pubchem_fp)
                            append!(fingerprints, temp_df)
                            append!(Smiles_list, not_precompiled[i:end])
                            
                    
                        catch 
                            for j = collect(i:length(not_precompiled))
                            
                                try
                                    pubchem_fp = DataFrame(pd[].from_smiles(not_precompiled[j], fingerprints=true, descriptors = false, timeout = 600, maxruntime = 600))
                                    append!(fingerprints, pubchem_fp)
                                    push!(Smiles_list, not_precompiled[j])
                                catch
                                    println("$(SMILES[j]) has failed")
                                    
                                end
                            end  
                        end
                    end

                    println("Progress: $(round(length(fingerprints.PubchemFP0)/length(not_precompiled)* 100, digits = 1))%")
            end
        end

        if !isempty(Smiles_list)
            # Extract column names from the DataFrame
            colnames = names(fingerprints)

            # Extract the numerical part of the fingerprint column names and sort them since the current order is: FP0, FP1, FP10, FP100, etc.
            sorted_fingerprint_cols = sort(colnames, by = x -> parse(Int, match(r"\d+", x).match))

            # Reorder the DataFrame
            ordered_pubchem_fp = select(fingerprints, sorted_fingerprint_cols)

            fingerprints = parse.(Int, Matrix(ordered_pubchem_fp))

            #Predict the class
            predicted_class_not_comp = ScikitLearn.predict(rf_cl[],fingerprints)

            #Get the probability
            predicted_probability_not_comp = vec(Int.(round.(maximum(ScikitLearn.predict_proba(rf_cl[],fingerprints), dims = 2)*100, digits = 0)))

            df_not_comp = DataFrame(SMILES = Smiles_list, Predicted_mobility = predicted_class_not_comp, Probability = predicted_probability_not_comp)

        end

    end

    if length(yes_precompiled) > 0 && length(Smiles_list) > 0

        df_results = vcat(df_comp, df_not_comp)

    elseif length(yes_precompiled) > 0 && length(Smiles_list) == 0
    
        df_results = df_comp
        
    else

        df_results = df_not_comp
        
    end

    CSV.write(joinpath(path,"Batch predicted mobility.csv"), df_results)

    return df_results

end

export interpolate_B_modifier, calculate_leverage, remove_missing_and_outliers, TPR_FDR, 
train_test_split_no_leakage_classifier, assign_labels, strat_labels, remove_outliers_IQR, remove_missing_identifiers,
append_occurrences_folds, smiles_to_mobility

