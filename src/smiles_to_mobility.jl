using CSV, DataFrames, PyCall, Statistics, LinearAlgebra, Conda, ScikitLearn


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

    Smiles_list = String[]
    
    #If there was a match with precomputed smiles, get the mobility from there
    if !isempty(yes_precompiled)
        #Find the indices of the precomputed ones
        precomp = [findfirst(x -> x == string, precompiled_smiles_df.SMILES) for string in yes_precompiled]

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
smiles_to_mobility("CCCCCCCCCCCCCCCCCCCCCCCCCCCOOCCCCCCC")
smiles_to_mobility(path, SMILES_precomputed) 
SMILES = SMILES_precomputed
path = mktempdir()
export smiles_to_mobility

