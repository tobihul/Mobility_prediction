using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots, MultivariateStats, Images
using GLM, TypedTables, LinearAlgebra, ScikitLearn, Random, MLDataUtils
using ScikitLearn: @sk_import
@sk_import ensemble: RandomForestClassifier
@sk_import model_selection: StratifiedKFold
@sk_import model_selection: train_test_split
folder_path = "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\processed_data"
Final_table = DataFrame(Inchikey = String[], LC_mode = String7[], MACCS = String[], Pubchem_fps = String[],MW = Float64[], XlogP = Float64[], Retention_factors = Float64[], Modifier = Float64[] )
MACCS_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\MACCS keys.txt")
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")
All_keys =  [MACCS_keys ;PubChem_keys]
All_keys_log_MW = [MACCS_keys ;PubChem_keys;"MW"]
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
    b = inv(X'*X)
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
        if rtfs[i] > 0.7 && pmbph[i] < 10
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
        index_first_occurrence[i] = findfirst(x-> x == unique_vector[i], original_data[:,1])

    end

    strat_labels = labels[index_first_occurrence]
    return strat_labels
end
function train_test_split_no_leakage_classifier(filtered_RPLC_data, split)
    #fingerprints = [PubChem_fps filtered_RPLC_data[:,5]]

    #fingerprints = fingerprints .- mean(fingerprints,dims = 1)

    #fingerprints = fingerprints ./ std(fingerprints,dims=1)

    #fingerprints = replace(fingerprints, NaN => 0.0)

    fingerprints = PubChem_fps

    unique_compounds = unique(filtered_RPLC_data[:,1])

    p_b = filtered_RPLC_data[:,end]./100

    y = assign_labels(p_b)

    first_labels = strat_labels(unique_compounds,filtered_RPLC_data[:,1], y)
    
    train, test = train_test_split(unique_compounds, test_size = split, random_state = 42, stratify = first_labels, shuffle = true)
    
    # Initialize variables
    train_indices = []
    test_indices = []

    # Find indices for train and test data
    @time for i in eachindex(unique_compounds)
        if (i % 100) == 0
            println("$(round(i/length(unique_compounds)*100, digits = 2))%")
        end
        occurrences = findall(x -> x == unique_compounds[i], filtered_RPLC_data[:, 1])
        if unique_compounds[i] in train
            append!(train_indices, occurrences)
        else
            append!(test_indices, occurrences)
        end
    end

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
         retention_factors = repeat(["No gradient info"], length(RT))
         Modifier = repeat(["No gradient info"], length(RT))
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

#Removing rows with missing retention factor and outliers
Final_table_unique = remove_missing_and_outliers(Final_table)

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

#############################################################################################################
#Train/test split without data leakage

X_train, X_test, y_train, y_test, train_indices, test_indices, y = train_test_split_no_leakage_classifier(filtered_RPLC_data, 0.1)

##modelling
rf_cl = RandomForestClassifier(n_estimators = 100, max_features = 0.5,
                                   max_depth = 100, min_samples_split = 2,
                                   min_samples_leaf = 2,
                                   random_state = 42, class_weight = "balanced", n_jobs = -1)


ScikitLearn.fit!(rf_cl, X_train, y_train)

#####Getting REACH FPs
using CSV, DataFrames

# Read REACH data
REACH_data = CSV.read("R:\\PHD2024TH-Q6813\\Models and other documents\\PubChem_REACH_list.csv", DataFrame)
non_missing_REACH_data = dropmissing!(REACH_data, :xlogp)
REACH_INCHI = Set(non_missing_REACH_data.inchi)  # Use a Set for efficient removal
Comptox_path = "C:\\Users\\uqthulle\\Downloads\\All Comptox FPS"
files = readdir(Comptox_path)

# Initialize results with the correct column types based on the first file
first_file_path = joinpath(Comptox_path, files[1])
first_file = CSV.read(first_file_path, DataFrame)
results = DataFrame(first_file[1:0, :])[:,[5;8;5829:6709]]  # Create an empty DataFrame with the same structure as the first file

# Track matched INCHIKEYs
matched_inchis = Set{String}()

for file in files
    @show file
    file_path = joinpath(Comptox_path, file)
    current_file = CSV.read(file_path, DataFrame)
    match_indices = findall(x -> x in REACH_INCHI, unique(current_file.INCHI_STRING_DTXCID))
    fingerprints = current_file[match_indices, :][:,[5;8;5829:6709]] 
    
    if !isempty(match_indices)
        fingerprints = current_file[match_indices, :][:, [5;8;5829:6709]]
        results = vcat(results, fingerprints)
        
        # Remove matched INCHIKEYs from REACH_INCHI
        matched_inchis = current_file.INCHI_STRING_DTXCID[match_indices]
        REACH_INCHI = setdiff!(REACH_INCHI, Set(matched_inchis))
    end

end
CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\REACH_fingerprints_PubChem.csv", results)
results

REACH_results = CSV.read("R:\\PHD2024TH-Q6813\\Models and other documents\\REACH_fingerprints_PubChem.csv", DataFrame)


REACH_P_FP = REACH_results[:,2:end]

# Function to extract the numeric part of the column names
function extract_numeric(col_name::String)
    return parse(Int, match(r"\d+", col_name).match)
end
extract_numeric(REACH_P_FP)
# Sort the column names based on the numeric part
sorted_columns = sort(names(REACH_P_FP), by=extract_numeric)

# Reorder the DataFrame columns
REACH_FPs = Matrix(dropmissing(REACH_P_FP[:, Symbol.(sorted_columns)]))

score_train = ScikitLearn.score(rf_cl, X_train, y_train)
score_test = ScikitLearn.score(rf_cl, X_test, y_test)
REACH_hat = ScikitLearn.predict(rf_cl, REACH_FPs)


Very_mobile = sum(REACH_hat.=="Very mobile")
Mobile = sum(REACH_hat.=="Mobile")
Non_mobile = sum(REACH_hat.=="Non-mobile")

bar([1,2,3],[Very_mobile,Mobile, Non_mobile],
xticks = (1:3,["Very mobile" "Mobile" "Non-mobile"]), formatter = :plain,
legend = false, yticks = false, yaxis = false, dpi = 300,
 grid = false, color =[:red, :darkorange, :green], alpha = 0.3)
 annotate!([1,2,3], [Very_mobile-1000, Mobile-1000, Non_mobile-1000],
 [text("$Very_mobile ($(round(Very_mobile/length(REACH_CAS)*100, digits=1))%)", 12, :center, :black),
  text("$Mobile ($(round(Mobile/length(REACH_CAS)*100, digits=1))%)", 12, :center, :black),
  text("$Non_mobile ($(round(Non_mobile/length(REACH_CAS)*100, digits=1))%)", 12, :center, :black)])

 round(Very_mobile/length(REACH_CAS)*100, digits = 2)



 function cid_to_cas(cid::Int64)
    pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$cid/xrefs/RN/TXT"
    response = HTTP.get(pubchem_url)
    cas_numbers = split(String(response.body), "\n")
    return cas_numbers[1]
end
function cas_gen(Cids::Vector{Int64})
    CAS = Vector{String}(undef, length(Cids))
    for i = 1:length(Cids)
        CAS[i] = cid_to_cas(Cids[i])
        @show i
    end
    return CAS
end
REACH_CID = non_missing_REACH_data.cid
REACH_CAS = cas_gen(REACH_CID)
csv.write("R:\\PHD2024TH-Q6813\\Models and other documents\\REACH_CAS.csv", REACH_CAS)