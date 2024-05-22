using CSV, Statistics, DataFrames, StatsPlots
folder_path = "C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data"
Final_table = DataFrame(Inchikey = String[], LC_mode = String7[], MACCS = String[], Pubchem_fps = String[], XlogP = Float64[], Retention_factors = Float64[])

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

Final_table_unique
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

using Random, DecisionTree, MLJ
Random.seed!(42)
X = [MACCS Pubchem_fps]

#Selecting either HILIC or RPLC data
y = Final_table_unique[indices_HILIC,6]
X = X[indices_HILIC,:]

#Shuffling the data
shuffle_indices = shuffle(collect(1:length(y)))
X = X[shuffle_indices,:]
y= y[shuffle_indices]


train, test = partition(eachindex(y), 0.8, rng=42)

df_results_2 = DataFrame(n_subfeatures = Int[],n_trees = Int[],partial_sampling=Float32[],max_depth = Int[],min_samples_leaf=Int[],min_samples_split = Int[],min_purity_increase = Float32[],score=Float32[])
df_results_2 = CSV.read("C:\\Users\\uqthulle\\Documents\\$file_name.csv ", DataFrame)

n_folds = 3
n_subfeatures_t= [1047]
n_trees_t=[100]
partial_sampling_t= [0.7]
max_depth_t=[-1]
min_samples_leaf_t=[2]
min_samples_split_t=[2]
min_purity_increase_t=[0.0]
seed=42
file_name = "DecisionTree.jl Results HILIC"   
iteration = 1
for n_subfeatures in n_subfeatures_t
    for n_trees in n_trees_t
        for partial_sampling in partial_sampling_t
            for max_depth in max_depth_t
                for min_samples_leaf in min_samples_leaf_t
                    for min_samples_split in min_samples_split_t
                        for min_purity_increase in min_purity_increase_t


                            r2 =  nfoldCV_forest(y, X,
                                                n_folds,
                                                n_subfeatures,
                                                n_trees,
                                                partial_sampling,
                                                max_depth,
                                                min_samples_leaf,
                                                min_samples_split,
                                                min_purity_increase;
                                                verbose = false,
                                                rng = seed)
                                            score = mean(r2)
                                            push!(df_results_2,[n_subfeatures,n_trees,partial_sampling,max_depth,min_samples_leaf,min_samples_split,min_purity_increase,score])
                        end
                    end
                end
            end
        end
    end
end    
                                            
                                    vscodedisplay(df_results_2)       
                                            CSV.write("C:\\Users\\uqthulle\\Documents\\$file_name.csv ", df_results_2)
