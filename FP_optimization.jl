using FP_optimization
using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots, MultivariateStats
folder_path = "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\raw_data"
Final_table = DataFrame(Compound_name = String[], canonical_smiles = String[], PubChem_cid = Int64[], LC_mode = String7[],Retention_factors = Float64[], Modifier = Float64[])
MACCS_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\MACCS keys.txt")
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")
function remove_overlapping(Results_table)

    # Separate entries with RP and HILIC modes
    rp_Inchikey = Results_table[Results_table.LC_mode .== "RP", :].canonical_smiles
    hilic_Inchikey = Results_table[Results_table.LC_mode .== "HILIC", :].canonical_smiles

    # Find the intersection of CIDss for both modes
    common_Inchikey = intersect(rp_Inchikey, hilic_Inchikey)

    #Removing the ones that have been separated by both as we already have them
    condition = in.(Results_table.canonical_smiles, Ref(common_Inchikey))
    Final_table_unique = Results_table[.!condition, :]

    return Final_table_unique
end
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
for i = 2:376
    #Locating the file
    current_file = joinpath("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\raw_data", readdir(folder_path)[i])
    println(i)
    #Getting the compounds name, canonical SMILES and retention time 
    rt_path = joinpath(current_file, join(["$(readdir(folder_path)[2:end][i-1])", "rtdata.tsv"],"_"))
    data_compounds = CSV.read(rt_path, DataFrame)
    compound_names = data_compounds.name
    compound_canonical_SMILES = data_compounds[:,7]
    compound_retention_times = data_compounds.rt
    compound_pchem_cid = data_compounds[:,5]

    column_path = joinpath(current_file, join(["$(readdir(folder_path)[2:end][i-1])", "metadata.tsv"],"_"))
    meta_data = CSV.read(column_path, DataFrame)

    #Getting the retention factors
    gradient_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "gradient.tsv"],"_"))
    gradient = CSV.read(gradient_path, DataFrame)
    local retention_factors
    Modifier = []
    try 
        run_time = gradient[end,1]
        retention_factors = compound_retention_times./run_time
    for rtfs in retention_factors
        push!(Modifier, interpolate_B_modifier(rtfs, gradient))
    end
    catch
        retention_factors = repeat(["Not gradient info"], length(compound_retention_times))
        Modifier = repeat(["Not gradient info"], length(compound_retention_times))
    end

    #Getting the LC mode info for the data set
    info_path = joinpath(current_file, join(["$(readdir(folder_path)[2:end][i-1])", "info.tsv"],"_"))
    info = CSV.read(info_path, DataFrame)
    #Check if the column info is missing in the meta data
    LC_mode = info[1,3]
    if ismissing(LC_mode) || LC_mode == "NA"

        LC_mode = repeat(["No mode info"], length(compound_names))

    else

        LC_mode = repeat([LC_mode], length(compound_names))

    end
 
    #Getting all the data in one table to append all the other datasets in this format
    current_table = DataFrame(Compound_name = compound_names, canonical_smiles = compound_canonical_SMILES, PubChem_cid = compound_pchem_cid, LC_mode = LC_mode, Retention_factors = retention_factors, Modifier = Modifier)
    Final_table = vcat(Final_table, current_table)
end

Final_table

#Removing missing smiles rows
Final_table = Final_table[.!ismissing.(Final_table.canonical_smiles), :]

#Removing compounds that occur both in RPLC and HILIC 
Final_table_unique = remove_overlapping(Final_table)

#Removing rows with missing retention factor and outliers
Final_table_unique = remove_missing_and_outliers(Final_table_unique)

#Getting indices for RPLC and HILIC rows
indices_RPLC = findall(row -> occursin("RP", row.LC_mode), eachrow(Final_table_unique))

RPLC_data = Final_table_unique[indices_RPLC,:]

#Getting data ready for package
SMILES = string.(RPLC_data.canonical_smiles)

labels = []
for i in eachindex(RPLC_data.Modifier)
    if RPLC_data.Modifier[i] >= 60
        push!(labels, "Non-mobile")
    elseif RPLC_data.Modifier[i] <= 20
        push!(labels, "Very mobile")
    else 
        push!(labels, "Mobile")
    end
end

sum(labels.== "Very mobile")

sum(labels.== "Non-mobile")

sum(labels.== "Mobile")

path = "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\Results fingerprint optimization"

calculate_nonhashed_fps(path, SMILES) 