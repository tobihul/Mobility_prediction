using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots, MultivariateStats
folder_path = "C:\\Users\\uqthulle\\Documents\\RepoRT-master\\raw_data"
Column_modes = CSV.read("C:\\Users\\uqthulle\\Documents\\Unique columns FINAL.csv", DataFrame)
Final_table = DataFrame(Compound_name = String[], canonical_smiles = String[], PubChem_cid = Int64[],  Retention_time = Float64[], Column_used = String31[], LC_mode = String7[])
#explore = CSV.read("C:\\Users\\uqthulle\\Documents\\RepoRT-master\\raw_data\\0001\\0001_gradient.tsv", DataFrame)

for i = 2:376
    #Locating the file
    current_file = joinpath("C:\\Users\\uqthulle\\Documents\\RepoRT-master\\raw_data", readdir(folder_path)[i])
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
    #Check if the column info is missing in the meta data
    column_name = meta_data[1,2]
        if ismissing(column_name) || column_name == "NA"

            column_name_table = repeat(["No column info"], length(compound_names))
        else

            column_name_table = repeat([column_name], length(compound_names))
            
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
    current_table = DataFrame(Compound_name = compound_names, canonical_smiles = compound_canonical_SMILES, PubChem_cid = compound_pchem_cid, Retention_time = compound_retention_times,  Column_used = column_name_table, LC_mode = LC_mode)
    Final_table = vcat(Final_table, current_table)
end

unique(Final_table[:,2])
#CSV.write("C:\\Users\\uqthulle\\Documents\\All compounds RepoRT.csv", Final_table)
Final_table = dropmissing(Final_table, :canonical_smiles)
#Look at the final table 
vscodedisplay(Final_table)

#Check how many of each LC mode there are
index_HILIC = findall(==("HILIC"), Final_table.LC_mode)
unique(Final_table[index_HILIC,2])
index_RPLC = findall(==("RP"), Final_table.LC_mode)
unique(Final_table[index_RPLC,2])
index_no_column = findall(==("No mode info"), Final_table.LC_mode)
#Track how many missing values there are for the categories
values = Final_table.PubChem_cid
p = 0
for i = 1:length(values)
    if ismissing(values[i])
        p+=1
    end
end
p

#Note the indexes of the non missing, non zero and integers of PubChem CIDs
values = Final_table.PubChem_cid
indexes = Int[]
for i = 1:length(values)
    if !ismissing(values[i]) && typeof(values[i]) == Int64 && values[i] != 0
        push!(indexes, i)
    end
end
indexes

#Select only the rows with non-missing CID's
non_missing_cids = values[indexes]
Non_missing_table = Final_table[indexes,:]

#Get the properties for the cids
Df_LogP_MW_pubchem = CSV.File(get_for_cids(non_missing_cids; properties="MolecularWeight,XLOGP", output="CSV")) |> DataFrame
Exact_masses = CSV.File(get_for_cids(non_missing_cids; properties="ExactMass", output="CSV")) |> DataFrame

#Put the properties and original table back together
Df_LogP_MW_pubchem = hcat(Non_missing_table[!, Not("PubChem_cid")],Df_LogP_MW_pubchem,Exact_masses[:,2])
Df_LogP_MW_pubchem = dropmissing(Df_LogP_MW_pubchem, :XLogP)

# Separate entries with RP and HILIC modes
rp_CID = unique(Final_table[Final_table.LC_mode .== "RP", :].canonical_smiles)
hilic_CID = unique(Final_table[Final_table.LC_mode .== "HILIC", :].canonical_smiles)



# Find the intersection of CIDss for both modes
common_CID = intersect(rp_CID, hilic_CID)

common_CID_Df = CSV.File(get_for_cids(common_CID; properties="MolecularWeight,XLOGP", output="CSV")) |> DataFrame

#Removing the ones that have been separated by both as we already have them
condition = in.(Final_table.canonical_smiles, Ref(common_CID))
Df_LogP_MW_pubchem = Final_table[.!condition, :]


######################################################################################################################################################################
#Now do some plotting
#Now separate
indices_RPLC = findall(row -> occursin("RP", row.LC_mode), eachrow(Df_LogP_MW_pubchem))

RPLC = scatter(Df_LogP_MW_pubchem[indices_RPLC,7], Df_LogP_MW_pubchem[indices_RPLC,8], label = "RPLC n = $(length(indices_RPLC))",
dpi = 300,  markerstrokewidth = 0.5, legend = false, xlims = (0,2500), ylims = (-12,25), xticks = ([1000, 2000]))

indices_HILIC = findall(row -> occursin("HILIC", row.LC_mode), eachrow(Df_LogP_MW_pubchem))

HILIC = scatter(Df_LogP_MW_pubchem[indices_HILIC,7], Df_LogP_MW_pubchem[indices_HILIC,8], label = "HILIC n = $(length(indices_HILIC))",
dpi = 300, markershape = :utriangle, c = :darkorange1, markerstrokewidth = 0.5, legend = false, xlims = (0,2500), ylims = (-12,25),xticks = ([1000, 2000]))

Common = scatter(common_CID_Df[:,2], common_CID_Df[:,3], label = "RPLC and HILIC n = $(length(common_CID))",
 dpi = 300,  markerstrokewidth = 0.5, xlims = (0,2500), ylims = (-12,25),
markershape = :star, c = palette(:tab10)[5], legend = false,xticks = ([1000, 2000]))

scatter(Df_LogP_MW_pubchem[indices_RPLC,7], Df_LogP_MW_pubchem[indices_RPLC,8], label = "Only RPLC n = $(length(indices_RPLC))",
xlabel = "MW", ylabel = "XLogP", dpi = 300,  markerstrokewidth = 0.5, xlims = (0,2500), ylims = (-12,25), alpha = 1)

scatter!(Df_LogP_MW_pubchem[indices_HILIC,7], Df_LogP_MW_pubchem[indices_HILIC,8], label = "Only HILIC n = $(length(indices_HILIC))",
dpi = 300, markershape = :utriangle, c = :darkorange1, markerstrokewidth = 0.5, xlims = (0,2500), ylims = (-12,25),xticks = ([1000, 2000]), alpha = 0.5)

p_together =scatter!(common_CID_Df[:,2], common_CID_Df[:,3], label = "RPLC and HILIC n = $(length(common_CID))",
xlabel = "MW", ylabel = "XLogP", dpi = 300,  markerstrokewidth = 0.5, xlims = (0,2500), ylims = (-12,25),
markershape = :star, c = palette(:tab10)[5], alpha = 0.3)

#############
indices_unknown = findall(row -> occursin("No mode info", row.LC_mode), eachrow(Non_missing_table))

Unknown = scatter!(Df_LogP_MW_pubchem[indices_unknown,7], Df_LogP_MW_pubchem[indices_HILIC,8], label = "No column info n = $(length(indices_unknown))",
xlabel = "MW", ylabel = "XLogP", dpi = 300)
#############
l = @layout [
    a{0.7w} [
        grid(3, 1)    
    ]
]

plot(p_together, RPLC, HILIC, Common, layout = l)

savefig("C:\\Users\\uqthulle\\Documents\\PCA with unique RPLC and HILIC entries.png")

#fingerprints
smiles = Final_table.canonical_smiles[indexes][1:10]

indices__over_50 = findall(Non_missing_table.Retention_time[indices_HILIC].<25)

scatter(Df_LogP_MW_pubchem[indices_HILIC,7][indices__over_50], Df_LogP_MW_pubchem[indices_HILIC,8][indices__over_50], 
zcolor = Non_missing_table.Retention_time[indices_HILIC][indices__over_50], color=:RdBu, label = false
, title = "HILIC compounds colored by RT (min)", dpi =300, zticks = 0:5:25, xlabel = "MW", ylabel = "LogP")



histogram(Non_missing_table.Retention_time[indices_RPLC], label = "RPLC RT")
histogram!(Non_missing_table.Retention_time[indices_HILIC], label = "HILIC RT",
 dpi = 300, xlabel = "Rt (min)", title = "Distribution of Rt values for RPLC and HILIC")


#Start by obtaining the mass deffects for the NORMAN database and the detected CEC'S
#This is the function used:
er_masses = [27.9949,46.9689,26.003,43.972,30.9984,13.0078]
function ER_calc_ext(mz_values,m_ru)
    ER=Array{Any}(undef,length(mz_values))
    for i=1:length(mz_values)
        KM=mz_values[i].*(round.(m_ru)./m_ru) # Kendrick mass
        ER[i]=round.(round.(KM) .- KM ; digits=3)
    end
    return ER
end
#NORMAN
Exact_masses_RepoRT = Df_LogP_MW_pubchem[:,9]
K_def = ER_calc_ext(Exact_masses_RepoRT, er_masses)
KMD_RepoRT = Matrix{Float64}(zeros(length(Exact_masses_RepoRT),6))
for j = 1:length(KMD_RepoRT[1,:])
    for i = 1:length(Exact_masses_RepoRT)
        KMD_RepoRT[i,j] = K_def[i][j]
    end
end
KMD_RepoRT

PCA_data = [Df_LogP_MW_pubchem[:,9] Df_LogP_MW_pubchem[:,8] KMD_RepoRT]

#Now we can perform PCA
#using ScikitLearn
#@sk_import decomposition: PCA

PCA_RepoRT = PCA_data
X = PCA_RepoRT
#Mean center and scale since variables all have different magnitudes
X = X .- mean(X,dims = 1)
X = X ./ std(X,dims=1)
#setup PCA model
M = fit(PCA, X; maxoutdim=3)
variance_explained = round.(principalvars(M) ./ tvar(M) * 100, digits = 2)
loadings = predict(M, X)

PCs = ["PC1", "PC2", "PC3"]
var_names = ["MolecularWeight", "XlogP3", "CO","CN","CCl","CS","CF","CH"]
legend_order = ["MolecularWeight", "XlogP3", "CO", "CN", "CF", "CS", "CCl", "CH"]


#Plot the loadings of each variable for each PC in a grouped bar plot:
b1 = bar(loadings[1,:], group = var_names, palette = :tab10, size = (1000,500),
left_margin = 7Plots.mm, bottom_margin = 7.5Plots.mm, right_margin = 5Plots.mm,
ylabel = "Loading", legendfont = font(7), guidefont = 15, xticks = false, title = "PC 1", legend = :topleft, dpi = 300)

b2 = bar(loadings[2,:], group = var_names, palette = :tab10, size = (1000,500),
left_margin = 7Plots.mm, bottom_margin = 7.5Plots.mm, right_margin = 5Plots.mm,
 legendfont = font(11), guidefont = 15, xticks = false, title = "PC 2", legend = false , dpi = 300)

b3 = bar(loadings[3,:], group = var_names, palette = :tab10, size = (1000,500),
left_margin = 7Plots.mm, bottom_margin = 7.5Plots.mm, right_margin = 5Plots.mm,
legendfont = font(5), guidefont = 15, xticks = false, title = "PC 3" , dpi = 300,legend_order = legend_order, legend = :false)

plot(b1,b2,b3, layout = (1,3), dpi = 300)

#Scores
scores = projection(M)
scatter(scores[:,1][indices_RPLC],scores[:,2][indices_RPLC], scores[:,3][indices_RPLC], label = "Only RPLC n = $(length(scores[:,1][indices_RPLC]))", 
seriesalpha = 0.15, markerstrokewidth = 0, ms = 2)
scatter!(scores[:,1][indices_HILIC],scores[:,2][indices_HILIC], scores[:,3][indices_HILIC], label = "Only HILIC n = $(length(scores[:,1][indices_HILIC]))", 
dpi = 300, xlabel = "PC1 ($(variance_explained[1])%)", ylabel = "PC2 ($(variance_explained[2])%)"
, zlabel = "PC3 ($(variance_explained[3])%)", markerstrokewidth = 0, seriesalpha = 0.15, ms = 2, 
title = "PCA of mass defects RPLC and HILIC", legend = :topright, margin=5Plots.mm)

scatter(scores[:,1],scores[:,2], scores[:,3], group = Df_LogP_MW_pubchem[:,5], seriesalpha = 0.05, markerstrokewidth = 0, ms = 2)


distance = sqrt.(X[:,1].^2 + X[:,2].^2 + X[:,3].^2 + X[:,4].^2 + X[:,5].^2 + X[:,6].^2 +X[:,7].^2 + X[:,8].^2) 

histogram(distance[indices_RPLC], label = "RPLC")
histogram!(distance[indices_HILIC], dpi = 300, xlabel = "Euclidean distance of XlogP, MW and KMDs",
label = "HILIC")


