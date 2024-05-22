using CSV, Statistics, DataFrames, PubChemCrawler, MultivariateStats, StatsBase, GLMakie
import AbstractPlotting: pixelarea

folder_path = "C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data"
Final_table = DataFrame(Name = String[],MW = Float64[], XLogP = Float64[],  LC_mode = String7[], Class = String[])

for i = 1:374
    #Locating the file
    current_file = joinpath("C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data", readdir(folder_path)[i])
    println(i)

    #Getting the class
    rt_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "rtdata_canonical_success.tsv"],"_"))
    data_compounds = CSV.read(rt_path, DataFrame)
    names = data_compounds[:,2]
    class = data_compounds[:,9]
    for c = 1:length(class)
        if ismissing(class[c]) || class == "NA"
            class[c] = "No class info"
        end
    end
    
    #Getting MW and XLogP
    descriptor_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "descriptors_canonical_success.tsv"],"_"))
    data_descriptors = CSV.read(descriptor_path, DataFrame)
    MW = data_descriptors.MW
    XLogP = data_descriptors.XLogP

    #Getting the LC mode info for the data set
    info_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "info.tsv"],"_"))
    info = CSV.read(info_path, DataFrame)
    #Check if the column info is missing in the meta data
    LC_mode = info[1,3]
    if ismissing(LC_mode) || LC_mode == "NA"

        LC_mode = repeat(["No mode info"], length(data_compounds[:,1]))

    else

        LC_mode = repeat([LC_mode], length(data_compounds[:,1]))

    end
   
  
    #Getting all the data in one table to append all the other datasets in this format
    current_table = DataFrame(Name = names, MW = MW, XLogP = XLogP, LC_mode = LC_mode, Class = class)
    Final_table = vcat(Final_table, current_table)
end

Final_table = dropmissing(Final_table)
unique(Final_table[:,5])
vscodedisplay(Final_table)

rp_Inchikey = Final_table[Final_table.LC_mode .== "RP", :].Name
hilic_Inchikey = Final_table[Final_table.LC_mode .== "HILIC", :].Name

# Find the intersection of CIDss for both modes
common_Inchikey = intersect(rp_Inchikey, hilic_Inchikey)

#Removing the ones that have been separated by both as we already have them
condition = in.(Final_table.Name, Ref(common_Inchikey))
Final_table = Final_table[.!condition, :]

scatter(Final_table[:,2], Final_table[:,3], group = Final_table[:,5], legend = false)

indices_RPLC = findall(row -> occursin("RP", row.LC_mode), eachrow(Final_table))
indices_HILIC = findall(row -> occursin("HILIC", row.LC_mode), eachrow(Final_table))

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
Exact_masses_RepoRT = Final_table[:,2]
K_def = ER_calc_ext(Exact_masses_RepoRT, er_masses)
KMD_RepoRT = Matrix{Float64}(zeros(length(Exact_masses_RepoRT),6))
for j = 1:length(KMD_RepoRT[1,:])
    for i = 1:length(Exact_masses_RepoRT)
        KMD_RepoRT[i,j] = K_def[i][j]
    end
end
KMD_RepoRT

X = [Final_table[:,2] Final_table[:,3] KMD_RepoRT] 
X = X .- mean(X,dims = 1)
X = X ./ std(X,dims=1)
#setup PCA model
M = fit(PCA, X; maxoutdim=8)

variance_explained = round.(principalvars(M) ./ tvar(M) * 100, digits = 2)
loadings = predict(M, X)
gorder = Final_table[:,5]
scores = projection(M)

groups = Final_table[:,4]
group_counts = countmap(groups)
sorted_groups = sort(collect(keys(group_counts)), by=x->group_counts[x], rev=true)
scatter(title = "PCA plot using RepoRT descriptor for KMD") # Use default plot attributes
for group in sorted_groups
    indices = findall(x -> x == group, groups)
    scatter!(scores[indices, 1], scores[indices, 2], label=group, markersize=3) # Adjust markersize as needed
end

# Show the plot
plot!(legend=:topright, dpi = 300, xlabel = "PC1 ($(variance_explained[1])%)", ylabel = "PC2 ($(variance_explained[2])%)")
savefig("C:\\Users\\uqthulle\\Documents\\Plots\\All classes plot 2.png")

groups = Final_table[:,5]
group_counts = countmap(groups)
sorted_groups = sort(collect(keys(group_counts)), by=x->group_counts[x], rev=true)
scatter(title = "PCA plot using RepoRT descriptor for KMD") # Use default plot attributes
for group in sorted_groups
    indices = findall(x -> x == group, groups)
    scatter!(scores[indices, 1], scores[indices, 2], label=group, markersize=3) # Adjust markersize as needed
end
display()
# Show the plot
plot!(legend=:false, dpi = 300, xlabel = "PC1 ($(variance_explained[1])%)", ylabel = "PC2 ($(variance_explained[2])%)")



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

groups = Final_table[:,5]

group_counts = countmap(groups)

sorted_groups = sort(collect(keys(group_counts)), by=x->(x == group_to_plot, group_counts[x]))
plot = Figure(Size = (12000,8000))
group_to_plot = sorted_groups[3]

indices_group = findall(x -> x == group_to_plot, groups)

indices_all = findall(x -> x != group_to_plot, groups)

ax = Axis(plot[1, 1],
    title = "$group_to_plot",
    titlesize = 10,
    xlabel = "PC1",
    ylabel = "PC2"
)
scatter!(ax,scores[indices_all, 1], scores[indices_all, 2], markersize=5)
scatter!(ax,scores[indices_group, 1], scores[indices_group, 2], markersize=5)

plot

save("C:\\Users\\uqthulle\\Documents\\Plots\\All classes plot.mki", plot, px_per_unit = 3)