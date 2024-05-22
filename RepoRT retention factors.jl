using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots, MultivariateStats

folder_path = "C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data"
Final_table = DataFrame(Names = String[], MW  = Float64[], XlogP = Float64[], LC_mode = String[], Retention_factor = Float64[], Modifier = Float64[])
function interpolate_B_modifier(time::Float64, gradient::DataFrame)
    times = gradient[:,1]./run_time
    idx = searchsortedlast(times, time)
    if idx == 0
        return times[1]
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

for i = 1:374
    #Locating the file
    current_file = joinpath("C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data", readdir(folder_path)[i])
    println(i)
    #Getting the SMILES and retention times
    rt_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "rtdata_canonical_success.tsv"],"_"))
    data_compounds = CSV.read(rt_path, DataFrame)
    RT = data_compounds[:,4]
    Namess = data_compounds[:,2]
    
    #Getting the retention factors
    gradient_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "gradient.tsv"],"_"))
    gradient = CSV.read(gradient_path, DataFrame)
    local retention_factors
    local Modifier
    try 
        run_time = gradient[end,1]
        retention_factors = RT./run_time
        for rtfs in retention_factors
            Modifier = interpolate_B_modifier(rtfs, gradient)
        end
    catch
        retention_factors = repeat(["Not gradient info"], length(RT))
        Modifier = repeat(["Not gradient info"], length(RT))
    end
    

    #Getting MW and XLogP
    descriptor_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "descriptors_canonical_success.tsv"],"_"))
    data_descriptors = CSV.read(descriptor_path, DataFrame)
    MWs = data_descriptors.MW
    XLogPs = data_descriptors.XLogP

    #Getting the LC mode info for the data set
    info_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "info.tsv"],"_"))
    info = CSV.read(info_path, DataFrame)
    #Check if the column info is missing in the meta data
    LC_mode = info[1,3]
    if ismissing(LC_mode) || LC_mode == "NA"

        LC_mode = repeat(["No mode info"], length(RT))

    else

        LC_mode = repeat([LC_mode], length(RT))

    end

    
  
    #Getting all the data in one table to append all the other datasets in this format
    current_table = DataFrame(Names = Namess, MW = MWs, XlogP = XLogPs, LC_mode = LC_mode, Retention_factor = retention_factors, Modifier = Modifier)
    Final_table = vcat(Final_table, current_table)
end



Final_table
Final_table = dropmissing(Final_table)
values = Final_table[:,5]
indexes = Int[]
for i = 1:length(values)
    if typeof(values[i]) == Float64
        push!(indexes, i)
    end
end
indexes
Final_table = Final_table[indexes,:]

rp_Inchikey = Final_table[Final_table.LC_mode .== "RP", :].Names
hilic_Inchikey = Final_table[Final_table.LC_mode .== "HILIC", :].Names

# Find the intersection of CIDss for both modes
common_Inchikey = intersect(rp_Inchikey, hilic_Inchikey)

#Removing the ones that have been separated by both as we already have them
condition = in.(Final_table.Names, Ref(common_Inchikey))
Final_table = Final_table[.!condition, :]


indices_RPLC = findall(row -> occursin("RP", row.LC_mode), eachrow(Final_table))
indices_HILIC = findall(row -> occursin("HILIC", row.LC_mode), eachrow(Final_table))

scatter(Final_table[indices_RPLC,end], bins = 50)

scatter(Final_table[indices_HILIC,2], Final_table[indices_HILIC,3], marker_z = Final_table[indices_HILIC,5],xlims = (0,1500), 
 color = :plasma, size = (1280, 720), grid = true, markerstrokewidth = 0.5, label = false,
xlabel = "Molecular weight (Da)", ylabel = "XLogP",
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm, xtickfont=font(13), 
ytickfont=font(13), 
guidefont=font(20), 
legendfont=font(13), markersize = 6, dpi = 300, title = "HILIC Retention factors",
clims = (0,1))

scatter(Final_table[indices_RPLC,2], Final_table[indices_RPLC,3], marker_z = Final_table[indices_RPLC,5],xlims = (0,1500), 
 color = :plasma, size = (1280, 720), grid = true, markerstrokewidth = 0.5, label = false,
xlabel = "Molecular weight (Da)", ylabel = "XLogP",
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm, xtickfont=font(13), 
ytickfont=font(13), 
guidefont=font(20), 
legendfont=font(13), markersize = 6, dpi = 300, title = "RPLC Retention factors",
clims = (0,1))

savefig("C:\\Users\\uqthulle\\Documents\\Plots\\HILIC retention factor PCA plot.png")

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


M = fit(PCA, X; maxoutdim=8)

variance_explained = round.(principalvars(M) ./ tvar(M) * 100, digits = 2)
loadings = MultivariateStats.predict(M, X)
gorder = Final_table[:,5]
scores = projection(M)

scatter(scores[indices_HILIC,1], scores[indices_HILIC,2], marker_z = Final_table[indices_HILIC,5], 
color = :plasma, size = (1280, 720), grid = true, markerstrokewidth = 0.5, label = false,
xlabel = "PC 1", ylabel = "PC 2",
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm, xtickfont=font(13), 
ytickfont=font(13), 
guidefont=font(20), 
legendfont=font(13), markersize = 6, dpi = 300, title = "PCA retention factors", clims = (0,1))