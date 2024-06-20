###Plot to show nonadecanoic acid
unique(filtered_RPLC_data[:,1])
nd_entries = findall(x-> x == "InChI=1S/C19H38O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18-19(20)21/h2-18H2,1H3,(H,20,21)",RPLC_data[:,1])
minimum((RPLC_data[nd_entries,end]./100)[2:end])
scatter((RPLC_data[nd_entries,end]./100), dpi = 300, 
title = "All nonadecanoic acid entries",
xlabel = "Entry n", ylabel = "Retention factor / Φ", label = "Φ")
scatter!(RPLC_data[nd_entries,end-1], label = "Retention factor")

#Plot to show excluded data_compounds
scatter(Float64.(Results_unique[:,end]./100),
 Float64.(Results_unique[:,end-1]), label = false,
dpi = 300, xlabel = "Φ", ylabel = "Retention factor", 
xlims = (0,1), ylims = (0,1.6), legend = :topleft)
plot!([0.1, 0.1], [0.7, 1], c = :red, ls = :dash, label = "excluded")
plot!([-0.5, 0.1], [0.7, 0.7], c = :red, ls = :dash, label = false)
plot!([-0.5, 0.1], [1, 1], c = :red, ls = :dash, label = false) 
hline!([1], c = :red, ls = :dash, label = false)

histogram(filtered_RPLC_data[:,end-1], xlims = (0,1), dpi = 300,
xlabel = "Retention factor",
title = "Distribution of Retention factor values", label = false)
vspan!([0, 0.2], color=:red, alpha=0.3, label = "Very mobile")
vspan!([0.2, 0.6], color=:yellow3, alpha=0.3, 
label = "Mobile")
p1 = vspan!([0.6, 1], color=:green, alpha=0.3,
label = "Non-Mobile",
bottom_margin = 5Plots.mm)


histogram(filtered_RPLC_data[:,end]./100, 
title = "Dsitribution and classification of chemcials according to Φ",
titlefont = 10,
formatter = :plain,
dpi = 300, label = false, xlims = (0,1), xlabel = "Φ",)
vspan!([0, 0.2], color=:red, alpha=0.3, label = "Very mobile")
vspan!([0.2, 0.6], color=:darkorange, alpha=0.3, label = "Mobile")
p2 = vspan!([0.6, 1], color=:green, alpha=0.3, label = "Non-mobile",
bottom_margin = 2.5Plots.mm)

plot(p1,p2, size = (800,400))

#####
#Chemical space of unique chemicals included with Φ
scatter(filtered_RPLC_data[:,5], filtered_RPLC_data[:,6],
zcolor = filtered_RPLC_data[:,end]./100, dpi = 300,
xlabel = "MW (Da)", ylabel = "XlogP", label = false,
colorbar_title  = "Φ", clims = (0,1))


cd("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\Plots")
savefig("Chemical space all data.png")
cd("R:\\PHD2024TH-Q6813\\Code\\Regression")