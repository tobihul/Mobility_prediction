include("All functions.jl")

###Load in the REACH data
REACH_data = CSV.read("R:\\PHD2024TH-Q6813\\Models and other documents\\REACH fingerprints final.csv", DataFrame)

###Looking at the data distribution
MW_REACH = histogram(REACH_data.mw,  
linecolor=:transparent,
 xlims = (0,1500),
dpi = 300,
xlabel = ("Molecular Weight (Da)"),
label = false,
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency")

XlogP_REACH = histogram(REACH_data.xlogp,  
bins = (200),
xlims = (-10,25),
linecolor=:transparent,
dpi = 300,
xlabel = ("XlogP"),
label = false,
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency")

plot(MW_REACH, XlogP_REACH, size = (800,600), dpi = 300)


#Prepare the data
REACH_fingerprints =Matrix(REACH_data[:,4:end])

#Load the optimized random forest model
rf_cl = joblib.load("optimized_random_forest_classifier_RepoRT.joblib")

#Predict the classes
REACH_mobility_classes = ScikitLearn.predict(rf_cl,REACH_fingerprints)

#Plot the data
very_mobile = sum(REACH_mobility_classes.=="Very mobile")
mobile = sum(REACH_mobility_classes.=="Mobile")
non_mobile = sum(REACH_mobility_classes.=="Non-mobile")

indices_very_mobile = findall(x-> x == "Very mobile", REACH_mobility_classes)
indices_mobile = findall(x-> x == "Mobile", REACH_mobility_classes)
indices_non_mobile = findall(x-> x == "Non-mobile", REACH_mobility_classes)

barplot = bar([very_mobile, mobile, non_mobile], yformatter = :plain, 
xticks = (1:3 ,["Very mobile", "Mobile", "Non-mobile"]),
color = [:red, :darkorange, :green], alpha = 0.3, label = false)
annotate!(barplot, 0.9, 5000, text("$(round((very_mobile/(length(REACH_mobility_classes))*100), digits = 1))%", :left, 15))
annotate!(barplot, 1.9, 10000, text("$(round((mobile/(length(REACH_mobility_classes))*100), digits = 1))%", :left, 15))
annotate!(barplot, 2.85, 25000, text("$(round((non_mobile/(length(REACH_mobility_classes))*100), digits = 1))%", :left, 15))

histogram(REACH_data[indices_non_mobile,2],  
linecolor=:transparent,
 xlims = (0,1500),
dpi = 300,
xlabel = ("Molecular Weight (Da)"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :green, alpha = 1,
label = "Non-mobile")

histogram!(REACH_data[indices_mobile,2],  
linecolor=:transparent,
 xlims = (0,1500),
dpi = 300,
xlabel = ("Molecular Weight (Da)"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :darkorange, alpha = 1,
label = "Mobile")

histogram!(REACH_data[indices_very_mobile,2],  
linecolor=:transparent,
 xlims = (0,1500),
dpi = 300,
xlabel = ("Molecular Weight (Da)"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :red, alpha = 1,
label = "Very mobile")

###############################################


histogram(REACH_data[indices_non_mobile,3],  
bins = (200),
xlims = (-10,25),
linecolor=:green,
dpi = 300,
xlabel = ("XlogP"),
label = "Non-mobile",
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :green, alpha = 1
)

histogram!(REACH_data[indices_mobile,3],  
bins = (200),
xlims = (-10,25),
linecolor=:darkorange,
dpi = 300,
xlabel = ("XlogP"),
label = "Mobile",
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :darkorange, alpha = 1
)

histogram!(REACH_data[indices_very_mobile,3],  
bins = (200),
xlims = (-10,25),
linecolor=:red,
dpi = 300,
xlabel = ("XlogP"),
label = "Very mobile",
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :red, alpha = 1
)