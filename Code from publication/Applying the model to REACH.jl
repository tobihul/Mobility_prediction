include("All functions.jl")

###Load in the REACH data
REACH_data = CSV.read("R:\\PHD2024TH-Q6813\\Models and other documents\\REACH fingerprints final.csv", DataFrame)

unique(REACH_data.smiles)
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
rf_class = joblib.load("optimized_random_forest_classifier_RepoRT.joblib")

#Predict the classes
REACH_mobility_classes = ScikitLearn.predict(rf_class,REACH_fingerprints)
REACH_probabilities = maximum(ScikitLearn.predict_proba(rf_class,REACH_fingerprints), dims = 2)
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

p_non_mobile = histogram(REACH_data[indices_non_mobile,2],  
linecolor=:transparent,
 xlims = (0,1500),
dpi = 300,
xlabel = ("Molecular Weight (Da)"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :green, alpha = 1,
label = "Non-mobile", linealpha = 0.3, legend = false)

p_mobile = histogram(REACH_data[indices_mobile,2],  
linecolor=:transparent,
 xlims = (0,1500),
dpi = 300,
xlabel = ("Molecular Weight (Da)"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :darkorange, alpha = 1,
label = "Mobile", linealpha = 0.3, legend = false)

p_very_mobile = histogram(REACH_data[indices_very_mobile,2],  
linecolor=:transparent,
 xlims = (0,1500),
dpi = 300,
xlabel = ("Molecular Weight (Da)"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :red, alpha = 1,
label = "Very mobile", linealpha = 0.3, legend = false)

l = @layout [
    a{0.7w} [
        grid(3, 1)
        

      
        
    ]
]
plot(p_all, p_non_mobile, p_mobile, p_very_mobile, layout=l, size = (1200,800),
bottom_margin = 5Plots.mm, right_margin = 5Plots.mm, left_margin = 5Plots.mm, dpi = 300)
###############################################
p_non_mobile = histogram(REACH_data[indices_non_mobile,3],  
linecolor=:transparent,
bins = 150,
 xlims = (-10,25),
dpi = 300,
xlabel = ("XlogP"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :green, alpha = 1,
label = "Non-mobile", linealpha = 0.3, legend = false)

p_mobile = histogram(REACH_data[indices_mobile,3],
bins = 150,
 xlims = (-10,25),  
linecolor=:transparent,
dpi = 300,
xlabel = ("XlogP"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :darkorange, alpha = 1,
label = "Mobile", linealpha = 0.3, legend = false)

p_very_mobile = histogram(REACH_data[indices_very_mobile,3],  
linecolor=:transparent,
bins = 150,
 xlims = (-10,25), 
dpi = 300,
xlabel = ("XlogP"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :red, alpha = 1,
label = "Very mobile", linealpha = 0.3, legend = false)

l = @layout [
    a{0.7w} [
        grid(3, 1)
        

      
        
    ]
]
plot(p_all, p_non_mobile, p_mobile, p_very_mobile, layout=l, size = (1200,800),
bottom_margin = 5Plots.mm, right_margin = 5Plots.mm, left_margin = 5Plots.mm, dpi = 300)
#Probability distribution
density(REACH_probabilities[REACH_mobility_classes.=="Non-mobile"],
 xlims = (0,1.1), c = :green, alpha = 1, linecolor = :green,
 label = "Non-mobile", xlabel = "Class probability", ylabel = "Probability density",
 legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11))

density!(REACH_probabilities[REACH_mobility_classes.=="Mobile"],
 c = :darkorange, alpha = 1, linecolor = :darkorange,
 label = "Mobile")

density!(REACH_probabilities[REACH_mobility_classes.=="Very mobile"]
, dpi = 300, c = :red, alpha = 1, label = "Very mobile", linecolor = :red,
 )



lev_REACH = calculate_leverage(X_train,REACH_fingerprints)

histogram(lev_REACH, dpi = 300,
xlabel = "hᵢ", label = false,
linecolor = :transparent, xlims =(0,0.6))
maximum(lev_REACH)
warning_h = (3 * length(X_train[1,:])+1)/length(X_train[:,1])
top95 = quantile(lev_REACH, 0.98)
vline!([warning_h], label = "warning (h*)")


#############################
#REACH OPERA data
OPERA_pred_data = CSV.read("R:\\PHD2024TH-Q6813\\Models and other documents\\OPERA predicted LogKoc and LogD.csv", DataFrame)

#Get the molecule ID to account for the missing ones
OPERA_mols = parse.(Int,[split(strings, "_")[2] for strings in OPERA_pred_data.MoleculeID])
OPERA_logKoc = OPERA_pred_data.LogKoc_pred
#Assign the mobility to the non missing entries

mobility_OPERA = []
for i in eachindex(REACH_mobility_classes)
    if i in OPERA_mols
        push!(mobility_OPERA, REACH_mobility_classes[i])
    else
        println(i)
    end
end

mobility_OPERA
indices_very_mobile = findall(x-> x == "Very mobile", mobility_OPERA)
indices_mobile = findall(x-> x == "Mobile", mobility_OPERA)
indices_non_mobile = findall(x-> x == "Non-mobile", mobility_OPERA)



##Plot logKoc

p_non_mobile = histogram(OPERA_logKoc[indices_non_mobile],  
linecolor=:transparent,
dpi = 300,
xlims = (0,7),
xlabel = ("Predicted logKoc"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :green, alpha = 1,
label = "Non-mobile", linealpha = 0.3, legend = true)

p_mobile = histogram!(OPERA_logKoc[indices_mobile],  
linecolor=:transparent,
dpi = 300,
xlims = (0,7),
xlabel = ("Predicted logKoc"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :darkorange, alpha = 1,
label = "Mobile", linealpha = 0.3, legend = true)

p_all = histogram!(OPERA_logKoc[indices_very_mobile],  
linecolor=:transparent,
dpi = 300,
xlims = (0,7),
xlabel = ("Predicted logKoc"),
legendfont = font(8), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :red, alpha = 1,
label = "Very mobile", linealpha = 0.3, legend = true)

vline!([3], c = :red, label = "EU very mobile criteria")
vline!([4], c = :darkorange, label = "EU mobile criteria")
l = @layout [
    a{0.7w} [
        grid(3, 1)
        

      
        
    ]
]
plot(p_all, p_non_mobile, p_mobile, p_very_mobile, layout=l, size = (1200,800),
bottom_margin = 5Plots.mm, right_margin = 5Plots.mm, left_margin = 5Plots.mm, dpi = 300)

(OPERA_pred_data.Conf_index_Koc, label = "log Koc Conf_index", dpi = 300)
density!(OPERA_pred_data.AD_index_Koc, label = "log Koc AD_index ", dpi = 300)

###Doing the same but for logD 5
OPERA_logD_55 = OPERA_pred_data.LogD55_pred

p_non_mobile = histogram(OPERA_logD_55[indices_non_mobile],  
linecolor=:transparent,
dpi = 300,
xlims = (-12,12),
xlabel = ("Predicted logD 5.5"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :green, alpha = 1,
label = "Non-mobile", linealpha = 0.3, legend = false)

p_mobile = histogram(OPERA_logD_55[indices_mobile],  
linecolor=:transparent,
dpi = 300,
xlims = (-12,12),
xlabel = ("Predicted logD 5.5"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :darkorange, alpha = 1,
label = "Mobile", linealpha = 0.3, legend = false)

p_very_mobile = histogram(OPERA_logD_55[indices_very_mobile],  
linecolor=:transparent,
dpi = 300,
xlims = (-12,12),
xlabel = ("Predicted logD 5.5"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :red, alpha = 1,
label = "Very mobile", linealpha = 0.3, legend = false)

l = @layout [
    a{0.7w} [
        grid(3, 1)
        

      
        
    ]
]
plot(p_all, p_non_mobile, p_mobile, p_very_mobile, layout=l, size = (1200,800),
bottom_margin = 5Plots.mm, right_margin = 5Plots.mm, left_margin = 5Plots.mm, dpi = 300)
density(OPERA_pred_data.Conf_index_LogD, label = "log D Conf_index", dpi = 300)
density!(OPERA_pred_data.AD_index_LogD, label = "log D AD_index ", dpi = 300)

####Log D 7.4
OPERA_logD_74 = OPERA_pred_data.LogD74_pred

p_non_mobile = histogram(OPERA_logD_74[indices_non_mobile],  
linecolor=:transparent,
dpi = 300,
xlims = (-15,12),
xlabel = ("Predicted logD 7.4"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :green, alpha = 1,
label = "Non-mobile", linealpha = 0.3, legend = false)

p_mobile = histogram(OPERA_logD_74[indices_mobile],  
linecolor=:transparent,
dpi = 300,
xlims = (-15,12),
xlabel = ("Predicted logD 7.4"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :darkorange, alpha = 1,
label = "Mobile", linealpha = 0.3, legend = false)

p_very_mobile = histogram(OPERA_logD_74[indices_very_mobile],  
linecolor=:transparent,
dpi = 300,
xlims = (-15,12),
xlabel = ("Predicted logD 7.4"),
legendfont = font(11), xtickfont=font(11), 
ytickfont=font(11), 
guidefont=font(11), ylabel = "Frequency", c = :red, alpha = 1,
label = "Very mobile", linealpha = 0.3, legend = false)

l = @layout [
    a{0.7w} [
        grid(3, 1)
        

      
        
    ]
]
plot(p_all, p_non_mobile, p_mobile, p_very_mobile, layout=l, size = (1200,800),
bottom_margin = 5Plots.mm, right_margin = 5Plots.mm, left_margin = 5Plots.mm, dpi = 300)

####Checking diffferent classifications of mobility both Models

OPERA_logkc_class_result = []

for logkoc in OPERA_logKoc
    if logkoc <= 2 
        push!(OPERA_logkc_class_result, "Very mobile")
    elseif logkoc <= 3 && logkoc > 2
        push!(OPERA_logkc_class_result, "Mobile")
    elseif isnan(logkoc)
        push!(OPERA_logkc_class_result, NaN)
    else
        push!(OPERA_logkc_class_result, "Non-mobile")
    end
end

sum(OPERA_logkc_class_result.=="Very mobile")/length(OPERA_logkc_class_result)*100
sum(OPERA_logkc_class_result.=="Mobile")/length(OPERA_logkc_class_result)*100
sum(OPERA_logkc_class_result.=="Non-mobile")/length(OPERA_logkc_class_result)*100


indices = []
for i in eachindex(mobility_OPERA)
    if OPERA_logkc_class_result[i] == "Non-mobile" && mobility_OPERA[i] == "Mobile"
        push!(indices, i)
    end
end
indices
count_RF_model = sum(mobility_OPERA.=="Mobile")

percentages = length(indices)/count_RF_model*100
indices

real_indices = vec(OPERA_mols[indices,:])

REACH_data[real_indices,1][20]