include("All functions.jl")
import ScikitLearn: RandomForestClassifier
merged_data = CSV.read("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_expanded_RPLC_pH3_data.csv", DataFrame)
fingerprints = CSV.read("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_expanded_entries_FPS.csv", DataFrame)
fingerprints = Matrix(fingerprints)
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")


X_train, X_test, y_train, y_test, train_indices, test_indices, y = train_test_split_no_leakage_classifier(merged_data, 0.2, fingerprints)

#Build the model with the optimal settings
rf_cl = RandomForestClassifier(n_estimators = 50, max_features = 0.5,
                                    min_samples_split = 2,
                                    max_depth = 100, 
                                   min_samples_leaf = 4,
                                   random_state = 42, class_weight = "balanced", n_jobs = -1,
                                   criterion = "log_loss", verbose = 1)


@time ScikitLearn.fit!(rf_cl, X_train, y_train)

#Saving the model

RandomForestClassifier = sklearn.ensemble.RandomForestClassifier

joblib.dump(rf_cl, "optimized_random_forest_classifier_RepoRT_improved and merged pH3 and 2.6 final.joblib", compress = 5)

############################################
##If you just want to use the model already trained start hyperparameter
#Load in the optimized random forest classifier
rf_cl = joblib.load("optimized_random_forest_classifier_RepoRT_improved pH3 and 2.6 final.joblib")

#Depths
depths = [maximum([tree.tree_.max_depth for tree in rf_cl.estimators_]) for _ in 1:length(rf_cl.estimators_)]

#Prediction of train and test
y_hat_train = ScikitLearn.predict(rf_cl,X_train)
y_hat_test = ScikitLearn.predict(rf_cl,X_test)

#Scores of train and test
score_train = ScikitLearn.score(rf_cl, X_train, y_train)
score_test = ScikitLearn.score(rf_cl, X_test, y_test)
proba_train= maximum(ScikitLearn.predict_proba(rf_cl,X_train), dims = 2)
proba_test = maximum(ScikitLearn.predict_proba(rf_cl,X_test), dims = 2)
#Confusion matrix for the train set
c_matrix = confusion_matrix(y_hat_train, y_train)

#TPR, FDR and F1-score for all classes for the train set
results = TPR_FDR(c_matrix)

scatter([1,2,3], results[:,2], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "TPR", shape = :diamond, xlims = (0.5,3.5), grid = :y, markersize = 4,
title = "Training data n = $(length(y_train))", titlefont = font(10),
legend = false)

scatter!([1,2,3], results[:,3], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "FDR", shape = :square, markersize = 4)

p_train = scatter!([1,2,3], results[:,4], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "F1 score", shape = :utriangle, dpi = 300, markersize = 4)

#Confusion matrix for the test set
c_matrix = confusion_matrix(y_hat_test, y_test)

#TPR, FDR and F1-score for all classes for the test set
results = TPR_FDR(c_matrix)
mean(results[:,4])
mean(results[:,3])
mean(results[:,2])
mean
scatter([1,2,3], results[:,2], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "TPR", shape = :diamond, xlims = (0.5,3.5), grid = :y, markersize = 4,
title = "Test data n = $(length(y_test))", titlefont = font(10))

scatter!([1,2,3], results[:,3], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "FDR", shape = :square, markersize = 4)

p_test = scatter!([1,2,3], results[:,4], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "F1 score", shape = :utriangle, dpi = 300, markersize = 4)

plot(p_train, p_test, size = (800,400), dpi = 300, legend = :left)

savefig("R:\\PHD2024TH-Q6813\\Research files\\Plots\\TPRFDRF1 scores plot.png")

#Checking feature importances
importance = rf_cl.feature_importances_.*100
sum(importance)
sorted_importance = sortperm(importance, rev = true)
labels = PubChem_keys[sorted_importance]
labels = labels = [join(split(labell)[2:end], " ") for labell in labels]



bar(labels[1:10],sort(importance, rev = true)[1:10],
xrotation=30, 
dpi = 300,
left_margin = 10Plots.mm,
bottom_margin = 13Plots.mm,
legend = false,
ylabel = "Importance (%)",
legendfont = font(11), xtickfont=font(8),
guidefont=font(15), ytickfont = font(10))

sum(sort(importance,rev = true)[1:10])
df_importances = DataFrame(Fingerprint = labels, Importance = (sort(importance, rev = true)))
CSV.write("R:\\PHD2024TH-Q6813\\Research files\\Models and other documents\\List of PubChem FPs and importances in RF model.csv", df_importances)
CSV.write("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Paper plots updated\\Importances.csv", df_importances)
savefig("R:\\PHD2024TH-Q6813\\Research files\\Models and other documents\\Important variables.png")
#Tracking misclassifications

#This is used to take a look at specific cases of misclassifications. Change y_test and y_hat_test as needed with: "Very mobile, Mobile or Non-Mobile
indices = []
for i in eachindex(y_test)
    if y_test[i] == "Very mobile" && y_hat_test[i] == "Non-mobile"
        push!(indices, i)
    end
end

# Indices of the compounds in the test set that were misclassified
cmp_indices = test_indices[indices]

#Search a random entry
search = rand(collect(1:length(cmp_indices)))
#search = findall(x-> x <=0, skipmissing(merged_data[cmp_indices,:].XlogP))[2]
sum(merged_data[cmp_indices, :][:,:].Modifier .<= 60)

misclassified_samples = merged_data[cmp_indices, :][search,:].InChi

SMILES = merged_data[cmp_indices, :][search,:].SMILES
probas = proba_test[indices][search]

#Their organic modifier
cmps = findall(x-> x == misclassified_samples,  merged_data.InChi)
mod = merged_data.Modifier[cmps]
sum(mod.<=60)
means = mean(merged_data.Modifier[cmps]./100)
stds = std(merged_data.Modifier[cmps]./100)

MW = merged_data.MW[cmps][1]
XlogP = merged_data.XlogP[cmps][1]

using StatsBase
mean(skipmissing(merged_data.XlogP[cmp_indices]))

studys = misclassified_samples = merged_data[cmp_indices, :][search,:].study
misclassified_inchi_all = merged_data[cmp_indices, :].study

counts = countmap(misclassified_inchi_all)
counts = countmap(merged_data[cmp_indices, :].XlogP)

histogram(merged_data[cmp_indices, :].XlogP)


InChI_search = "InChI=1S/C17H33NO4/c1-5-6-7-8-9-10-11-12-17(21)22-15(13-16(19)20)14-18(2,3)4/h15H,5-14H2,1-4H3"

misclassified_InChIs = misclassified_samples = merged_data[cmp_indices, :].InChi
indexxx = findall(x-> x == InChI_search, misclassified_InChIs)















indices = []
for i in eachindex(y_test)
    if y_test[i] == "Mobile" && y_hat_test[i] == "Mobile"
        push!(indices, i)
    end
end

# Indices of the compounds in the test set that were misclassified
cmp_indices = test_indices[indices]
probas = (proba_test[indices])
p_bs = merged_data.Modifier[cmp_indices]

histogram_plot = histogram(p_bs ./ 100,
    xlabel = "Fraction of organic modifier / class probability",
    xlims = (0, 1),
    ylims = (0,2000),
    label = "Correctly classified as mobile",
    legend = :none,  # Remove legend temporarily
    xtickfont = font(12),
    ytickfont = font(12),
    guidefont = font(12),
    left_margin = 5Plots.mm,
    right_margin = 5Plots.mm,
    dpi = 300,
    bins = 10,
    ylabel = "Frequency",
    formatter = :plain
)

plot_with_twin = twinx()

# Add the density plot
density!(plot_with_twin, probas,
    linecolor = :red,
    linewidth = 2,
    alpha = 0.65,
    dpi =300,
    ylims = (0,6),
    ytickfont = font(12),
    label = "",  # Remove label here
    ylabel = "Normalized probability",
    legend = :none # Remove legend here too
)

# Add the vertical line
vline!([0.2,0.6],
    linestyle = :dash,
    label = "",  # Remove label here
    color = :red
)

# Manually add all legend entries
plot!(histogram_plot,
    legend = :topleft,
    legendfont = font(12),
    foreground_color_legend = nothing,
    background_color_legend = :white,
    series_annotations = [""],
    dpi = 300
)
plot!(Float64[], Float64[], label = "Class probability distribution", seriestype = :line, color = :red, linewidth = 2, alpha = 0.65)
plot!(Float64[], Float64[], label = "Mobile label threshold", seriestype = :line, color = :red, linestyle = :dash)

savefig("R:\\PHD2024TH-Q6813\\Research files\\Models and other documents\\Plots\\Mobile correct.png")