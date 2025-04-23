include("All functions.jl")
import ScikitLearn: RandomForestClassifier
RPLC_pH3_pH26_data = CSV.read("C:\\Users\\uqthulle\\Documents\\RPLC_data_pH3and2.6_updated.csv", DataFrame)
fingerprints = CSV.read("C:\\Users\\uqthulle\\Documents\\pH 3 and 2.6 RepoRT fingerprints final.csv", DataFrame)
fingerprints = Matrix(fingerprints)
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")

#X_train, X_test, y_train, y_test, train_indices, test_indices, y = train_test_split_no_leakage_classifier(RPLC_pH3_pH26_data, 0.1, fingerprints)


X_train = Matrix(CSV.read("C:\\Users\\uqthulle\\Documents\\X_train_final.csv", DataFrame))
X_test = Matrix(CSV.read("C:\\Users\\uqthulle\\Documents\\X_test_final.csv", DataFrame))
y_train = vec(CSV.read("C:\\Users\\uqthulle\\Documents\\y_train_final.csv", DataFrame)[:,1])
y_test = vec(CSV.read("C:\\Users\\uqthulle\\Documents\\y_test_final.csv", DataFrame)[:,1])
train_indices = vec(CSV.read("C:\\Users\\uqthulle\\Documents\\train_indices_final.csv", DataFrame)[:,1])
test_indices = vec(CSV.read("C:\\Users\\uqthulle\\Documents\\test_indices_final.csv", DataFrame)[:,1])



#Build the model with the optimal settings
rf_cl = RandomForestClassifier(n_estimators = 200, max_features = 0.5,
                                    min_samples_split = 2,
                                    max_depth = 100, 
                                   min_samples_leaf = 2,
                                   random_state = 42, class_weight = "balanced", n_jobs = -1)


ScikitLearn.fit!(rf_cl, X_train, y_train)

#Saving the model

RandomForestClassifier = sklearn.ensemble.RandomForestClassifier

joblib.dump(rf_cl, "optimized_random_forest_classifier_RepoRT_improved pH3 and 2.6 final.joblib", compress = 5)

############################################
##If you just want to use the model already trained start hyperparameter
#Load in the optimized random forest classifier
rf_cl = joblib.load("optimized_random_forest_classifier_RepoRT_improved pH3 and 2.6.joblib")

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
scatter([1,2,3], results[:,2], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "TPR", shape = :diamond, xlims = (0.5,3.5), grid = :y, markersize = 4,
title = "Test data n = $(length(y_test))", titlefont = font(10))

scatter!([1,2,3], results[:,3], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "FDR", shape = :square, markersize = 4)

p_test = scatter!([1,2,3], results[:,4], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "F1 score", shape = :utriangle, dpi = 300, markersize = 4)

plot(p_train, p_test, size = (800,400), dpi = 300)



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

df_importances = DataFrame(Fingerprint = labels, Importance = (sort(importance, rev = true)))
CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\List of PubChem FPs and importances in RF model.csv", df_importances)
CSV.write("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Paper plots updated\\Importances.csv", df_importances)
savefig("Plots\\Important variables.png")
#Tracking misclassifications

#This is used to take a look at specific cases of misclassifications. Change y_test and y_hat_test as needed with: "Very mobile, Mobile or Non-Mobile
indices = []
for i in eachindex(y_test)
    if y_test[i] == "Non-mobile" && y_hat_test[i] == "Mobile"
        push!(indices, i)
    end
end

#Indices of the compounds
cmp_indices = test_indices[indices]

probas = proba_test[indices]
mean(probas)
std(probas)
#Their organic modifier
p_bs = RPLC_pH3_pH26_data.Modifier[cmp_indices]



histogram_plot = histogram(p_bs ./ 100,
    xlabel = "Fraction of organic modifier / class probability",
    label = "Non-mobile misclassified as mobile",
    xlims = (0, 1),
    ylims = (0,300),
    legend = :none,  # Remove legend temporarily
    xtickfont = font(12),
    ytickfont = font(12),
    guidefont = font(12),
    left_margin = 5Plots.mm,
    right_margin = 5Plots.mm,
    dpi = 300,
    bins = 10,
    ylabel = "Frequency"
)

# Create a twin axis for the density plot
plot_with_twin = twinx()

# Add the density plot
density!(plot_with_twin, probas,
    linecolor = :red,
    linewidth = 2,
    alpha = 0.65,
    dpi =300,
    ylims = (0,4),
    ytickfont = font(12),
    label = "",  # Remove label here
    ylabel = "Normalized probability",
    legend = :none # Remove legend here too
)

# Add the vertical line
vline!([0.6],
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

# Add invisible plots to create legend entries

plot!(Float64[], Float64[], label = "Class probability distribution", seriestype = :line, color = :red, linewidth = 2, alpha = 0.65)
plot!(Float64[], Float64[], label = "Non-mobile threshold", seriestype = :line, color = :red, linestyle = :dash)
savefig("Plots//Non-mobile as mobile.png")
#Their InChI
Inchi = filtered_RPLC_data[cmp_indices,2][62]

cmps = findall(x-> x == Inchi, filtered_RPLC_data[:,2])
means = mean(filtered_RPLC_data.Modifier[cmps]./100)
stds = std(filtered_RPLC_data.Modifier[cmps]./100)

MW = filtered_RPLC_data.MW[cmps][1]
XlogP = filtered_RPLC_data.XlogP[cmps][1]

 scatter(filtered_RPLC_data.Modifier[cmps]./100)

#Applicability domain of the test set
lev = calculate_leverage(X_train, X_test)

histogram(lev, dpi = 300,
xlabel = "hᵢ", label = false, xlims = (0,0.2),
linecolor = :transparent)

warning_h = (3 * length(X_train[1,:])+1)/length(X_train[:,1])

vline!([warning_h], label = "warning (h*)")


######Permutation Testing

#Build the model with the optimal settings
rf_cl = RandomForestClassifier(n_estimators = 25, max_features = 0.5,
                                    min_samples_split = 2,
                                    max_depth = 100, 
                                   min_samples_leaf = 2,
                                   random_state = 42, class_weight = "balanced", n_jobs = -1)


@time ScikitLearn.fit!(rf_cl, X_train, y_train)

y_hat_test = ScikitLearn.predict(rf_cl,X_test)

c_matrix = confusion_matrix(y_hat_test, y_test)

results = TPR_FDR(c_matrix)

mean_F1_not_random = mean(results[:,4])

n_permutation = 100
permutation_F1 = zeros(n_permutation)

for i = 1:n_permutation
    @show i

    y_train_perm = shuffle(y_train)

    ScikitLearn.fit!(rf_cl, X_train, y_train_perm)

    y_hat_test_perm = ScikitLearn.predict(rf_cl,X_test)

    c_matrix = confusion_matrix(y_hat_test_perm, y_test)

    results = TPR_FDR(c_matrix)

    @show permutation_F1[i] = mean(results[:,4])

end

df_permutation_score = DataFrame(score = permutation_F1)
CSV.write("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Paper plots updated\\Permutation test_score.csv", df_permutation_score)
al_results = vcat(mean_F1_not_random, permutation_F1)

histogram(al_results, xlims = (0,1),
xlabel = "Mean F1 score all classes", 
ylabel = "permutations",
dpi = 300, label = false, ylims = (0,50))
vline!([mean_F1_not_random], linestyle = :dash, label = 
"Test data")
vline!([mean(permutation_F1)], c =:black, linestyle = :dash,
label = "Shuffled permutations")

m = mean(al_results)
s = std(al_results)

using HypothesisTests

t_result = OneSampleTTest(permutation_F1, mean_F1_not_random)

pval = pvalue(t_result)
savefig("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Paper plots updated\\Permutation test.png")


data = CSV.read("C:\\Users\\uqthulle\\Documents\\CV3 results new model", DataFrame)

vscodedisplay(data)