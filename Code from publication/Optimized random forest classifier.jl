include("All functions.jl")

#Random forest modelling after performing hyperparameter optimization using hyperparameter optimization.jl
filtered_RPLC_data = CSV.read("R:\\PHD2024TH-Q6813\\Models and other documents\\filtered_RPLC_data.csv", DataFrame)
fingerprints = CSV.read("R:\\PHD2024TH-Q6813\\Models and other documents\\RepoRT fingerprints final.csv", DataFrame)
fingerprints = Matrix(fingerprints)

#Splitting data into train and test by avoiding data leakage



X_train, X_test, y_train, y_test, train_indices, test_indices, y = train_test_split_no_leakage_classifier(filtered_RPLC_data, 0.1, fingerprints)

#Build the model with the optimal settings
rf_cl = RandomForestClassifier(n_estimators = 50, max_features = 0.2,
                                    min_samples_split = 2,
                                   min_samples_leaf = 4,
                                   random_state = 42, class_weight = "balanced", n_jobs = -1)


ScikitLearn.fit!(rf_cl, X_train, y_train)

#Saving the model
sklearn = pyimport("sklearn")
RandomForestClassifier = sklearn.ensemble.RandomForestClassifier
joblib = pyimport("joblib")
np = pyimport("numpy")

joblib.dump(rf_cl, "optimized_random_forest_classifier_RepoRT.joblib")

############################################
##If you just want to use the model already trained start hyperparameter
#Load in the optimized random forest classifier
rf_cl = joblib.load("optimized_random_forest_classifier_RepoRT.joblib")

#Depths
depths = [maximum([tree.tree_.max_depth for tree in rf_cl.estimators_]) for _ in 1:length(rf_cl.estimators_)]

#Prediction of train and test
y_hat_train = ScikitLearn.predict(rf_cl,X_train)
y_hat_test = ScikitLearn.predict(rf_cl,X_test)

#Scores of train and test
score_train = ScikitLearn.score(rf_cl, X_train, y_train)
score_test = ScikitLearn.score(rf_cl, X_test, y_test)

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
labels_plot_1 = join(split(labels[1])[2:end])
label_plot_rest = [split(labell)[2] for labell in labels[2:10]]
labels = [labels_plot_1; label_plot_rest]
bar(labels[1:10],sort(importance, rev = true)[1:10],
xrotation=30, 
dpi = 300,
left_margin = 10Plots.mm,
bottom_margin = 7.5Plots.mm,
legend = false,
ylabel = "Importance (%)",
legendfont = font(11), xtickfont=font(10),
guidefont=font(15), ytickfont = font(10))

df_importances = DataFrame(Fingerprint = labels, Importance = (sort(importance, rev = true)))
CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\List of PubChem FPs and importances in RF model.csv", df_importances)
#Tracking misclassifications

#This is used to take a look at specific cases of misclassifications. Change y_test and y_hat_test as needed with: "Very mobile, Mobile or Non-Mobile
indices = []
for i in eachindex(y_test)
    if y_test[i] == "Non-mobile" && y_hat_test[i] == "Non-mobile"
        push!(indices, i)
    end
end

#Indices of the compounds
cmp_indices = test_indices[indices]


#Their organic modifier
p_bs = filtered_RPLC_data.Modifier[cmp_indices][62]

histogram(filtered_RPLC_data[cmp_indices,6], xlims = (0,1))

histogram(p_bs./100, bins = 10, dpi = 300, xlims = (0,1),
ylabel = "Frequency",
xlabel = "Φ",
label = "Non-mobile misclassified as very mobile",
legend = :topright, legendfont = font(12), xtickfont=font(12), 
ytickfont=font(12), 
guidefont=font(18),
ylims = (0,50))

vline!([0.6], linestyle = :dash, label = "Non-mobile label threshold",
c = :red)


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

