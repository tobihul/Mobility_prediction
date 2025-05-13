include("All functions.jl")
using EvoTrees
#Random forest modelling after performing hyperparameter optimization using hyperparameter optimization.jl
RPLC_pH3_pH26_data = CSV.read("C:\\Users\\uqthulle\\Documents\\RPLC_data_pH3and2.6_updated.csv", DataFrame)
fingerprints = CSV.read("C:\\Users\\uqthulle\\Documents\\pH 3 and 2.6 RepoRT fingerprints final.csv", DataFrame)
fingerprints = Matrix(fingerprints)
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")

#Splitting data into train and test by avoiding data leakage

X_train, X_test, y_train, y_test, train_indices, test_indices, y = train_test_split_no_leakage_classifier(RPLC_pH3_data, 0.1, fingerprints)

X_train, X_eval, y_train, y_eval, train_indices, eval_indices, y = train_test_split_no_leakage_classifier(RPLC_pH3_pH26_data[train_indices,:], 0.111, fingerprints[train_indices,:])

x_train = X_train 

#Build the model with the optimal settings using EvoTrees
using EvoTrees: fit

#category_map = Dict("Mobile" => 1, "Non-mobile" => 2, "Very mobile" => 3)
#y_train = [category_map[x] for x in y_train]
#cat_y_test = [category_map[x] for x in y_test]

config = EvoTreeClassifier(
    early_stopping_rounds = 10,
    nrounds = 50,
    eta = 0.1,
    L2 = 0.1,
    lambda = 0.1,
    gamma = 0.1,
    max_depth = 5,
    min_weight = 1,
    rowsample = 1,
    colsample = 1,
    nbins = 64,
    rng = 1,
    device = :cpu

)
x_train = X_train
@time model = EvoTrees.fit_evotree(config; x_train, y_train, feature_names = PubChem_keys)
plot(model,2)
#Prediction of train and test
y_hat_train_probas = EvoTrees.predict(model, x_train)
y_hat_test_probas = EvoTrees.predict(model, X_test)

labels = ["Mobile", "Non-mobile", "Very mobile"]

y_hat_train = labels[argmax.(eachrow(y_hat_train_probas))]
y_hat_test = labels[argmax.(eachrow(y_hat_test_probas))]

#Scores of train and test
score_train = sum(y_hat_train .== y_train) / length(y_train)

score_test = sum(y_hat_test .== y_test) / length(y_test)

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
importance = EvoTrees.importance(model)
importance_matrix = hcat(string.(first.(importance)), last.(importance))

sorted_features = [parse(Int, v[2]) for v in split.(importance_matrix[:,1], "_")]
feature_contributions = importance_matrix[:,2].*100

labels = PubChem_keys[sorted_features]

bar(labels[1:10],feature_contributions[1:10],
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
    if y_test[i] == "Mobile" && y_hat_test[i] == "Very mobile"
        push!(indices, i)
    end
end

#Indices of the compounds
cmp_indices = test_indices[indices]

probas = proba_test[indices]
mean()
std(probas)
#Their organic modifier
p_bs = filtered_RPLC_data.Modifier[cmp_indices]


histogram_plot = histogram(p_bs ./ 100,
    xlabel = "Fraction of organic modifier / class probability",
    xlims = (0, 1),
    ylims = (0,250),
    label = "Non-mobile misclassified as mobile",
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
    ylims = (0,7),
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
plot!(Float64[], Float64[], label = "Non-mobile label threshold", seriestype = :line, color = :red, linestyle = :dash)
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
