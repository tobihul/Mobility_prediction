include("All functions.jl")
import ScikitLearn: RandomForestClassifier
merged_data = CSV.read("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_pH3_RPLC_data.csv", DataFrame)
fingerprints = CSV.read("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_expanded_entries_FPS.csv", DataFrame)
fingerprints = Matrix(fingerprints)
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")

all_indices = collect(1:nrow(merged_data))

# Then split *indices*, not the data directly
train_indices, test_indices = train_test_split(all_indices,
    test_size=0.2,
    stratify=merged_data.Class,
    shuffle=true,
    random_state=42
)

# Now index everything using those:
X_train = fingerprints[train_indices, :]
X_test = fingerprints[test_indices, :]
y_train = merged_data.Class[train_indices]
y_test = merged_data.Class[test_indices]





sum(y_train .== "Very mobile")
sum(y_train .== "Mobile")
sum(y_train .== "Non-mobile")

sum(y_test .== "Very mobile")
sum(y_test .== "Mobile")
sum(y_test .== "Non-mobile")
#Build the model with the optimal settings
rf_cl = RandomForestClassifier(n_estimators = 200, max_features = 0.5,
                                    min_samples_split = 2,
                                    max_depth = 100, 
                                   min_samples_leaf = 2,
                                   random_state = 42, class_weight = "balanced", n_jobs = -1)


ScikitLearn.fit!(rf_cl, X_train, y_train)

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
left_margin = 55Plots.mm,
bottom_margin = 20Plots.mm,
legend = false,
ylabel = "Importance (%)",
legendfont = font(11), xtickfont=font(8),
guidefont=font(15), ytickfont = font(10))

indices = []
for i in eachindex(y_test)
    if y_test[i] == "Non-mobile" && y_hat_test[i] == "Non-mobile"
        push!(indices, i)
    end
end

# Indices of the compounds in the test set that were misclassified
cmp_indices_non_mobile = test_indices[indices]
comp_indices_mobile = test_indices[indices]
comp_indices_very_mobile = test_indices[indices]

scatter(merged_data.MW[cmp_indices_non_mobile], merged_data.XlogP[cmp_indices_non_mobile], label = "Non-mobile test",
dpi = 300, xlabel = "MW", ylabel = "XlogP")
scatter!(merged_data.MW[cmp_indices_mobile], merged_data.XlogP[cmp_indices_mobile], label = "Mobile test")
scatter!(merged_data.MW[comp_indices_very_mobile], merged_data.XlogP[comp_indices_very_mobile], label = "Very mobile test")
savefig("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\MW xlogp plot test classes.png")
# Get the corresponding rows from the original filtered_RPLC_data
search = rand(collect(1:length(cmp_indices)))


misclassified_samples = merged_data[cmp_indices, :][search,:].SMILES
SMILES = merged_data[cmp_indices, :][search,:].SMILES
probas = mean(proba_test[indices])

#Their organic modifier
cmps = findall(x-> x == misclassified_samples,  merged_data.SMILES)
mod = merged_data.Modifier[cmps]
sum(mod.<=60)
means = mean(merged_data.Modifier[cmps]./100)
stds = std(merged_data.Modifier[cmps]./100)

MW = merged_data.MW[cmps][1]
XlogP = merged_data.XlogP[cmps][1]