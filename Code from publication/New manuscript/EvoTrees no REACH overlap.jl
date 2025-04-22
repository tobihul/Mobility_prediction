include("All functions.jl")
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")

#Random forest modelling after performing hyperparameter optimization using hyperparameter optimization.jl
filtered_RPLC_data = CSV.read("C:\\Users\\uqthulle\\Documents\\RPLC_no_REACH.csv", DataFrame)
filtered_RPLC_data_test = CSV.read("C:\\Users\\uqthulle\\Documents\\RPLC_with_REACH.csv", DataFrame)

fingerprints_train = CSV.read("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT no REACH overlap.csv", DataFrame)
fingerprints_train = Matrix(fingerprints_train)

fingerprints_test = CSV.read("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT with REACH overlap.csv", DataFrame)

fingerprints_test = Matrix(fingerprints_test)
#These are the fingerprint names for each corresponding column
#Splitting data into train and test by avoiding data leakage

shuffle_indices = shuffle(collect(1:length(filtered_RPLC_data[:,1])))
shuffle_indices_test = shuffle(collect(1:length(filtered_RPLC_data_test[:,1])))

x_train = fingerprints_train[shuffle_indices,:]
y_train = assign_labels(filtered_RPLC_data.Modifier./100)[shuffle_indices]

x_eval = fingerprints_test[shuffle_indices_test,:][1:7826,:]
y_eval = assign_labels(filtered_RPLC_data_test.Modifier./100)[shuffle_indices_test][1:7826,:][:]

x_test = fingerprints_test[shuffle_indices_test,:][7827:end,:]
y_test = assign_labels(filtered_RPLC_data_test.Modifier./100)[shuffle_indices_test][7827:end,:][:]


using EvoTrees: fit

config = EvoTreeClassifier(
    early_stopping_rounds = 5,
    nrounds = 200,
    eta = 0.1,
    L2 = 0.1,
    lambda = 0,
    gamma = 0,
    max_depth = 9,
    min_weight = 1,
    rowsample = 0.8,
    colsample = 1,
    nbins = 32,
    rng = 1,
    #device = :cpu

)
@time model = EvoTrees.fit_evotree(config; x_train, y_train, x_eval, y_eval, feature_names = PubChem_keys, print_every_n = 1)
plot(model,2)
#Prediction of train and test
y_hat_train_probas = EvoTrees.predict(model, x_train)
y_hat_eval_probas = EvoTrees.predict(model, x_eval)
y_hat_test_probas = EvoTrees.predict(model, x_test)

labels = ["Mobile", "Non-mobile", "Very mobile"]

y_hat_train = labels[argmax.(eachrow(y_hat_train_probas))]
y_hat_eval = labels[argmax.(eachrow(y_hat_eval_probas))]
y_hat_test = labels[argmax.(eachrow(y_hat_test_probas))]

#Scores of train and test
score_train = sum(y_hat_train .== y_train) / length(y_train)
score_eval = sum(y_hat_eval .== y_eval) / length(y_eval)
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