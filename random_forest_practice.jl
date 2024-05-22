using CSV, DataFrames, Random, StatsPlots

data_raw = CSV.read("C:\\Users\\uqthulle\\Documents\\LCvsGC dataset.csv", DataFrame)

y = data_raw.Response

X = Matrix(data_raw[:,10:end])

function perclass_splits(y, percent)
    uniq_class = unique(y)
    keep_index = []
    for class in uniq_class
        class_index = findall(y .== class)
        row_index = randsubseq(class_index, percent)
        push!(keep_index, row_index...)
    end
    return keep_index
end


train_index = perclass_splits(y, 0.85)

test_index = setdiff(1:length(y), train_index)

# split features

X_train = X[train_index, :]

X_test = X[test_index, :]

# split classes

y_train = y[train_index]

y_test = y[test_index]

#Running the model: Decision tree

model = DecisionTreeClassifier(max_depth = 100)


DecisionTree.fit!(model, X_train, y_train)

print_tree(model)

#view training data
train = [X_train y_train]

vscodedisplay(train)

#Prediciting

y_hat = DecisionTree.predict(model, X_test)

#Check the accuracy

accuracy = mean(y_hat .==y_test)

scatter(y_test, y_hat, xlabel = "y_test", ylabel = "y_hat", title = "Accuracy = $(round(accuracy,digits = 2))%", legend = false)

#Confusion matrix

DecisionTree.confusion_matrix(y_test, y_hat)

#Running the model: Random forest

model  = RandomForestClassifier(n_trees = 10)

DecisionTree.fit!(model, X_train, y_train)

## Prediciton

y_hat = DecisionTree.predict(model,X_test)

accuracy = mean(y_hat .==y_test)

DecisionTree.confusion_matrix(y_test, y_hat)

#ROC curve

FPR = 

TPR = 






#Optimizing number of n_trees

n_trees = collect(1:100)
scores = zeros(100)
for i in eachindex(n_trees)

    train_index = perclass_splits(y, 0.85)

    test_index = setdiff(1:length(y), train_index)

    X_train = X[train_index, :]

    X_test = X[test_index, :]

    y_train = y[train_index]

    y_test = y[test_index]

    model  = RandomForestClassifier(n_trees = trees_test[i])

    DecisionTree.fit!(model, X_train, y_train)

    y_hat = DecisionTree.predict(model,X_test)

    scores[i] = mean(y_hat .==y_test)
    @show i
end


model  = RandomForestClassifier(n_trees =10, min_samples_leaf = 3, n_subfeatures = 10)
DecisionTree.fit!(model, X_train, y_train)

y_hat = DecisionTree.predict(model,X_test)

    scores[i] = mean(y_hat .==y_test)
scatter(trees_test, scores)
vline!([10])


X[:,1445]
y

all = scatter(X[:,1436], X[:,1445], group = y, grid = false, xlabel = "MW", ylabel = "logP", title = "Compounds by analysis technique", dpi = 300)
indices_LC_GC = findall(x->x =="GC;LC", y)
indices_LC = findall(x->x =="LC", y)
indices_GC = findall(x->x =="GC", y)

LC_GC = scatter(X[indices_LC_GC,1436], X[indices_LC_GC,1445])
LC = scatter(X[indices_LC,1436], X[indices_LC,1445])
GC = scatter(X[indices_GC,1436], X[indices_GC,1445])

savefig("C:\\Users\\uqthulle\\Downloads\\LC GC dataset.png")



accuracy = nfoldCV_forest(y, X,
                          3,
                          n_subfeatures,
                          n_trees,
                          partial_sampling,
                          max_depth,
                          min_samples_leaf,
                          min_samples_split,
                          min_purity_increase;
                          verbose = true,
                          rng = seed)

train_index = perclass_splits(y, 0.85)

test_index = setdiff(1:length(y), train_index)

X_train = X[train_index, :]

X_test = X[test_index, :]

y_train = y[train_index]

y_test = y[test_index]

model  = RandomForestClassifier(n_trees =10, min_samples_leaf = 10, n_subfeatures = round(sqrt(length(y))))

DecisionTree.fit!(model, X_train, y_train)

y_hat = DecisionTree.predict(model,X_test)

score = mean(y_hat .==y_test)

# Define your grid of hyperparameters
n_trees_range = [10, 50, 100]
min_samples_leaf_range = [5, 10, 15]

best_score = -1
best_params = (n_trees = -1, min_samples_leaf = -1)

# Perform grid search
for n_trees in n_trees_range
    for min_samples_leaf in min_samples_leaf_range
        scores = Float64[]
        for _ in 1:5  # 5-fold cross-validation
            Random.seed!(42)  # Ensure reproducibility
            train_index = perclass_splits(y, 0.85)
            test_index = setdiff(1:length(y), train_index)

            X_train = X[train_index, :]
            X_test = X[test_index, :]
            y_train = y[train_index]
            y_test = y[test_index]

            model = RandomForestClassifier(n_trees = n_trees, min_samples_leaf = min_samples_leaf, n_subfeatures = round(sqrt(length(y_train))))
            DecisionTree.fit!(model, X_train, y_train)

            y_hat = DecisionTree.predict(model, X_test)
            score = mean(y_hat .== y_test)
            push!(scores, score)
        end
        avg_score = mean(scores)
        println("Params: n_trees=$n_trees, min_samples_leaf=$min_samples_leaf, Avg Score: $avg_score")
        if avg_score > best_score
            best_score = avg_score
            best_params = (n_trees = n_trees, min_samples_leaf = min_samples_leaf)
        end
    end
end

println("Best parameters: ", best_params)
println("Best score: ", best_score)


#
using ScikitLearn
using ScikitLearn.CrossValidation: cross_val_score
using ScikitLearn: @sk_import
@sk_import ensemble: RandomForestRegressor
@sk_import model_selection: train_test_split
using ScikitLearn.GridSearch: GridSearchCV
X = [MACCS Pubchem_fps]
X = X[indices_HILIC,:]
y = Final_table_unique[:,5][indices_HILIC]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

rf_regressor = RandomForestRegressor()

fit!(rf_regressor, X_train, y_train)

y_pred = rf_regressor.predict(X_test)
scatter(y_test,y_pred)

# Assuming y_true are the true target values and y_pred are the predicted values

# Calculate the mean of the true target values
y_mean = mean(y_test)

# Calculate the total sum of squares (TSS)
tss = sum((y_test .- mean(y_test)).^2)

# Calculate the residual sum of squares (RSS)
rss = sum((y_test .- y_pred).^2)

# Calculate the R-squared value
r_squared = 1 - rss / tss

# Print the R-squared value
println("R-squared: ", r_squared)

# Define the range of parameters to search
n_estimators = collect(100:100:600)
        min_samples_leaf = collect(2:2:8)
        max_features = [1.0, "sqrt", "log2"]
        param_grid = Dict("n_estimators" => n_estimators, "min_samples_leaf" => min_samples_leaf, "max_features" => max_features)

# Create a random forest regressor
rf_regressor = RandomForestRegressor()

# Perform grid search with cross-validation
grid_search = GridSearchCV(rf_regressor, param_grid, cv = 3)

# Fit the grid search to the data
fit!(grid_search, X, y)

# Get the best parameters and best score
best_params = grid_search.best_params_
best_score = grid_search.best_score_

# Print the best parameters and best score
println("Best Parameters: ", best_params)
println("Best Score: ", best_score)

predict(rf_regresso)