using Statistics, CSV, DataFrames, PyCall, Conda, PubChemCrawler, DecisionTree, Random, ScikitLearn
import HTTP
const prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
pcp = pyimport("pubchempy")
pd = pyimport("padelpy")
function my_get_cid(; name=nothing, smiles=nothing, kwargs...)
    input = "compound/"
    name !== nothing && (input *= "name/$(HTTP.escapeuri(name))/")
    smiles !== nothing && (input *= "smiles/$((smiles))/")
    url = prolog * input * "cids/TXT"
    r = HTTP.request("GET", url; kwargs...)
    cids_string = String(r.body)
    cids = split(cids_string, "\n")
    cids = [cid for cid in cids if !isempty(cid) && !isspace(cid[1])]
    return parse(Int, cids[1])
end
data = CSV.read("C:\\Users\\uqthulle\\Downloads\\HILIC measured compounds with RI.csv", DataFrame)

compounds = data[:,1]

cids::Vector{Int32} = zeros(length(data[:,1]))
for i = 1:length(data[:,1])
    cids[i] = my_get_cid(name=data[i,1])
    @show i
end

properties = CSV.File(get_for_cids(cids; properties="CanonicalSMILES", output="CSV")) |> DataFrame
smiles = properties[:,2]
fps = pd.from_smiles(smiles, fingerprints = false)
set = DataFrame(SMILES = smiles)
desc_p = DataFrame(pd.from_smiles(set[i,"SMILES"],fingerprints=true, descriptors = false))

compounds_names = keys(fps)
compound_values = values(fps)

fpss = DataFrame(pd.from_smiles(set[1,"SMILES"],fingerprints=false))
Emptyrows = Vector{Int}()
for i = 2:size(set,1)
        fingerprints = DataFrame(pd.from_smiles(set[i,"SMILES"],fingerprints=false))
        if any(x -> any(col -> col == "", x), eachcol(fingerprints))
            push!(Emptyrows,i)
            continue
        else
            fpss = vcat(fpss,fingerprints) 
        end
end
X = Matrix(fpss)

vscodedisplay(fpss)
X = parse.(Float64,X)
X = X .- mean(X, dims=1)
X = X./std(X, dims = 1)
y = (data[:,2])
y = vec(y[setdiff(1:size(y, 1), Emptyrows), :])

using ScikitLearn.CrossValidation: train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

model = RandomForestRegressor(n_trees=100, min_samples_leaf=8)

fit!(model, X_train, y_train)

predict(model,X_test)

y_test