# Mobility_prediction.jl

**Mobility_prediction.jl** is a tool that can be used to obtain the predicted mobility of any chemical. For more than 200,000 unqiue SMILES, the fingerprints have already been pre-computed. If your query chemicals are included, the tool should be almost instant. If they are not it may take some time to calculate the fingerprints depending on your hardware.

To read more about how it works: https://pubs.acs.org/doi/10.1021/acs.est.5c07274

## Installation

In order to install the **LC_MS_Resolved_Peaks.jl** package in Julia, run the follwing: "]" to enter package manager and then "add https://github.com/tobihul/Mobility_prediction"

Alternatively: 

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/tobihul/Mobility_prediction"))

```

## Usage

There are several options to query chemicals

Firstly, a single SMILES can be run without saving it

```julia
using Mobility_prediction

SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

smiles_to_mobility(SMILES)
```

If you wish to save the result as a csv:

```julia
using Mobility_prediction

path = "path\\to\\save\\to"

SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

smiles_to_mobility(path, SMILES)
```
Finally, to run a list of SMILES and save it as a csv:

```julia
using Mobility_prediction

path = "path\\to\\save\\to"

SMILES = ["CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CC(=O)OC1C=CC2C3CC4=C5C2(C1OC5=C(C=C4)OC(=O)C)CCN3C"]

smiles_to_mobility(path, SMILES)
```
### Input

If running a single entry:

* **SMILES::String**

If running a single entry and saving as csv:

* **path::String**
* **SMILES::String**

If running a batch and saving as csv:

* **path::String**
* **SMILES::Vector{String}**

### Output

If running a single entry:

* **Tuple{String15, Int64}** Where the string is the mobility class and the Int is the class probability

If running a single entry and saving as csv

* **Tuple{String15, Int64}** Where the string is the mobility class and the Int is the class probability

If running a batch and saving as csv:

* **DataFrame** Where column 1 is the SMILES, column 2 is the predicted mobility class and column 3 is the probability of it being that class


## Notes

* **If you run a Batch file and it is saved as a csv and then run it again with different data but do not change the path, the file will be overwritten.**

* **If any of the SMILES in your list have spelling mistakes or are erroneous in any other way, the calculation will take significantly longer as the batch will run the errored region one-by-one.**

### Application to all RepoRT entries and REACH dataset chemicals n = 207,000+

The model was applied to all unique SMILES found both in the RepoRT dataset and all unique smiles in the REACH dataset. In the figure below the distribution of all chemicals in the MW and XlogP space can be seen grouped by predicted mobility class.

<p align="center">
  <img src="Precompiled FPS plot.png" alt="Plot description" width="600"/>
</p>

