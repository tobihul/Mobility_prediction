using Statistics, PubChemCrawler, CSV, DataFrames, StatsPlots, Distributions,LinearAlgebra, PyCall, Conda, Distributed
import HTTP
const prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
pcp = pyimport("pubchempy")
pd = pyimport("padelpy")
#This function converts the chemical names to PubChem CID's and is implemeted in the next two functions
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

my_get_cid(name ="657-24-9")
Compounds = CSV.read("C:\\Users\\uqthulle\\Documents\\Schulze PMOC's.csv", DataFrame)
cids::Vector{Int32} = zeros(length(Compounds[:,1]))
for i = 1:length(Compounds[:,1])
    cids[i] = my_get_cid(name=Compounds[i,1])
    @show i
end
cidss::Vector{Int32} = trunc.(Int32, cids)
cidssf::Vector{Int32} = filter!(e->e≠0,cidss)

Df_LogP_MW = CSV.File(get_for_cids(cidssf; properties="MolecularWeight,XLogP", output="CSV")) |> DataFrame
Df_LogP_MW_no_missing = dropmissing(Df_LogP_MW)

scatter(Df_LogP_MW_no_missing[:,2],Df_LogP_MW_no_missing[:,3], xlims = (0,3500), ylims = (-35,75), palette = :default,
size = (1280, 720), grid = false, xlabel = "Molecular weight", ylabel = "XLogP3", 
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
legend = :topleft,
legendfont=font(13), markersize = 5, markerstrokewidth = 0.75, dpi = 300,xtickfont=font(13), ztickfont = font(13),
ytickfont=font(13), 
guidefont=font(20))


X = Df_LogP_MW_no_missing[:,2]
Y = Df_LogP_MW_no_missing[:,3]
scatter(X,Y)

Compounds = CSV.read("R:\\PHD2024TH-Q6813\\Literature review\\Detected compounds and methods.csv", DataFrame)
cids::Vector{Int32} = zeros(length(Compounds[:,3]))
for i = 1:length(Compounds[:,3])
    cids[i] = my_get_cid(name=Compounds[i,3])
    @show i
end
cidss::Vector{Int32} = trunc.(Int32, cids)
cidssf::Vector{Int32} = filter!(e->e≠0,cidss)
Df_LogP_MW_pubchem = CSV.File(get_for_cids(cidssf; properties="CanonicalSMILES,MolecularWeight,XLOGP", output="CSV")) |> DataFrame
Df_LogP_MW_no_missing = dropmissing(Df_LogP_MW_pubchem)
XLOGP3 = [Df_LogP_MW_no_missing[:,3] Df_LogP_MW_no_missing[:,4]]
smiles = Df_LogP_MW_pubchem[:,2]
smiles_no_missing = Df_LogP_MW_no_missing[:,2]

function smiles_to_logp(smiles::Vector{String})
    XLogPs::Vector{Float32} = zeros(length(smiles))
    for i = 1:length(smiles)
        fp::Dict{Any,Any} = pd.from_smiles(smiles[i], fingerprints = true, descriptors = false)
        XLogPs[i] = parse(Float64, (fp["XLogP"]))
        @show i
    end
    return XLogPs
end

pd.from_smiles(smiles[1], fingerprints = true, descriptors = true)


@time Xlogpss = smiles_to_logp(smiles)

Properties = [Df_LogP_MW_no_missing[:,3] XLogPs]

scatter(Properties[:,1], Properties[:,2], xlims = (0,3500), ylims = (-35,75), palette = :default,
size = (1280, 720), grid = false, xlabel = "Molecular weight", ylabel = "XLogP", 
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
legend = :topleft,
legendfont=font(13), markersize = 5, markerstrokewidth = 0.75, dpi = 300,xtickfont=font(13), ztickfont = font(13),
ytickfont=font(13), 
guidefont=font(20), label = "XlogP from PaDEL n = 64")

scatter!(XLOGP3[:,1], XLOGP3[:,2], xlims = (0,3500), ylims = (-35,75), palette = :default,
size = (1280, 720), grid = false, xlabel = "Molecular weight", ylabel = "XLogP", 
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
legend = :topleft,
legendfont=font(13), markersize = 5, markerstrokewidth = 0.75, dpi = 300,xtickfont=font(13), ztickfont = font(13),
ytickfont=font(13), 
guidefont=font(20), label = "XlogP from PubChem n = 49")



Df_LogP_MW_no_missing[:,3]

XLogPs

df_test=DataFrame(PubChemCID = Df_LogP_MW_no_missing[:,1], SMILES = Df_LogP_MW_no_missing[:,2], MW = Df_LogP_MW_no_missing[:,3], PubChem_XlogP3 = Df_LogP_MW_no_missing[:,4], PaDEL_XlogP = XLogPs)



scatter(XLOGP3[:,2], XLogPs, legend = false, dpi = 300, xlabel = "PubChem XlogP3", ylabel = "PaDEL XlogP", title = "PMOCs from Neuwald et al 2021")
cor(XLOGP3[:,2], XLogPs)




CSV.write("R:\\PHD2024TH-Q6813\\Literature review\\XlogP comparisons.csv", df_test)

Df_LogP_MW = CSV.File(get_for_cids(31389; properties="MolecularWeight,XLogP,CanonicalSMILES", output="CSV")) |> DataFrame


scatter(Df_LogP_MW_no_missing[:,2],Df_LogP_MW_no_missing[:,3], xlims = (0,3500), ylims = (-35,75), palette = :default,
size = (1280, 720), grid = false, xlabel = "Molecular weight", ylabel = "XLogP3", 
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
legend = :topleft,
legendfont=font(13), markersize = 5, markerstrokewidth = 0.75, dpi = 300,xtickfont=font(13), ztickfont = font(13),
ytickfont=font(13), 
guidefont=font(20))

savefig("C:\\Users\\uqthulle\\Downloads\\Pubchem vs PaDEL XlogP.png")

using QRCode

qrcode("https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/65ae77639138d231617c37ef/original/exploring-the-chemical-subspace-of-rplc-a-data-driven-approach.pdf")
exportqrcode("https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/65ae77639138d231617c37ef/original/exploring-the-chemical-subspace-of-rplc-a-data-driven-approach.pdf")
exportqrcode("https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/65ae77639138d231617c37ef/original/exploring-the-chemical-subspace-of-rplc-a-data-driven-approach.pdf", "C:\\Users\\uqthulle\\Downloads\\QR code Denise paper.png", Medium(), targetsize = 10, compact = true)