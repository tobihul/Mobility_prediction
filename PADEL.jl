using Statistics
using PyCall
using CSV
using DataFrames
using Conda
using ProgressBars
using RIprediction
pcp = pyimport("pubchempy")
pd = pyimport("padelpy")


################################################################

################################################################
# Buiding the functions
#
#
################################################################

function smiles2ds(smiles,fp::Int64=1)

    desc_p = []
    if fp == 1
        try
            desc_p = DataFrame(pd.from_smiles(smiles,fingerprints=true))
        catch
            println("Something has gone wrong with this entry!")
        end
    elseif fp == 0
        try
            desc_p = DataFrame(pd.from_smiles(smiles,fingerprints=false))
        catch
            println("Something has gone wrong with this entry!")
        end
    end

    try
        desc_p[!,"AATS0e"]

    catch
        println("Something has gone wrong with this entry!")
        desc_p = []

    end


    return desc_p


end


smiles = "CC1=C(C(=CC=C1)S(=O)(=O)O)C"

logp = pd.from_smiles(smiles, fingerprints = false)


##################################################################
# Collect pubchem info based on pcid numbers


function pubchem_qury_pcid(pcid)

    rep = Array{Any}(undef,length(pcid),7)

    for i=1:length(pcid)
        #println(i)
        if pcid[i] != "NULL"
            n = parse(Int64,pcid[i])
            c = pcp.Compound.from_cid(n)
            rep[i,1] = c.iupac_name
            rep[i,2] = c.inchikey
            rep[i,3] = c.inchi
            rep[i,4] = c.canonical_smiles
            rep[i,5] = c.molecular_formula
            rep[i,6] = c.monoisotopic_mass
            rep[i,7] = c.xlogp


        else
            rep[i,:] = zeros(1,7)

        end

    end

    rep_ = DataFrame(NAME = rep[:,1],INCHIKEY = rep[:,2],INCHI = rep[:,3],SMILES = rep[:,4], FORMULA = rep[:,5],
     MONOISOMASS = rep[:,6],XLOGP = rep[:,7])

    return rep_

end



# pcid = df[!,"PubChem_CID"][1:5]
# rep = pubchem_qury_pcid(pcid)

##################################################################
# Collect pubchem info based on smiles


function pubchem_qury_smiles(smiles)

    rep = Array{Any}(undef,length(smiles),7)

    for i=1:length(smiles)
        #println(i)
        if smiles[i] != "NULL"
            n = pcp.get_cids(smiles[i,1],"smiles")
            c = pcp.Compound.from_cid(n)
            rep[i,1] = c.iupac_name
            rep[i,2] = c.inchikey
            rep[i,3] = c.inchi
            rep[i,4] = c.canonical_smiles
            rep[i,5] = c.molecular_formula
            rep[i,6] = c.monoisotopic_mass
            rep[i,7] = c.xlogp

        else
            rep[i,:] = zeros(1,7)

        end

    end

    rep_ = DataFrame(NAME = rep[:,1],INCHIKEY = rep[:,2],INCHI = rep[:,3],SMILES = rep[:,4], FORMULA = rep[:,5],
     MONOISOMASS = rep[:,6],XLOGP = rep[:,7])

    return rep_

end


# pubchem_qury_smiles(smiles)


##################################################################
# Batch collecation of descriptors and PubChem info


function smiles2desc(smiles)

    pc_inf = pubchem_qury_smiles([smiles[1]])

    desc_p = hcat(pc_inf,smiles2ds(smiles[1]))

    for i in ProgressBar(2:size(smiles,1))

        try
            pc_inf_t = pubchem_qury_smiles([smiles[i]])
            desc_p_temp = hcat(pc_inf_t,smiles2ds(smiles[i]))
            append!(desc_p,desc_p_temp)
            #println(i)

        catch
            continue
        end
    end


    return desc_p

end

# desc_p = smiles2desc(smiles)


function smiles2fingerprints(smiles)
    #set active path to where the settings file is
    cd(joinpath(pathof(RIprediction)[1:end-19],"Models"))
    # cd("Models")
    set = DataFrame(SMILES = smiles)

    pubchemFPs = "PubchemFP" .* string.(collect(115:262))
    indStart = size(set,2)+1
    indFPs = []

    for i = 1:size(set,1)
        println(i)
        desc_p = []
        try
            desc_p = DataFrame(pd.from_smiles(set[i,"SMILES"],fingerprints=true, descriptors = false))
        catch
            continue
        end
        if isempty(desc_p)
            if i==1
                println("Error at first iteration")
                return
            end
            continue
        end
        #get FPs
        if i == 1
            #expand dataframe
            for f = 1:size(desc_p,2)
                if contains(names(desc_p)[f],"APC2") || any(names(desc_p)[f] .== pubchemFPs)
                        indFPs = [indFPs;f]
                        set[!,names(desc_p)[f]] = fill("",size(set,1))
                end
            end
        set[i,indStart:end] = desc_p[1,indFPs]
        else
        set[i,indStart:end] = desc_p[1,indFPs]
        end
    end
    return set
end


pubchem_qury_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

smiles = ["C(N(CP(=O)(O)O)CP(=O)(O)O)P(=O)(O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]
out = RIpredictionFP(smiles)