using Mobility_prediction
using Test

#Test data 
@testset "smiles_to_mobility single" begin
    ##Test a single run already pre-compiled
    @test smiles_to_mobility("CN1C=NC2=C1C(=O)N(C(=O)N2C)C") == ("Mobile", 82)
    ##Test a single run that needs to be calculated
    @test smiles_to_mobility("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC") == (["Non-mobile"], 88)
end

SMILES_precomputed = ["CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CCO", "C1=CC(=C(C=C1CCN)O)O"]
SMILES_non_precomputed = ["CCCCCCCCCCCCCCCCCCCCCCCCCCCOOCCCCCCC", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "CCCCCCCCCCCCCCCCCCCCOCCCCCCCCCCCCCCC"]
SMILES_erroneous = ["CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "C56"]
#@testset "smiles_to_mobility batch run"

    path = mktempdir()

    #Precomputed test
    result_precomp = smiles_to_mobility(path, SMILES_precomputed) 
    @test result_precomp.Predicted_mobility  == ["Mobile", "Very mobile","Very mobile"]
    @test result_precomp.Probability  == [82,49,94]

    #Non precomputed test
    result_non_precomp = smiles_to_mobility(path, SMILES_non_precomputed) 
    @test result_non_precomp.Predicted_mobility  == ["Non-mobile", "Non-mobile","Non-mobile"]
    @test result_non_precomp.Probability  == [95,88,96]

    #Test with erroneous SMILES
    result_wrong_smiles = smiles_to_mobility(path, SMILES_erroneous)
    @test result_non_precomp.Predicted_mobility  == ["Mobile", "Non-mobile","Non-mobile"]
    @test result_non_precomp.Probability  == [95,88,96]