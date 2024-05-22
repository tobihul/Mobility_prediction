using MS_Import, SAFD, StatsPlots, CompCreate, ULSA

pathin = "C:\\Users\\uqthulle\\Documents\\Test data" 
filenames = ["2710 4pm_01_inj01_32.mzXML"]
mz_thresh = [0, 0] #Sets threshold for mz values
int_thresh = 500 #Remove most of the noise

GC.gc()

mz_val, mz_int, t0, t_end, m, path, msModel, msIonisation, msManufacturer,
polarity, Rt = import_files_MS1(pathin, filenames, mz_thresh, int_thresh)
FileName = m[1]

###############
# Parameters

path2comp="/path/to/component.csv"
weight_f=[1,1,1,1,1,1,1]
mode="POSITIVE" # or "NEGATIVE"
source="ESI" # or "EI" please follow the MassBank classification of sources
mm=pathof(ULSA)
path2aux=joinpath(mm[1:end-7],"MassBankJulia.jld")
DB=load(path2aux,"MassBankJulia") # loading the MassBank
AccuMass = 194.0789 # The accurate mass of the parent in neutral mode
mass_tol = 0.02 # The mass tellerance
ms2val = [138.068, 195.091, 110.073, 123.043] # The mz values of the fragments
ms2int = [3404.4, 2898.2, 873.0, 426.4] # The intensity of fragments


################
#

ids = featureID_comp(mode,source,path2comp,weight_f) # for identification

featureID_comp_batch(mode,source,path2comps,weight_f) # The identification in batch mode

mfs = featureMF_comp(mode,source,path2comp) # for molecular formula assignment

table = featureID_ext(DB,mode,source,AccuMass,weight_f,mass_tol,ms2val,ms2int) # for idenetification of indivdual features generated externally
