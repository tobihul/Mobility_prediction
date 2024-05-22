using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots, MultivariateStats

folder_path = "C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data"
Final_table = DataFrame(Inchikey = String[], LC_mode = String7[], pH_A = Float64[], pH_B = Float64[], Retention_factors = Float64[])
Descriptors = []
for i = 1:374
    #Locating the file
    current_file = joinpath("C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data", readdir(folder_path)[i])
    println(i)
    #Getting the SMILES and retention times
    rt_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "rtdata_canonical_success.tsv"],"_"))
    data_compounds = CSV.read(rt_path, DataFrame)
    Compound_Inchi = data_compounds[:,6]
    RT = data_compounds[:,4]

     #Getting the retention factors
     gradient_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "gradient.tsv"],"_"))
     gradient = CSV.read(gradient_path, DataFrame)
     local retention_factors
     try 
         run_time = gradient[end,1]
         retention_factors = RT./run_time
     catch
         retention_factors = repeat(["No gradient info"], length(RT))
     end

      #Getting the method pH 
      metadata_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "metadata.tsv"],"_"))
      metadata = CSV.read(metadata_path, DataFrame)
      pH_A = ones(length(RT)).*select(metadata, "eluent.A.pH")[1,1]
      pH_B = ones(length(RT)).*select(metadata, "eluent.B.pH")[1,1]
      
     #Getting the descriptors and adding the eluent A and B pH 
     #descriptors_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "descriptors_canonical_success.tsv"],"_"))
     #descriptors_info = CSV.read(descriptors_path, DataFrame)
     #df_merged = innerjoin(data_compounds, descriptors_info, on = :id)
     #selected_columns = String[]
     #for col in names(df_merged)
         # Check if the column has no missing values
         #if !any(ismissing, df_merged[!, col])
             # If the column has no missing values, add its name to the selected_columns vector
             #push!(selected_columns, col)
         #end
     #end
     #df_no_missing_current = select(df_merged, selected_columns)
   
 
     #Descriptors = vcat(Descriptors,df_merged)

    #Getting the LC mode info for the data set
    info_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "info.tsv"],"_"))
    info = CSV.read(info_path, DataFrame)
    #Check if the column info is missing in the meta data
    LC_mode = info[1,3]
    if ismissing(LC_mode) || LC_mode == "NA"

        LC_mode = repeat(["No mode info"], length(RT))

    else

        LC_mode = repeat([LC_mode], length(RT))

    end

    
  
    #Getting all the data in one table to append all the other datasets in this format
    current_table = DataFrame(Inchikey = Compound_Inchi, LC_mode = LC_mode, pH_A = pH_A, pH_B = pH_B, Retention_factors = retention_factors)
    Final_table = vcat(Final_table, current_table)
end

Final_table
Descriptors
Final_table[:,3][1]

histogram(Final_table.pH_B)


mean(Final_table.pH_A)