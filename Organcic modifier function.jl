using CSV, Statistics, DataFrames, StatsPlots, LinearAlgebra

folder_path = "C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data"
Final_table = DataFrame(Inchikey = String[], LC_mode = String7[], MACCS = String[], Pubchem_fps = String[], XlogP = Float64[], Retention_factors = Float64[])

#Locating the file
current_file = joinpath("C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data", readdir(folder_path)[1])
#Getting the retention factors
gradient_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][1])", "gradient.tsv"],"_"))
gradient = CSV.read(gradient_path, DataFrame)

timee = retention_factors[15]
gradient
run_time = gradient[end,1]

function interpolate_B_modifier(time::Float64, gradient::DataFrame)
    times = gradient[:,1]./run_time
    idx = searchsortedlast(times, time)
    if idx == 0
        return times[1]
    elseif idx == length(times)
        return gradient[end,3]
    elseif time > 1
        return gradient[end,3]
    else
        t1, t2 = times[idx], times[idx+1]
        B1, B2 = gradient[:,3][idx], gradient[:,3][idx+1]
        return B1 + (B2 - B1) * (time - t1) / (t2 - t1)
    end
end


