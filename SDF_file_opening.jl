using CSV, DataFrames, StatsPlots
function parse_text_file(filename)
    data = Dict{String, Any}()
    current_key = "trash"

   
    for line in eachline(filename)
        
        if startswith(line, ">")
            current_key = strip(split(line, '>')[2])
            continue
        end

        if isempty(line)
            continue
        end

        if startswith(line,"\$\$\$\$")
            current_key = "trash"
            continue
        end

        if haskey(data, current_key)
                push!(data[current_key], line)
            else
                data[current_key] = [line]
            end
        end
    
    return data
end
# Function to convert parsed data to CSV
function convert_to_csv(data, output_file)
    # Extract relevant information
    namess = get(data, "<NAME", [""])
    mws = get(data, "<MW", [""])
    precursor_mzs = get(data, "<PRECURSOR M/Z", [""])
    collision_energies = get(data, "<COLLISION ENERGY", [""])
    ion_modes = get(data, "<ION MODE", [""])
    instruments = get(data, "<INSTRUMENT", [""])
    num_peaks = get(data, "<NUM PEAKS", [""])
    mass_spectral_peaks = get(data, "<MASS SPECTRAL PEAKS", [[""]])  # Default value in case not found
    mass = Vector{Vector{String}}(undef, length(namess))
    mass[1] = mass_spectral_peaks[1:parse(Int, num_peaks[1])]

    current_index = parse(Int, num_peaks[1])  

    for i = 2:length(namess)
        start_index = current_index + 1
        current_index += parse(Int, num_peaks[i])  
        end_index = current_index
        mass[i] = mass_spectral_peaks[start_index:end_index]
    end
    mass[125]
    # Create a DataFrame
    df = DataFrame(
        NAME = namess,
        MW = mws,
        PRECURSOR_M_Z = precursor_mzs,
        COLLISION_ENERGY = collision_energies,
        ION_MODE = ion_modes,
        INSTRUMENT = instruments,
        NUM_PEAKS = num_peaks,
        MASS = mass
    )

    # Write DataFrame to CSV
    CSV.write(output_file, df)
end
function sdf_to_data(sdf_file)
    # Parse the text file
    data = parse_text_file(sdf_file)

    # Convert parsed data to CSV
    convert_to_csv(data, replace(sdf_file, ".txt" => ".csv"))

    df_results = CSV.read(replace(sdf_file, ".txt" => ".csv"), DataFrame)
    
    # If compound_name and collision_energy are not provided, return df_results only
    return df_results
end
function sdf_to_data(sdf_file, compound_name, collision_energy)
    # Parse the text file
    data = parse_text_file(sdf_file)

    # Convert parsed data to CSV
    convert_to_csv(data, replace(sdf_file, ".txt" => ".csv"))

    df_results = CSV.read(replace(sdf_file, ".txt" => ".csv"), DataFrame)
    
        index = findall(x -> x == compound_name, df_results.NAME)
        if isempty(index)
            throw("There are no compounds with this name or it has been written incorrectly")
        end
        
        found_index = nothing
        for i in index
            if df_results.COLLISION_ENERGY[i] == collision_energy
                found_index = i
                break
            end
        end
        
        if found_index === nothing
            throw("There are no spectra with this collision energy")
        end

        Mass_info = df_results[found_index,:MASS]
        
        output_string = replace(Mass_info, r"[\[\]\"]" => "")

        # Replace ", " with "\n"
        output_string = replace(output_string, ", " => "\n")
        
        mass = Float64[]
        intensity = Float64[]

        # Regular expression to match mass and intensity pairs
        regex = r"(\d+\.\d+)\s(\d+(?:\.\d+)?)"

        # Extract mass and intensity pairs using regular expression
        matches = eachmatch(regex, Mass_info)

        # Iterate over matches and extract mass and intensity
        for match in matches
            push!(mass, parse(Float64, match.captures[1]))
            push!(intensity, parse(Float64, match.captures[2]))
        end
        top_3 = sortperm(intensity, rev = true)
        spectrum = bar(mass, intensity, title = "$(results.NAME[found_index]); Collision Energy = $(results.COLLISION_ENERGY[found_index]); Polarity: $(results.ION_MODE[found_index])", xlabel = "m/z",
        ylabel = "intensity", dpi = 300, grid = false, legend = false, titlefontsize = 10,
        right_margin = 20Plots.mm)
        annotate!(mass[top_3[1]]+5, intensity[top_3][1], ("$(round(mass[top_3[1]], digits = 4))", :10, :left))
        annotate!(mass[top_3[2]]+5, intensity[top_3][2], ("$(round(mass[top_3[2]], digits = 4))", :10, :left))
        annotate!(mass[top_3[3]]+5, intensity[top_3][3], ("$(round(mass[top_3[3]], digits = 4))", :10, :left))
        display(spectrum)

        return df_results, mass, intensity, output_string
    
    
   
end

##########################################################################################################################
#Input can be either just the sdf file as a .txt only which will return all the compounds 
#in the file with the data provided in the sdf in a csv in the same file path or you can 
#add a compound name and collision energy to get the spectrum and raw data.

#Example of only generating the csv file with the results, make sure the csv file is close if you want to generate the same file again or it will error: "SystemError: opening file"
sdf_file = "C:\\Users\\uqthulle\\Downloads\\AntiMicro_LibView_20240207_150354 1.txt"
sdf_file = "C:\\Users\\uqthulle\\Downloads\\Clan_Lab_LibView_20240207_150319 1.txt"
results = sdf_to_data(sdf_file)

#Example of generating the results and also the spectrum  and raw data for a compound
compound_name = "3,4-Methylenedioxyamphetamine"
collision_energy = 55
results, mass, intensity, output_string = sdf_to_data(sdf_file, compound_name, collision_energy)

#This will show the raw data of the spectrum in the REPL in the MASBANK format
print(output_string)
    
