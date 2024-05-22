using CSV

folder_path = "C:\\Users\\uqthulle\\Downloads\\MassBank-data-2023.11\\MassBank-data-2023.11"

using DataFrames

function extract_info(filename::String)
    compound_name = ""
    smiles = ""
    column_name = "Unknown"

    for line in eachline(filename)
        if startswith(line, "CH\$NAME:")
            compound_name = split(strip(line), ": ")[2]
        elseif startswith(line, "CH\$SMILES:")
            smiles = split(strip(line), ": ")[2]
        elseif startswith(line, "AC\$CHROMATOGRAPHY: COLUMN_NAME")
            column_name = split(strip(line), "COLUMN_NAME ")[2]
        end
    end
    df = DataFrame(Compound_Name = compound_name, SMILES = smiles, Column_Name = column_name)
    return df
end
function process_folder(folder_path::String)
    filenames = readdir(folder_path)
    data = DataFrame(Compound_Name = String[], SMILES = String[], Column_Name = String[])
    for filename in filenames
        endswith(filename, ".txt")
            filepath = joinpath(folder_path, filename)
            info = extract_info(filepath)
            push!(data, info)
        
    end
    return data
end


@time extract_info(filename)

CSV.write("C:\\Users\\uqthulle\\Downloads\\data_massbank test.csv", data)

# Example usage:
file_path = "/path/to/your/text/file.txt"
df = parse_text_file(file_path)
