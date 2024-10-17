# # Step 1: Read the contents of the file
# file_content = read("new_output.txt", String)
# elements = split(file_content)
# numbers = parse.(Float64, elements)
# comma_separated = join(numbers, ",")
# write("comma_separated.txt", comma_separated)

a = "/home/towlej/.julia/dev/TJLF.jl/tjlf-ep/new_output.txt"
b ="/home/towlej/.julia/dev/TJLF.jl/new_output.csv"
function convert_txt_to_csv(input_file, output_file)
    # Open the input file for reading
    open(input_file, "r") do infile
        # Initialize variables to store the data
        label = ""
        data = Dict{String, Vector{Float64}}()

        # Read each line in the file
        for line in eachline(infile)
            # Trim whitespace from the line
            line = strip(line)

            # Check if the line starts with a label (e.g., omegaGAM, gammap)
            if startswith(line, "omegaGAM") || startswith(line, "gammap") || startswith(line, "gammaE")
                # Extract the label and the initial numbers
                parts = split(line)
                label = parts[1]  # The first part is the label
                values = parse.(Float64, parts[2:end])  # Convert the remaining parts to numbers
                data[label] = values  # Store the initial data
            else
                # Continue reading numbers for the current label
                if label != ""
                    # Convert the line into numbers and append to the current label data
                    append!(data[label], parse.(Float64, split(line)))
                end
            end
        end

        # Open the output file for writing
        open(output_file, "w") do outfile
            for (key, values) in data
                # Write the label
                write(outfile, "$key, ")

                # Write the values, joined by commas
                write(outfile, join(values, ", "))

                # Newline after each array
                write(outfile, "\n")
            end
        end
    end
end


convert_txt_to_csv(a,b)

using CSV, DataFrames

# Load the CSV data into a DataFrame
df = CSV.read("/home/towlej/.julia/dev/TJLF.jl/new_output.csv", DataFrame, header=false)

# Convert each line to an array
gammap = df[1, :][2:end]  # Skips the first element, which is the label
omegaGAM = df[2, :][2:end] 
gammaE = df[3, :][2:end]

# Convert to Julia arrays and print
gammap_array = Array{Float64}(gammap)
omegaGAM_array = Array{Float64}(omegaGAM)
gammaE_array = Array{Float64}(gammaE)

# Print the arrays
# println("gammap = ", gammap_array)
# println("omegaGAM = ", omegaGAM_array)
#println("gammaE = ", gammaE_array)