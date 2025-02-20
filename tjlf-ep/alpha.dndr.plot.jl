using Plots
gr()
# Function to read the data from a file where each value is on a separate line (after skipping the first line)
function read_column_data(filename::String)
    file = open(filename, "r")  # Open the file in read mode
    lines = readlines(file)    # Read all lines into an array of strings
    close(file)                # Close the file after reading

    # Ignore the first line, and read the remaining lines as individual numbers
    data_lines = lines[2:end]  # Skip the first line

    # Convert each line to a float and store it in an array
    data_array = parse.(Float64, data_lines)  # Convert each line into a Float64

    return data_array
end

# Function to read the data from a file where the numbers are comma-separated in a single line
function read_comma_separated_data(filename::String)
    file = open(filename, "r")  # Open the file in read mode
    lines = readlines(file)    # Read all lines into an array of strings
    close(file)                # Close the file after reading

    # Ignore the first line, and read the second line as a string
    data_str = lines[2]        # Get the second line (which is the array of numbers)

    # Remove any unwanted characters such as brackets
    clean_str = replace(data_str, r"[\[\]]" => "")  # Remove brackets

    # Split the string by commas and convert it into an array of floats
    data_array = parse.(Float64, split(clean_str, ','))  # Split by comma

    return data_array
end

# Read the two files into arrays
data1 = read_column_data("/home/towlej/.julia/dev/TJLF/tjlf-ep/alpha_dndr_crit.fortran.input")  # File 1: Column format
data2 = read_comma_separated_data("/home/towlej/.julia/dev/TJLF/tjlf-ep/alpha_dndr_crit.input")  # File 2: Comma-separated format

# Generate an x-axis scale from 2 to 202 with the same number of points as the data arrays
x_values = range(2, stop=202, length=length(data1))

# Check if both data arrays have the same length
if length(data1) != length(x_values) || length(data2) != length(x_values)
    error("Data arrays and x-values must have the same length!")
end

# Plot the two datasets on the same plot
plot(x_values, data1, label="Fortran", xlabel="Radius", ylabel="dn/dr", title="Fortran vs Julia dn/dr critical", 
     legend=:topright, linewidth=2, linestyle=:solid)  # Solid line for the first dataset

plot!(x_values, data2, label="Julia", xlabel="Radius", ylabel="dn/dr", legend=:topright, 
      linewidth=2, linestyle=:dash)  # Dashed line for the second dataset
savefig("plot.png")