using Oscar
input_file = "./fanlike seeds"  # Change to your actual file path
output_folder = "./Oscar_fanlike_seeds"

open(input_file, "r") do f
    for (i, line) in enumerate(eachline(f))
        # Parse the line into a list of lists of integers
        # Example line: "[[1,2],[2,3],[3,4]]"
        faces = eval(Meta.parse(line))
        K = simplicial_complex(faces)
        dim_K = dim(K)
        save("$(output_folder)/K_$(i-1)_$(dim_K)", K)
    end
end

input_file = "./minimally non-fanlike seeds"  # Change to your actual file path
output_folder = "./Oscar_minimally_non-fanlike_seeds"

open(input_file, "r") do f
    for (i, line) in enumerate(eachline(f))
        # Parse the line into a list of lists of integers
        # Example line: "[[1,2],[2,3],[3,4]]"
        faces = eval(Meta.parse(line))
        K = simplicial_complex(faces)
        dim_K = dim(K)
        save("$(output_folder)/L_$(i-1)_$(dim_K)", K)
    end
end