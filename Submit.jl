include("Include.jl")

# loading file
path_to_json_file = joinpath(_PATH_TO_DATA,"e_coli_core.json");
data = readreactionfile(path_to_json_file);
(L,X,list_of_measurements) = load_ecoli_data_file(path_to_data_file);

