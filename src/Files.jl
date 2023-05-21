# build data matrix and label vectors for plot
function load_ecoli_data_file(path::String)::Tuple{Array{Int,1}, Array{Float64,2}, Array{String,1}}

    # load the data file -
    df = json.read(path,DataFrame);
    list_of_compounds = compounds(df);

    # grab the labels, and build the label vector
    Y = Int.(Vector(df[:,:visitid]));

    # grab the data, and build the data matrix -
    X = Matrix(df[:,3:end]) |> standardize

    # return -
    return (Y, X, list_of_names[3:end])
end

# read reaction file for stoichiometric matrix
function readreactionfile(path::String)::Dict{String,Any}
    return JSON.parsefile(path)
end