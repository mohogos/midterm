# Creating encryption/decryption key

function build(type::Type{EcoliEncryptionKey})::EcoliEncryptionKey

    model = EcoliEncryptionKey();

    keydict = Dict{Char, String}();
    keydict['A'] = "C6H12O6"
    keydict['B'] = "C5H10N2O3"
    keydict['C'] = "C5H8NO4"
    keydict['D'] = "C2H1O3"
    keydict['E'] = "H2O"
    keydict['F'] = "H"
    keydict['G'] = "C6H5O7"
    keydict['H'] = "C3H5O3"
    keydict['I'] = "C4H4O5"
    keydict['J'] = "C21H26N7O14P2"
    keydict['K'] = "C21H27N7O14P2"
    keydict['L'] = "C21H25N7O17P3"
    keydict['M'] = "C21H26N7O17P3"
    keydict['N'] = "H4N"
    keydict['O'] = "C3H4O10P2"
    keydict['P'] = "O2"
    keydict['Q'] = "C3H4O7P"
    keydict['R'] = "C4H2O5"
    keydict['S'] = "C3H2O6P"
    keydict['T'] = "C6H10O10P"
    keydict['U'] = "HO4P"
    keydict['V'] = "C6H9O9P"
    keydict['W'] = "HO4P"
    keydict['X'] = "C2H3O2"
    keydict['Y'] = "C3H3O3"
    keydict['a'] = "C49H74O4"
    keydict['b'] = "C49H76O4"
    keydict['c'] = "C5H9O8P"
    keydict['d'] = "C2H3O2"
    keydict['e'] = "C2H4O"
    keydict['f'] = "C7H13O10P"
    keydict['g'] = "C23H34N7O17P3S"
    keydict['h'] = "C4H4O4"
    keydict['i'] = "C25H35N7O19P3S"
    keydict['j'] = "C6H3O6"
    keydict['k'] = "C5H9O8P"
    keydict['l'] = "C2H3O5P"
    keydict['m'] = "C10H12N5O10P2"
    keydict['n'] = "C5H4O5"
    keydict['o'] = "C10H12N5O7P"
    keydict['p'] = "C10H12N5O13P3"
    keydict['q'] = "C6H5O7"
    keydict['r'] = "CO2"
    keydict['s'] = "C21H32N7O16P3S"
    keydict['t'] = "C3H5O6P"
    keydict['u'] = "C4H7O7P"
    keydict['v'] = "C2H6O"
    keydict['w'] = "C6H11O9P"
    keydict['x'] = "C6H10O12P2"
    keydict['y'] = "CH1O2"
    keydict['z'] = "C4H2O4"
    keydict['1'] = "C3H5O6P"
    keydict['2'] = "C6H11O9P"


# Building stoichiometric model

function _build_stoichiometric_matrix(data::Dict{String,Any})::Array{Float64,2}
    
    # initialize -
    list_of_metabolites = data["metabolites"];
    list_of_reactions = data["reactions"];
    number_of_reactions = length(list_of_reactions);
    number_of_metabolites = length(list_of_metabolites);
    S = Array{Float64,2}(undef, number_of_metabolites, number_of_reactions);
    
    for i ∈ 1:number_of_metabolites

        metabolite_dict = list_of_metabolites[i];
        key_id = metabolite_dict["id"];

        for j ∈ 1:number_of_reactions
            reaction_dict = list_of_reactions[i];
            reaction_metabolites_dict = reaction_dict["metabolites"];
            if (haskey(reaction_metabolites_dict, key_id) == true)
                st_value = reaction_metabolites_dict[key_id];
            else
                S[i,j] = 0.0;
            end
        end
    end

    # return -
    return S
end

"""
Fill me in.
"""
function _build_metabolite_id_array(data::Dict{String,Any})::Array{String,1}

    # initialize -
    metabolite_id_array = Array{String,1}()

    metabolites = data["metabolites"];
    for metabolite ∈ metabolites
        id_value = metabolite["id"];
        push!(metabolite_id_array, id_value);
    end

    # return -
    return metabolite_id_array;
end

"""
Fill me in.
"""
function _build_reaction_id_array(data::Dict{String,Any})::Array{String,1}
    
    # initialize -
    reaction_id_array = Array{String,1}()

    reactions = data["reactions"];
    for reaction ∈ reactions
        id_value = reaction["id"];
        push!(reaction_id_array, id_value);
    end

    # return -
    return reaction_id_array;
end

"""
Fill me in.
"""
function _build_bounds_array(data::Dict{String,Any})::Array{Float64,2}

    # initialize -
    list_of_reactions = data["reactions"];
    number_of_reactions = length(list_of_reactions)
    bounds_array = Array{Float64,2}(undef,number_of_reactions,2)

    lower_bound = 0.0
    upper_bound = 1000.0

    for i ∈ 1:number_of_reactions
        reaction = list_of_reactions[i];
        L = reaction["lower_bound"]
        U = reaction["upper_bound"]
        bounds_array[i,1]
    end
    # return -
    return bounds_array
end

"""
Fill me in
"""
function build(type::Type{MyStoichiometricNetworkModel}, data::Dict{String,Any})::MyStoichiometricNetworkModel

    # build an empty instance of our model -
    model = MyStoichiometricNetworkModel();

    # construct model elements -
    model.species = _build_metabolite_id_array(data);
    model.reactions = _build_reaction_id_array(data);
    model.bounds = _build_bounds_array(data);
    model.S = _build_stoichiometric_matrix(data);

    # return -
    return model;
end


