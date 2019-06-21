# --- PRIVATE METHODS ---------------------------------------------------------- #
function check_stoichiometric_matrix(original_dictionary,generated_dictionary)

    # get the STMs -
    old_stm = original_dictionary["S"]
    new_stm = generated_dictionary["S"]

    (NRO,NCO) = size(old_stm)
    (NRG,NCG) = size(new_stm)
    if (NRO != NRG || NCO != NCG)
        return false
    end

    # element by element comparision -
    for row_index = 1:NRO
        for col_index = 1:NCO

            if (old_stm[row_index,col_index] != new_stm[row_index, col_index])
                return false
            end
        end
    end

    return true
end

function check_metabolite_order(original_dictionary,generated_dictionary)

    # get mets -
    mets_old = original_dictionary["mets"]
    mets_new = generated_dictionary["mets"]

    if (length(mets_old) != length(mets_new))
        return false
    end

    N = length(mets_old)
    for index = 1:N

        msg = "comparing: $(mets_old[index]) == $(mets_new[index])"
        println(msg)

        if (mets_old[index] != mets_new[index])
            return false
        end
    end

    return true
end

function check_reaction_order(original_dictionary,generated_dictionary)

    # get rxns -
    rxns_old = original_dictionary["rxns"]
    rxns_new = generated_dictionary["rxns"]

    if (length(rxns_old) != length(rxns_new))
        return false
    end

    N = length(rxns_old)
    for index = 1:N

        msg = "comparing: $(rxns_old[index]) == $(rxns_new[index])"
        println(msg)

        if (rxns_old[index] != rxns_new[index])
            return false
        end
    end

    return true
end

function check_reaction_reversibility(original_dictionary,generated_dictionary)

    # get rxns -
    rev_old = original_dictionary["rev"]
    rev_new = generated_dictionary["rev"]

    if (length(rev_old) != length(rev_new))
        return false
    end

    N = length(rev_old)
    for index = 1:N

        msg = "comparing: $(rev_old[index]) == $(rev_new[index])"
        println(msg)

        if (rev_old[index] != rev_new[index])
            return false
        end
    end

    return true
end

# --- PUBLIC METHODS ----------------------------------------------------------- #
function is_generated_cobra_mat_equal(path_to_original_cobra_mat_file::String,path_to_generated_cobra_mat_file::String)::Bool

    # load the original file -
    original_model_name = splitdir(splitext(path_to_original_cobra_mat_file)[1])[2]
    original_dictionary = load_cobra_model_file(path_to_original_cobra_mat_file,original_model_name)

    # load generated file -
    generated_model_name = splitdir(splitext(path_to_generated_cobra_mat_file)[1])[2]
    generated_dictionary = load_cobra_model_file(path_to_generated_cobra_mat_file,generated_model_name)

    # checks -
    if (check_stoichiometric_matrix(original_dictionary,generated_dictionary) == false)
        msg = "Ooops! stoichiometric matrix check failed ..."
        println(msg)
        return false
    end

    if (check_metabolite_order(original_dictionary, generated_dictionary) == false)
        msg = "Ooops! metabolite order check failed ..."
        println(msg)
        return false
    end

    if (check_reaction_order(original_dictionary, generated_dictionary) == false)
        msg = "Ooops! reaction order check failed ..."
        println(msg)
        return false
    end

    if (check_reaction_reversibility(original_dictionary, generated_dictionary) == false)
        msg = "Ooops! reaction reversibility check failed ..."
        println(msg)
        return false
    end

    return true
end
