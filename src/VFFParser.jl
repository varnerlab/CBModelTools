# --- PRIVATE METHODS ----------------------------------------------------------- #
function parse_vff_reaction_section(reaction_section_buffer::Array{String,1})

    # initialize -
    reaction_dictionary = Dict{String,CBMetabolicReaction}()
    reaction_order_array = Array{String,1}()

    # process each line -
    for reaction_line in reaction_section_buffer

        # skip comments and empty lines -
        if (occursin("//", reaction_line) == false && isempty(reaction_line) == false)

            # ok, record: name,{ec},LHS,RHS,R,F

            # initialize -
            reaction_wrapper = CBMetabolicReaction()

            # split around ,
            reaction_record_components_array = split(reaction_line,",")

            # get the name/key -
            reaction_key = reaction_record_components_array[1]

            # cache this (we need to keep the order that appears in the VFF) -
            push!(reaction_order_array, reaction_key)

            # populate -
            reaction_wrapper.record = reaction_line
            reaction_wrapper.reaction_name = reaction_key
            reaction_wrapper.ec_number = reaction_record_components_array[2]
            reaction_wrapper.left_phrase = reaction_record_components_array[3]
            reaction_wrapper.right_phrase = reaction_record_components_array[4]
            reaction_wrapper.reverse = reaction_record_components_array[5]
            reaction_wrapper.forward = reaction_record_components_array[6]

            # cache -
            reaction_dictionary[reaction_key] = reaction_wrapper
        end
    end

    # return -
    return (reaction_order_array, reaction_dictionary)
end

function parse_vff_rules_section(rules_section_buffer::Array{String,1})::Dict{String,CBRule}

    # initialize -
    rules_dictionary = Dict{String,CBRule}()

    # parse rule structure -
    for rule in rules_section_buffer

        # build new wrapper -
        rule_wrapper = CBRule()

        # split around =
        rule_components_array = split(rule,"=")

        # populate the wrapper -
        reaction_name = rule_components_array[1]
        rule_wrapper.reaction_name = reaction_name
        rule_wrapper.rule_text = rule_components_array[2]

        # cache -
        rules_dictionary[reaction_name] = rule_wrapper
    end

    # return -
    return rules_dictionary
end

function parse_vff_metabolite_section(metabolite_section_buffer::Array{String,1})::Dict{String,CBMetabolite}

    # initialize -
    metabolite_dictionary = Dict{String,CBMetabolite}()

    # iterate through each metabolite -
    for metabolite_record in metabolite_section_buffer

        # split around =
        metabolite_record_component_array = split(metabolite_record,"=")

        # metabolite symbol is [1]
        metabolite_symbol = metabolite_record_component_array[1]
        metabolite_name_field = metabolite_record_component_array[2]
        if (isempty(metabolite_name_field) == false)

            # split along :: -
            metabolite_name_symbol_array = split(metabolite_name_field,"::")
            full_metabolite_name = metabolite_name_symbol_array[1]
            kegg_id = metabolite_name_symbol_array[2]

            # build wrapper -
            metabolite_wrapper = CBMetabolite()
            metabolite_wrapper.symbol = metabolite_symbol
            metabolite_wrapper.name = full_metabolite_name
            metabolite_wrapper.kegg_id = kegg_id

            # cache -
            metabolite_dictionary[metabolite_symbol] = metabolite_wrapper
        end
    end

    # return -
    return metabolite_dictionary
end

function build_metabolite_symbol_array(reaction_array::Array{CBMetabolicReaction,1})

  # Method variables -
  species_symbol_array = Array{CBMetabolite,1}()

  # Helper function to parse the reaction phrases, split out the symbols
  function _parse_phrase(reaction_phrase::String)

    # Method variables -
    local_species_array = Array{CBMetabolite,1}()

    # Split around + -
    fragment_array = split(reaction_phrase,"+")
    for fragment in fragment_array

      if (contains(fragment,"*"))

          local_fragment_array = split(fragment,"*");
          species_symbol = CBMetabolite();
          species_symbol.symbol = local_fragment_array[2];

      else

        # Build -
        species_symbol = CBMetabolite()
        species_symbol.symbol = fragment
      end

      # check - is this the []?
      if (species_symbol.symbol != "[]")
          # grab -
          push!(local_species_array,species_symbol)
      end
    end

    # return -
    return local_species_array
  end

  function _isequal(species_model_1::CBMetabolite,species_model_2::CBMetabolite)
    if (species_model_1.symbol == species_model_2.symbol)
      return true
    end
    return false
  end

  function _add_symbol!(species_symbol_array,species_symbol)

    contains_species_already = false
    for cached_species_model in species_symbol_array

      if (_isequal(cached_species_model,species_symbol))
        contains_species_already = true
        break
      end
    end

    if (contains_species_already == false)
      push!(species_symbol_array,species_symbol)
    end
  end

  # iterate through and get the symbols -
  for reaction in reaction_array

    tmp_species_array_left = _parse_phrase(reaction.left_phrase)
    tmp_species_array_right = _parse_phrase(reaction.right_phrase)
    append!(tmp_species_array_left,tmp_species_array_right)

    for species_model in tmp_species_array_left
      _add_symbol!(species_symbol_array,species_model)
    end
  end

  # ok, so the species symbol array is *not* sorted)
  # let's sort the species -
  tmp_array = String[]
  for species_symbol::CBMetabolite in species_symbol_array
      push!(tmp_array,species_symbol.symbol)
  end

  # generate permutation array -
  idxa_sorted = sortperm(tmp_array)
  sorted_symbol_array = species_symbol_array[idxa_sorted]
  sorted_tmp_array = tmp_array[idxa_sorted]

  # ok, if a species contains [], then put it at the end -
  c_symbol_array = CBMetabolite[]
  e_symbol_array = CBMetabolite[]
  m_symbol_array = CBMetabolite[]
  r_symbol_array = CBMetabolite[]
  x_symbol_array = CBMetabolite[]
  t_symbol_array = CBMetabolite[]
  no_compartment_array = CBMetabolite[]
  partitioned_symbol_array = CBMetabolite[]

  # find all the species in [e]
  idx_all_e = findall(x->occursin("[e]",x)==true,sorted_tmp_array)
  for index in idx_all_e
      push!(partitioned_symbol_array,sorted_symbol_array[index])
  end

  # find all the species in [c]
  idx_all_c = findall(x->occursin("[c]",x)==true,sorted_tmp_array)
  for index in idx_all_c
      push!(partitioned_symbol_array,sorted_symbol_array[index])
  end

  # find all the species in [m]
  idx_all_m = findall(x->occursin("[m]",x)==true,sorted_tmp_array)
  for index in idx_all_m
      push!(partitioned_symbol_array,sorted_symbol_array[index])
  end

  # find all the species in [r]
  idx_all_r = findall(x->occursin("[r]",x)==true, sorted_tmp_array)
  for index in idx_all_r
      push!(partitioned_symbol_array,sorted_symbol_array[index])
  end

  # find all the species in [x]
  idx_all_x = findall(x->occursin("[x]",x)==true,sorted_tmp_array)
  for index in idx_all_x
      push!(partitioned_symbol_array,sorted_symbol_array[index])
  end

  # find all the species in [t]
  idx_all_t = findall(x->occursin("[t]",x)==true,sorted_tmp_array)
  for index in idx_all_t
      push!(partitioned_symbol_array,sorted_symbol_array[index])
  end

  # find all species w/no compartment -
  idx_all_nc = findall(x->occursin(r"\[*\]",x)==false,sorted_tmp_array)
  for index in idx_all_nc
      push!(partitioned_symbol_array,sorted_symbol_array[index])
  end

  # return -
  return partitioned_symbol_array
end

function build_stoichiometric_matrix(species_symbol_array::Array{CBMetabolite,1},reaction_array::Array{CBMetabolicReaction,1})

  # Method variables -
  number_of_species = length(species_symbol_array)
  number_of_reactions = length(reaction_array);
  stoichiometric_matrix = zeros(number_of_species,number_of_reactions);

  function _parse_reaction_phrase(lexeme,reaction_phrase)

    # Split around + -
    coefficient = 0.0;
    fragment_array = split(reaction_phrase,"+")
    for fragment in fragment_array

      if (contains(fragment,"*"))
          local_fragment_array = split(fragment,"*");
          test_lexeme = local_fragment_array[2];

          if (lexeme == test_lexeme)
            coefficient = parse(Float64,local_fragment_array[1]);
            break
          end

      else

        # Build -
        test_lexeme = fragment;
        if (lexeme == test_lexeme)
          coefficient = 1.0;
          break
        end
      end
    end

    return coefficient;
  end

  function _find_stoichiometric_coefficient(species_model::CBMetabolite,reaction::CBMetabolicReaction)

    # Method variables -
    stoichiometric_coefficient = 0.0

    # Check the left and right phrase -
    stoichiometric_coefficient += -1.0*(_parse_reaction_phrase(species_model.symbol,reaction.left_phrase))
    stoichiometric_coefficient += _parse_reaction_phrase(species_model.symbol,reaction.right_phrase)
    return stoichiometric_coefficient;
  end


  # setup counters -
  for (row_index,species_symbol) in enumerate(species_symbol_array)
    for (col_index,reaction) in enumerate(reaction_array)

      # Is this species involved in this reaction?
      stoichiometric_matrix[row_index,col_index] = _find_stoichiometric_coefficient(species_symbol,reaction);

    end
  end

  # return -
  return stoichiometric_matrix
end
# ------------------------------------------------------------------------------ #

# --- PUBLIC METHODS ----------------------------------------------------------- #
function load_vff_model_file(path_to_vff_file::String)::Dict{String,Any}

    # check the file path -
    is_file_path_ok(path_to_vff_file)

    # initialize -
    problem_dictionary = Dict{String,Any}()

    # load the file into an array -
    file_buffer_array = read_file_from_path(path_to_vff_file)

    # parse the reaction section -
    reaction_section = extract_section(file_buffer_array,"#REACTION::START","#REACTION::END")
    (reaction_order_array, reaction_dictionary) = parse_vff_reaction_section(reaction_section)
    problem_dictionary["list_of_reactions"] = reaction_dictionary
    problem_dictionary["reaction_order_array"] = reaction_order_array

    # parse the rules section -
    rule_section = extract_section(file_buffer_array,"#RULES::START","#RULES::END")
    problem_dictionary["list_of_rules"] = parse_vff_rules_section(rule_section)

    # parse the metabolite section -
    metabolite_section = extract_section(file_buffer_array,"#METABOLITES::START","#METABOLITES::END")
    problem_dictionary["list_of_metabolites"] = parse_vff_metabolite_section(metabolite_section)

    # return -
    return problem_dictionary
end

function generate_vff_from_cobra_dictionary(cobra_dictionary::Dict{String,Any}, path_to_vff_file::String; path_to_model_mapping_files::Union{String,Nothing}=nothing)

    # checks -
    is_dir_path_ok(path_to_vff_file)

    # initialize -
    vff_buffer = String[]

    # get the stoichiometric array, and some other stuff
    stoichiometric_matrix = Matrix(cobra_dictionary["S"])
    list_of_reaction_tags = cobra_dictionary["rxns"]
    list_of_species = cobra_dictionary["mets"]
    list_of_reversible_reactions = cobra_dictionary["rev"]
    list_of_rules = cobra_dictionary["rules"]
    list_of_species_names = cobra_dictionary["metNames"]

    # do we have a reaction_to_ec_mapping file?
    reaction_to_ec_mapping_dict = Dict{String,Mapping}()
    if (path_to_model_mapping_files != nothing)

        # create a path to the reaction ec mapping file -
        path_to_reaction_to_ec_mapping_file = joinpath(path_to_model_mapping_files,"ReactionECNumberMap.dat")

        # check path -
        is_file_path_ok(path_to_reaction_to_ec_mapping_file)

        # Ok, build the dictionary -
        reaction_to_ec_mapping_dict = build_mapping_dictionary(path_to_reaction_to_ec_mapping_file)
    end

    # what is the size?
    (number_of_species,number_of_reactions) = size(stoichiometric_matrix)

    # Declare a progress meter for user feedback -
    p = Progress(number_of_reactions,color=:yellow)

    # Add a header -
    filename = splitdir(path_to_vff_file)[2]
    +(vff_buffer,"// --------------------------------------------------------- //")
    +(vff_buffer,"// $(filename)")
    +(vff_buffer,"// GENERATED BY: CBModelTools")
    +(vff_buffer,"// GENERATED ON: $(Dates.now())")
    +(vff_buffer,"// SOURCE: https://github.com/varnerlab/CBModelTools")
    +(vff_buffer,"//")
    +(vff_buffer,"// VFF (Varner Flat File) format is a simple text based format")
    +(vff_buffer,"// for encoding metabolic (or other) types of biological models.")
    +(vff_buffer,"// Metabolic VFFs have three sections:")
    +(vff_buffer,"// RULES:")
    +(vff_buffer,"// METABOLITES:")
    +(vff_buffer,"// REACTIONS:")
    +(vff_buffer,"// --------------------------------------------------------- //")
    +(vff_buffer,"")

    # rules section -
    +(vff_buffer,"#RULES::START ---------------------------------------------- //")
    for (index,rule) in enumerate(list_of_rules)

        # init empty line -
        line = ""

        # get tag and reaction text -
        reaction_tag = list_of_reaction_tags[index]
        rule_text = list_of_rules[index]
        if (isempty(rule_text) == true)
            # build the line -
            line *= "$(reaction_tag)=[]"
        else
            # build the line -
            line *= "$(reaction_tag)=$(rule_text)"
        end

        # add rule to the buffer -
        +(vff_buffer,line)
    end

    +(vff_buffer,"#RULES::END ------------------------------------------------ //")
    +(vff_buffer,"")

    # Write the metabolites section -
    +(vff_buffer,"#METABOLITES::START ---------------------------------------- //")

    # build a kegg-id metabolite name map -
    mmt = build_metabolite_matching_table()
    for (index,metabolite) in enumerate(list_of_species)

        # init empty line -
        line = ""

        # get the metabolite name -
        metabolite_name = list_of_species_names[index]

        # lookup KEGG ID -
        kegg_id = "[]"
        if (haskey(mmt,metabolite_name) == true)
            kegg_id = mmt[metabolite_name]
        end

        # write -
        line*="$(metabolite)=$(metabolite_name)::$(kegg_id)"

        # add metabolite to the buffer -
        +(vff_buffer,line)
    end

    +(vff_buffer,"#METABOLITES::END ------------------------------------------ //")
    +(vff_buffer,"")

    # write reactions -
    +(vff_buffer,"#REACTION::START ------------------------------------------- //")
    for reaction_index = 1:number_of_reactions

        # initialize empty buffer -
        line = ""

        # get the reaction tag -
        reaction_tag_string = list_of_reaction_tags[reaction_index]

        # user message -
        msg = "Completed $(reaction_tag_string) ($(reaction_index) of $(number_of_reactions))"

        # update the progress bar -
        ProgressMeter.next!(p; showvalues = [(:status,msg)])

        # add the tag to the buffer -
        line *= "$(reaction_tag_string),"

        # ok, so I need to check, do we have is reaction key in my ec mapping?
        if (haskey(reaction_to_ec_mapping_dict,reaction_tag_string) == true)

            # we have a reaction tag = ec number record
            mapping_record = reaction_to_ec_mapping_dict[reaction_tag_string]

            # get the set of values -
            value_set = mapping_record.value
            for ec_number in value_set
                line*="$(ec_number)::"
            end

            # add [], -
            line*="[],"
        else
            # We have no ecnumber information -
            line *= "[],"
        end

        # find the reactants -
        idx_reactants = findall(stoichiometric_matrix[:,reaction_index].<0.0)
        if (isempty(idx_reactants) == true)
            line *= "[],"
        else

            # how many species do we have?
            number_of_species = length(idx_reactants)
            counter = 1
            for index in idx_reactants

                # get the metabolite -
                metabolite_string = list_of_species[index]
                stcoeff = stoichiometric_matrix[index,reaction_index]

                if (stcoeff != -1.0)
                    # add data to the buffer -
                    line *= "$(abs(stcoeff))*$(metabolite_string)"
                else
                    # add data to the buffer -
                    line *= "$(metabolite_string)"
                end

                # do we have more?
                if (counter < number_of_species)
                    line *= "+"
                else
                    line *= ","
                end

                counter = counter + 1
            end
        end

        # find the products -
        idx_products = findall(stoichiometric_matrix[:,reaction_index].>0.0)
        if (isempty(idx_products) == true)
            line *= "[],"
        else

            # how many species do we have?
            number_of_species = length(idx_products)
            counter = 1
            for index in idx_products

                # get the metabolite -
                metabolite_string = list_of_species[index]
                stcoeff = stoichiometric_matrix[index,reaction_index]

                if (stcoeff != 1.0)
                    # add data to the buffer -
                    line *= "$(stcoeff)*$(metabolite_string)"
                else
                    # add data to the buffer -
                    line *= "$(metabolite_string)"
                end

                # do we have more?
                if (counter < number_of_species)
                    line *= "+"
                else
                    line *= ","
                end

                counter = counter + 1
            end
        end

        # is this reaction reversible?
        rev_value = list_of_reversible_reactions[reaction_index]
        if (rev_value == 1.0)
            line *= "-inf,inf"
        else
            line *= "0,inf"
        end

        # add buffer to list of strings -
        +(vff_buffer,line)
    end
    +(vff_buffer,"#REACTION::END --------------------------------------------- //")

    # Write out the vff file -
    open("$(path_to_vff_file)", "w") do f

        for line_item in vff_buffer
            write(f,"$(line_item)\n")
        end
    end
end
# ------------------------------------------------------------------------------ #
