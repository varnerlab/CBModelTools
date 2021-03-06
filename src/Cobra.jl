# --- PRIVATE METHODS ---------------------------------------------------------- #
function is_kegg_organism_code_ok(kegg_organism_code::String)
end

function is_cobra_dictionary_ok(cobra_dictionary::Dict{String,Any})
end

function lookup_gene_ec_mapping_record(kegg_gene_id::String)::Union{Set{String},Nothing}

    # initialize -
    record_set = Set{String}()

    # call the KEGG API to get the ec number -
    list_of_ec_numbers = get_ec_number_for_gene(kegg_gene_id)
    if list_of_ec_numbers != nothing

        # create a record, add to the mapping buffer -
        for ec_number in list_of_ec_numbers
            push!(record_set,ec_number)
        end

        # return -
        return record_set
    end

    return nothing
end
# ------------------------------------------------------------------------------ #

# --- PUBLIC METHODS ----------------------------------------------------------- #
function load_cobra_model_file(path_to_cobra_mat_file::String, model_name::String)::Dict{String,Any}

    # check file ok?
    is_file_path_ok(path_to_cobra_mat_file)

    # open the cobra file -
    file = matopen(path_to_cobra_mat_file)
    cobra_dictionary = read(file,model_name)
    close(file)

    # return -
    return cobra_dictionary
end

function extract_files_from_cobra_model(cobra_dictionary::Dict{String,Any}, kegg_organism_code::String, path_to_model_directory::String)

    # check - does this path exist?
    is_dir_path_ok(path_to_model_directory)

    # check the cobra dictionary - is it legit?
    is_cobra_dictionary_ok(cobra_dictionary)

    # ok, looks ok if we get here - write out a bunch of files that we'll need

    # reaction mapping -
    path_to_reaction_mapping_file = joinpath(path_to_model_directory,"ReactionNameMap.dat")
    export_reaction_mapping_file(cobra_dictionary,path_to_reaction_mapping_file)

    # gene order -
    path_to_gene_order_file = joinpath(path_to_model_directory,"GeneOrder.dat")
    export_gene_order_file(cobra_dictionary, kegg_organism_code, path_to_gene_order_file)

    # reaction tags to gene mapping -
    path_reaction_gene_mapping_file = joinpath(path_to_model_directory,"ReactionGeneMap.dat")
    export_reaction_tag_to_gene_mapping_file(cobra_dictionary, kegg_organism_code, path_reaction_gene_mapping_file)

    # reaction tag to ec mapping -
    path_reaction_to_ec_number_mapping = joinpath(path_to_model_directory, "ReactionECNumberMap.dat")
    export_reaction_tag_to_ec_mapping_file(cobra_dictionary, kegg_organism_code, path_reaction_to_ec_number_mapping)
end

function write_cobra_model_file(path_to_cobra_mat_file::String,cobra_dictionary::Dict{String,Any})

    # check the dir -
    is_dir_path_ok(path_to_cobra_mat_file)

    # get the model name -
    model_name = splitdir(splitext(path_to_cobra_mat_file)[1])[2]

    # ok, the dir is legit, so let's write out the cobra file -
    file = matopen(path_to_cobra_mat_file, "w")
    write(file, "$(model_name)", cobra_dictionary)
    close(file)
end

function generate_cobra_dictionary_from_vff(path_to_vff_file::String, path_to_mapping_file_dir::String)::Dict{String,Any}

    # check the path -
    is_file_path_ok(path_to_vff_file)

    # does the mapping file dir exist?
    is_dir_path_ok(path_to_mapping_file_dir)

    # initialize -
    cobra_dictionary = Dict{String,Any}()

    # parse the sections of the vff -
    problem_dictionary = load_vff_model_file(path_to_vff_file)

    # get some stuff from the problem dictionary -
    list_of_reactions = problem_dictionary["list_of_reactions"]
    reaction_order_array = problem_dictionary["reaction_order_array"]
    gene_order_array = problem_dictionary["gene_order_array"]
    list_of_genes = problem_dictionary["list_of_genes"]

    # let's fill in the fields for the cobra dictionary -
    # rxns: reaction tags -
    cobra_dictionary["rxns"] = reaction_order_array

    # rxnNames: build reaction name array -
    reaction_name_array = Array{String,1}()
    path_reaction_name_mapping_file = joinpath(path_to_mapping_file_dir,"ReactionNameMap.dat")
    is_file_path_ok(path_reaction_name_mapping_file)
    reaction_tag_name_map = build_mapping_dictionary(path_reaction_name_mapping_file)
    for reaction_tag in reaction_order_array

        # check, do we have this tag in the map?
        if (haskey(reaction_tag_name_map,reaction_tag) == true)
            reaction_name_mapping_wrapper = reaction_tag_name_map[reaction_tag]
            value_set = reaction_name_mapping_wrapper.value
            push!(reaction_name_array, pop!(value_set))
        else
            push!(reaction_name_array,reaction_tag)
        end
    end
    cobra_dictionary["rxnNames"] = reaction_name_array

    # genes: list of genes -
    cobra_dictionary["genes"] = gene_order_array

    # rules: rules -
    rule_text_array = Array{String,1}()
    list_of_rules = problem_dictionary["list_of_rules"]
    for reaction_tag in reaction_order_array

        # do we have this reaction_tag in our rule list?
        if (haskey(list_of_rules,reaction_tag) == true)
            rule_text = list_of_rules[reaction_tag].rule_text
            if (rule_text != "[]")
                push!(rule_text_array,rule_text)
            else
                push!(rule_text_array,"")
            end
        else
            push!(rule_text_array,"")
        end
    end
    cobra_dictionary["rules"] = rule_text_array

    # rxnECNumbers: Add the ec numbers of the cobra dictionary -
    ec_number_array = Array{String,1}()
    for reaction_tag in reaction_order_array

        # check - do we have this reaction tag?
        if (haskey(list_of_reactions,reaction_tag) == true)

            # get the wrapper -
            reaction_wrapper = list_of_reactions[reaction_tag]

            # get the ec_number field -
            ec_number_record = reaction_wrapper.ec_number

            # cache -
            push!(ec_number_array, ec_number_record)
        else
            # for some reason, we don't have this reaction tag?
            # if this happens, we have some crazy sh*tnitz happening, better throw an exception
            throw(ErrorException("reaction tag not in the list of reactions."))
        end
    end
    cobra_dictionary["rxnECNumbers"] = ec_number_array

    # mets: add metabolite symbols -
    # get list of reaction wrappers in the proper order -
    reaction_wrapper_array = Array{CBMetabolicReaction,1}()
    for reaction_tag in reaction_order_array

        # do we have this reaction?
        reaction_wrapper = list_of_reactions[reaction_tag]
        push!(reaction_wrapper_array, reaction_wrapper)
    end

    # Get the metabolites from the reaction list -
    metabolite_wrapper_array = build_metabolite_symbol_array(reaction_wrapper_array)

    # lb: lower flux bound (from the vff) -
    lb_array = Array{Float64,1}()
    rev_array = Float64[]
    for reaction_wrapper in reaction_wrapper_array

        # get lower bound -
        reverse_bound_value = reaction_wrapper.reverse
        if (reverse_bound_value == "-inf")
            push!(lb_array,-1000.0)
        else
            push!(lb_array,parse(Float64,reverse_bound_value))
        end

        if (parse(Float64,reverse_bound_value)<0.0)
            push!(rev_array,1.0)
        else
            push!(rev_array,0.0)
        end
    end
    cobra_dictionary["lb"] = lb_array
    cobra_dictionary["rev"] = rev_array

    # ub: upper flux bound (from the vff) -
    ub_array = Array{Float64,1}()
    for reaction_wrapper in reaction_wrapper_array

        # get upper bound -
        forward_bound_value = reaction_wrapper.forward
        if (forward_bound_value == "inf")
            push!(ub_array,1000.0)
        else
            push!(ub_array,parse(Float64,forward_bound_value))
        end
    end
    cobra_dictionary["ub"] = ub_array

    # build metabolite name array -
    metabolite_symbol_order_array = problem_dictionary["metabolite_symbol_order_array"]
    list_of_metabolite_wrappers = problem_dictionary["list_of_metabolites"]
    met_names_array = Array{String,1}()
    met_kegg_id_array = Array{String,1}()
    met_symbol_array = Array{String,1}()
    for metabolite_symbol in metabolite_symbol_order_array

        # do we have this metabolite wrapper?
        if (haskey(list_of_metabolite_wrappers,metabolite_symbol) == true)

            # get the met name -
            met_wrapper = list_of_metabolite_wrappers[metabolite_symbol]
            met_symbol = met_wrapper.symbol
            met_name = met_wrapper.name
            met_kegg_id = met_wrapper.kegg_id

            # cache -
            push!(met_symbol_array, met_symbol)
            push!(met_names_array,met_name)
            push!(met_kegg_id_array,met_kegg_id)
        else

            # cache -
            push!(met_symbol_array,"")
            push!(met_names_array,"")
            push!(met_kegg_id_array,"")
        end
    end
    cobra_dictionary["metNames"] = met_names_array
    cobra_dictionary["metKEGGID"] = met_kegg_id_array
    cobra_dictionary["mets"] = met_symbol_array

    # S: generate the stoichiometric matrix -
    S = build_stoichiometric_matrix(met_symbol_array,reaction_wrapper_array)
    cobra_dictionary["S"] = S

    # c: objective vector - default is 0's
    # how many reactions do we have?
    (number_of_species,number_of_reactions) = size(S)
    c_vector = zeros(number_of_reactions)
    cobra_dictionary["c"] = c_vector

    # rxnGeneMat - reaction to gene mapping array -
    number_of_genes = length(gene_order_array)
    rxnGeneMat = zeros(number_of_reactions, number_of_genes)
    for (row_index, reaction_tag) in enumerate(reaction_order_array)
        for (col_index, gene_symbol) in enumerate(gene_order_array)

            # ok, get the reation-gene map for this reaction -
            if (haskey(list_of_genes, reaction_tag) == true)

                # grab the genes -
                gene_set = list_of_genes[reaction_tag].genes

                # does this set contain the gene?
                if (in(gene_symbol,gene_set) == true)
                    rxnGeneMat[row_index,col_index] = 1.0
                end
            end
        end
    end
    cobra_dictionary["rxnGeneMat"] = rxnGeneMat

    # return -
    return cobra_dictionary
end

function export_reaction_mapping_file(cobra_dictionary::Dict{String,Any}, path_to_mapping_file::String)

    # Check the dictionary -
    is_cobra_dictionary_ok(cobra_dictionary)

    # initalize -
    mapping_buffer = String[]

    # add header text to mapping file -
    filename = splitdir(path_to_mapping_file)[2]
    +(mapping_buffer,"// ------------------------------------------------------------- //")
    +(mapping_buffer,"// $(filename)")
    +(mapping_buffer,"// GENERATED BY: CBModelTools")
    +(mapping_buffer,"// GENERATED ON: $(Dates.now())")
    +(mapping_buffer,"// SOURCE: https://github.com/varnerlab/CBModelTools")
    +(mapping_buffer,"//")
    +(mapping_buffer,"// Reaction->reaction name map")
    +(mapping_buffer,"// record: reacton_tag=reaction_name")
    +(mapping_buffer,"// reaction_tag: from the rxns field of the COBRA mat file")
    +(mapping_buffer,"// reaction_name: from the rxnNames field of the COBRA mat file")
    +(mapping_buffer,"//")
    +(mapping_buffer,"// The order of the reactions in this file is used to order the")
    +(mapping_buffer,"// columns of the stoichiometric matrix.")
    +(mapping_buffer,"// ------------------------------------------------------------ //")

    list_of_reactions = cobra_dictionary["rxns"]
    list_of_reaction_names = cobra_dictionary["rxnNames"]

    # process each gene -
    for (index,reaction_tag) in enumerate(list_of_reactions)

        # init -
        line = ""

        # get content -
        reaction_name = list_of_reaction_names[index]

        # build line -
        line *= "$(reaction_tag)=$(reaction_name)"

        # cache -
        +(mapping_buffer, line)
    end

    # Write out the mapping file
    write_file_to_path(path_to_mapping_file,mapping_buffer)
end

function export_gene_order_file(cobra_dictionary::Dict{String,Any}, kegg_organism_code::String, path_to_mapping_file::String)

    # Check organism code -
    is_kegg_organism_code_ok(kegg_organism_code)

    # Check the dictionary -
    is_cobra_dictionary_ok(cobra_dictionary)

    # initalize -
    mapping_buffer = String[]
    list_of_genes = cobra_dictionary["genes"]

    # process each gene -
    for gene_id in list_of_genes

        # make the KEGG gene id -
        kegg_gene_id = "$(kegg_organism_code):$(gene_id)"

        # cache -
        +(mapping_buffer, kegg_gene_id)
    end

    # Write out the vff file -
    open("$(path_to_mapping_file)", "w") do f

        for line_item in mapping_buffer
            write(f,"$(line_item)\n")
        end
    end
end

function export_reaction_tag_to_gene_mapping_file(cobra_dictionary::Dict{String,Any}, kegg_organism_code::String, path_to_mapping_file::String)

    # Check organism code -
    is_kegg_organism_code_ok(kegg_organism_code)

    # Check the dictionary -
    is_cobra_dictionary_ok(cobra_dictionary)

    # initalize -
    mapping_buffer = String[]

    # add header text to mapping file -
    filename = splitdir(path_to_mapping_file)[2]
    +(mapping_buffer,"// --------------------------------------------------------- //")
    +(mapping_buffer,"// $(filename)")
    +(mapping_buffer,"// GENERATED BY: CBModelTools")
    +(mapping_buffer,"// GENERATED ON: $(Dates.now())")
    +(mapping_buffer,"// SOURCE: https://github.com/varnerlab/CBModelTools")
    +(mapping_buffer,"")
    +(mapping_buffer,"// Reaction->gene mapping - ")
    +(mapping_buffer,"// record: reacton_tag={gene_id}")
    +(mapping_buffer,"// reaction_tag: from the rxns field of the COBRA mat file")
    +(mapping_buffer,"// gene_id: (kegg organism code):(gene location code)")
    +(mapping_buffer,"// --------------------------------------------------------- //")

    # get the list of reaction tags, and genes from the cobra dictionary -
    list_of_reaction_tags = cobra_dictionary["rxns"]
    list_of_genes = cobra_dictionary["genes"]

    # get the rxn gene mapping matrix -
    RGM = Matrix(cobra_dictionary["rxnGeneMat"])

    # What is the size of the system?
    (number_of_reactions, number_of_genes) = size(RGM)

    # main loop -
    for reaction_index = 1:number_of_reactions

        # what is the tag for this reaction?
        reaction_tag = list_of_reaction_tags[reaction_index]

        # initialize -
        record = ""
        record *= "$(reaction_tag)="

        # does this reaction have associated genes?
        idx_genes = findall(x->x==1,RGM[reaction_index,:])
        local_gene_set = Set{String}()
        if (!isempty(idx_genes))

            # ok, we have genes, grab them -
            gene_id_array = list_of_genes[idx_genes]

            # process each gene -
            for gene_id in gene_id_array

                # make the KEGG gene id -
                kegg_gene_id = gene_id
                if (occursin(".",gene_id) == true)
                    # need to cutoff the trailing *.1
                    kegg_gene_id = "$(kegg_organism_code):$(gene_id[1:end-2])"
                else
                    kegg_gene_id = "$(kegg_organism_code):$(gene_id)"
                end

                # cache -
                push!(local_gene_set,kegg_gene_id)
            end
        end

        if (isempty(local_gene_set) == false)

            # make an ec record -
            gene_record = ""
            for gene_number in local_gene_set
                gene_record*="$(gene_number),"
            end

            # remove trailing ,
            gene_record = gene_record[1:end-1]

            # add the record -
            record*="$(gene_record)"

            # push into buffer and go around again -
            +(mapping_buffer,record)
        end
    end

    # Write out the file -
    write_file_to_path(path_to_mapping_file, mapping_buffer)
end

function export_reaction_tag_to_ec_mapping_file(cobra_dictionary::Dict{String,Any}, kegg_organism_code::String, path_to_mapping_file::String)

    # Check organism code -
    is_kegg_organism_code_ok(kegg_organism_code)

    # Check the dictionary -
    is_cobra_dictionary_ok(cobra_dictionary)

    # initalize -
    mapping_buffer = String[]

    # add header text to mapping file -
    filename = splitdir(path_to_mapping_file)[2]
    +(mapping_buffer,"// --------------------------------------------------------- //")
    +(mapping_buffer,"// $(filename)")
    +(mapping_buffer,"// GENERATED BY: CBModelTools")
    +(mapping_buffer,"// GENERATED ON: $(Dates.now())")
    +(mapping_buffer,"// SOURCE: https://github.com/varnerlab/CBModelTools")
    +(mapping_buffer,"//")
    +(mapping_buffer,"// Reaction->ecnumber mapping - ")
    +(mapping_buffer,"// record: reacton_tag={ecnumber}")
    +(mapping_buffer,"// reaction_tag: from the rxns field of the COBRA mat file")
    +(mapping_buffer,"// ec_number: possible ecnumber estimated from KEGG")
    +(mapping_buffer,"// --------------------------------------------------------- //")


    # get the list of reaction tags, and genes from the cobra dictionary -
    list_of_reaction_tags = cobra_dictionary["rxns"]
    list_of_genes = cobra_dictionary["genes"]

    # get the rxn gene mapping matrix -
    RGM = Matrix(cobra_dictionary["rxnGeneMat"])

    # What is the size of the system?
    (number_of_reactions, number_of_genes) = size(RGM)

    # Declare a progress meter for user feedback -
    p = Progress(number_of_reactions,color=:yellow)

    # main loop -
    for reaction_index = 1:number_of_reactions

        # what is the tag for this reaction?
        reaction_tag = list_of_reaction_tags[reaction_index]

        # user message -
        msg = "Starting $(reaction_tag) ($(reaction_index) of $(number_of_reactions))"

        # update the progress bar -
        ProgressMeter.next!(p; showvalues = [(:status,msg)])

        # initialize -
        record = ""
        record *= "$(reaction_tag)="

        # does this reaction have associated genes?
        idx_genes = findall(x->x==1,RGM[reaction_index,:])

        # init -
        ec_record_set = Set{String}()
        if (!isempty(idx_genes))

            # ok, we have genes, grab them -
            gene_id_array = list_of_genes[idx_genes]

            # process each gene -
            for gene_id in gene_id_array

                # make the KEGG gene id -
                kegg_gene_id = gene_id
                if (occursin(".",gene_id) == true)
                    # need to cutoff the trailing *.1
                    kegg_gene_id = "$(kegg_organism_code):$(gene_id[1:end-2])"
                else
                    kegg_gene_id = "$(kegg_organism_code):$(gene_id)"
                end

                # lookup -
                local_ec_record_set = lookup_gene_ec_mapping_record(kegg_gene_id)
                if local_ec_record_set != nothing

                    # push into ec_record_set -
                    for ec_number in local_ec_record_set
                        push!(ec_record_set,ec_number)
                    end
                end
            end
        end

        if (isempty(ec_record_set) == false)

            # make an ec record -
            ec_record = ""
            for ec_number in ec_record_set
                ec_record*="$(ec_number),"
            end

            # remove trailing ,
            ec_record = ec_record[1:end-1]

            # add the record -
            record*="$(ec_record)"

            # push into buffer and go around again -
            +(mapping_buffer,record)
        end
    end

    # Write out the file -
    write_file_to_path(path_to_mapping_file,mapping_buffer)
end
# ------------------------------------------------------------------------------ #
