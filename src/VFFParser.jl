function extract_section(file_buffer_array::Array{String,1},start_section::String,end_section::String)

    # initialize -
    section_buffer = String[]

    # find the SECTION START AND END -
    section_line_start = 1
    section_line_end = 1
    for (index,line) in enumerate(file_buffer_array)
        if (line == start_section)
            section_line_start = index
        elseif (line == end_section)
            section_line_end = index
        end
    end

    for line_index = (section_line_start+1):(section_line_end-1)
        line_item = file_buffer_array[line_index]
        push!(section_buffer,line_item)
    end

    # return -
    return section_buffer
end

function load_vff_model_file(path_to_vff_file::String)::Dict{String,Any}

    # check the file path -
    is_file_path_ok(path_to_vff_file)

    # initialize -
    problem_dictionary = Dict{String,Any}()

    # load the file into an array -
    file_buffer_array = read_file_from_path(path_to_vff_file)

    # parse the reaction section -
    reaction_section = extract_section(file_buffer_array,"#REACTION::START","#REACTION::END")

    @show reaction_section

    # return -
    return problem_dictionary
end
