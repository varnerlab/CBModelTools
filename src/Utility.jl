function is_file_path_ok(path_to_file::String)

    if (isfile(path_to_file) == false)
        throw(ArgumentError("$(path_to_file) does not exist."))
    end
end

function is_dir_path_ok(path_to_file::String)

    # get the dir of the path -
    dir_name = dirname(path_to_file)

    # check, is this a legit dir?
    if (isdir(dir_name) == false)
        throw(ArgumentError("$(dir_name) does not exist."))
    end
end

function write_file_to_path(path_to_file::String,buffer::Array{String,1})

    # check, is this a legit dir?
    is_dir_path_ok(path_to_file)

    # write -
    open("$(path_to_file)", "w") do f

        for line_item in buffer
            write(f,"$(line_item)\n")
        end
    end
end

function read_file_from_path(path_to_file::String)::Array{String,1}

    # is this mapping file path legit?
    is_file_path_ok(path_to_file)

    # initialize -
    buffer = String[]

    # Read in the file -
    open("$(path_to_file)", "r") do file
        for line in eachline(file)
            +(buffer,line)
        end
    end

    # return -
    return buffer
end

function extract_section(file_buffer_array::Array{String,1},start_section::String,end_section::String)

    # initialize -
    section_buffer = String[]

    # find the SECTION START AND END -
    section_line_start = 1
    section_line_end = 1
    for (index,line) in enumerate(file_buffer_array)
        if (occursin(start_section,line) == true)
            section_line_start = index
        elseif (occursin(end_section,line) == true)
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

function build_mapping_dictionary(path_to_mapping_file::String;record_delim::String="=",field_delim::String=",")::Dict{String,Mapping}

    # check the file -
    is_file_path_ok(path_to_mapping_file)

    # initalize -
    mapping_dictionary = Dict{String,Mapping}()

    # ok, read -
    file_record_array = read_file_from_path(path_to_mapping_file)

    # process the record -
    for record in file_record_array

        # ignore comments -
        if (occursin("//",record) == false)

            # create map entry -
            mapping = Mapping()

            # create set -
            token_set = Set{String}()

            # record -
            token_array = split(record,record_delim)

            # the key is the first entry -
            key = token_array[1]
            mapping.key = key

            # push the tokens (genes, ecnumbers etc) into set -
            token_array = split(token_array[2],field_delim)
            for token in token_array

                # grab -
                push!(token_set,token)
            end

            # set the value -
            mapping.value = token_set

            # add entry to dictionary -
            mapping_dictionary[key] = mapping
        end
    end

    # return -
    return mapping_dictionary
end
