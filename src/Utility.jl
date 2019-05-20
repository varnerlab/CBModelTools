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
