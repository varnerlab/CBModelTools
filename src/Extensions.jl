import Base.+

function +(buffer::Array{String,1},content::String)
    push!(buffer,content)
end
