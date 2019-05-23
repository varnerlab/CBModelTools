mutable struct Mapping

    key::String
    value::Set{String}

    function Mapping()
		this = new()
	end
end

mutable struct CBMetabolite

    symbol::String
    name::String
    kegg_id::String

    function CBMetabolite()
		this = new()
	end
end

mutable struct CBRule

    reaction_name::String
    rule::String

    function CBRule()
		this = new()
	end
end

mutable struct CBMetabolicReaction

    record::String
    reaction_name::String
    ec_number::String
    left_phrase::String
    right_phrase::String
    reverse::String
    forward::String

    function CBMetabolicReaction()
		this = new()
	end
end
