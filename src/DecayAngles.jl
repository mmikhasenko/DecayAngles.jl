module DecayAngles

using AbstractTrees

export DecayNode
export map_tree
export map_with_parent
export to_table
export nodevalue
export children
include("decaynode.jl")

export flatten_nested_tuple
export flatten_sort_nested_tuple
include("utils.jl")

end # module DecayAngles
