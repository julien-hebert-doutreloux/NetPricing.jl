struct BBLogEntry
    parent::BBNode
    branch_arc::Union{Int, Nothing}
    left::Union{BBNode, Nothing}
    right::Union{BBNode, Nothing}
    bestsolution::FeasiblePoint
end