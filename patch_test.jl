
using ApproxOperator
include("input.jl")

elements, d = import_tri3("./msh/square_1.msh")

nₚ = length(d)

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])

prescribe!(elements["Γᵍ"],:g=>(x,y,z)->1.0+2x+3y)

ops = [
    Operator{:∫∫∇v∇udxdy}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫vgdΓ}(:α=>1e7),
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω"],k)
ops[4](elements["Γᵍ"],k,f)

d .= k\f
