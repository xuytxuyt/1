
using ApproxOperator, YAML

config = YAML.load_file("./yml/patch_test.yml")
elements, nodes = importmsh("./msh/test.msh",config)
nₚ = length(nodes)

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])

prescribe!(elements["Γᵍ"],:g=>(x,y,z)->1.0+2x+3y)

ops = [
    Operator{:∫∇v∇udΩ}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫vgdΓ}(:α=>1e7),
    Operator{:∫λgdΓ}(),
    Operator{:H₁}()
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω"],k)
ops[4](elements["Γᵍ"],k,f)

d = k\f

push!(nodes,:d=>d)

set𝓖!(elements["Ω"],:TriGI13)
set_memory_𝝭!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
prescribe!(elements["Ω"],:u=>(x,y,z)->1.0+2x+3y)
prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->2.0)
prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->3.0)
h1,l2 = ops[6](elements["Ω"])