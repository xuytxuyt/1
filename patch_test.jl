
using ApproxOperator, YAML

config = YAML.load_file("./yml/patch_test.yml")
config_λ = YAML.load_file("./yml/patch_test_lm.yml")
elements, nodes = importmsh("./msh/test.msh",config)
elements_λ, nodes_λ = importmsh("./msh/test_lm.msh",config_λ)
nₚ = length(nodes)
nₗ = length(nodes_λ)

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements_λ["Γᵍ"])

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
g = zeros(nₚ,nₗ)
q = zeros(nₗ)

ops[1](elements["Ω"],k)
# ops[4](elements["Γᵍ"],k,f)
ops[5](elements["Γᵍ"],elements_λ["Γᵍ"],g,q)

d = [k g;g' zeros(nₗ,nₗ)]\[f;q]

push!(nodes,:d=>d)

set𝓖!(elements["Ω"],:TriGI13)
set_memory_𝝭!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
prescribe!(elements["Ω"],:u=>(x,y,z)->1.0+2x+3y)
prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->2.0)
prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->3.0)
h1,l2 = ops[6](elements["Ω"])