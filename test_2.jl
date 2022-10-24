using Revise, ApproxOperator

elements, nodes = importmsh("./msh/test.msh",Val(:test))
# elements_lam, nodes_λ = importmsh("./msh/test.msh",config)

set_memory_𝝭!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
set_memory_𝝭!(elements["Γ"],:𝝭)
set_memory_𝝭!(elements["Γ_λ"],:𝝭)
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ"])
set𝝭!(elements["Γ_λ"])

prescribe!(elements["Γ"],:g=>(x,y,z)->x*y)
prescribe!(elements["∂Ω"],:u=>(x,y,z)->x*y)

ops = [
    Operator{:∫∇v∇udΩ}(:k=>1.0),
    Operator{:∫λₙgdΓ}(),
    Operator{:∫udΓ}()
]

k = zeros(5,5)
f = zeros(5)
g = zeros(5,4)
q = zeros(4)

ops[1](elements["Ω"],k)
ops[2](elements["Γ"],elements["Γ_λ"],g,q)

d = [k g;g' zeros(4,4)]\[f;q]

dext = ops[3](elements["∂Ω"])