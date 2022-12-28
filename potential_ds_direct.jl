using Revise, TOML, ApproxOperator

config = TOML.parsefile("./toml/ds_direct.toml")
elements,nodes = importmsh("./msh/square_1.msh",config)

nₚ = length(elements["∂Ω"])

ApproxOperator.set𝓖_DB!(elements["Γ"],:SegGI5)
set_memory_𝝭!(elements["Γ"],:𝝭,:𝝭̄,:∂𝝭∂x,:∂𝝭∂y)
set𝝭!(elements["Γ"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])

# prescribing
r = 5
u(x,y,z) = (x+y)^r
∂u∂x(x,y,z) = r*(x+y)^abs(r-1)
∂u∂y(x,y,z) = r*(x+y)^abs(r-1)
∂²u∂x²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
∂²u∂y²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
t(x,y,z,n₁,n₂) = ∂u∂x(x,y,z)*n₁+∂u∂y(x,y,z)*n₂
b(x,y,z) = -(∂²u∂x²(x,y,z)+∂²u∂y²(x,y,z))

prescribe!(elements["Ω"],:b=>b)
prescribe!(elements["Γᵍ"],:u=>u)
prescribe!(elements["Γ"],:t=>t)
prescribe!(elements["∂Ω"],:u=>u)

ops = [
    Operator{:∫∫∇v∇udxdy}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫udΓ}(),
    Operator{:H₁}(),
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)
# ops[3](elements["Γ"],f)
ops[4](elements["Γᵍ"],k,f)

d = k\f

dex = ops[4](elements["∂Ω"])

push!(getfield(nodes[1],:data),:d=>(2,d))

prescribe!(elements["Ω"],:u=>u)
prescribe!(elements["Ω"],:∂u∂x=>∂u∂x)
prescribe!(elements["Ω"],:∂u∂y=>∂u∂y)
H₁,L₂ = ops[5](elements["Ω"])
