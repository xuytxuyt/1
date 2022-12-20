using Revise, TOML, ApproxOperator

config = TOML.parsefile("./toml/ds_penalty.toml")
elements,nodes = importmsh("./msh/square_1.msh",config)

nₚ = getnₚ(elements["Ω"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])

# prescribing
r = 2
u(x,y,z) = (x+y)^r
∂u∂x(x,y,z) = r*(x+y)^abs(r-1)
∂u∂y(x,y,z) = r*(x+y)^abs(r-1)
∂²u∂x²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
∂²u∂y²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
t(x,y,z,n₁,n₂) = ∂u∂x(x,y,z)*n₁+∂u∂y(x,y,z)*n₂
b(x,y,z) = -(∂²u∂x²(x,y,z)+∂²u∂y²(x,y,z))

prescribe!(elements["Ω"],:b=>b)
prescribe!(elements["Γᵍ"],:g=>u)
prescribe!(elements["∂Ω"],:u=>u)

ops = [
    Operator{:∫∫∇v∇udxdy}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫vgdΓ}(:α=>1e9),
    Operator{:H₁}(),
    Operator{:∫udΓ}()
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)
ops[4](elements["Γᵍ"],k,f)

# d = k\f
# push!(getfield(nodes[1],:data),:d=>(2,d))

# dex = ops[6](elements["∂Ω"])

# prescribe!(elements["Ω"],:u=>u)
# prescribe!(elements["Ω"],:∂u∂x=>∂u∂x)
# prescribe!(elements["Ω"],:∂u∂y=>∂u∂y)
# H₁,L₂ = ops[5](elements["Ω"])
