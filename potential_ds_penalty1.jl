using Revise, TOML, ApproxOperator

config = TOML.parsefile("./toml/ds_penalty.toml")
elements,nodes = importmsh("./msh/square_1.msh",config)

ApproxOperator.set𝓖_DB!(elements["Ω_"],:SegGI5)
set_memory_𝝭!(elements["Ω_"],:𝝭)
set𝝭!(elements["Ω_"])
elements["Ω_"][1].𝓖[1].𝑤=0.
elements["Ω_"][1].𝓖[4].𝑤=0.
elements["Ω_"][1].𝓖[7].𝑤=0.
elements["Ω_"][1].𝓖[10].𝑤=0.
elements["Ω_"][1].𝓖[13].𝑤=0.
elements["Ω_"][2].𝓖[3].𝑤=0.
elements["Ω_"][2].𝓖[6].𝑤=0.
elements["Ω_"][2].𝓖[9].𝑤=0.
elements["Ω_"][2].𝓖[12].𝑤=0.
elements["Ω_"][2].𝓖[15].𝑤=0.

# elements["Ω_"][1].𝓖[2].𝑤=0.
# elements["Ω_"][1].𝓖[3].𝑤=0.
# elements["Ω_"][1].𝓖[5].𝑤=0.
# elements["Ω_"][1].𝓖[6].𝑤=0.
# elements["Ω_"][1].𝓖[8].𝑤=0.
# elements["Ω_"][1].𝓖[9].𝑤=0.
# elements["Ω_"][1].𝓖[11].𝑤=0.
# elements["Ω_"][1].𝓖[12].𝑤=0.
# elements["Ω_"][1].𝓖[14].𝑤=0.
# elements["Ω_"][1].𝓖[15].𝑤=0.
# elements["Ω_"][2].𝓖[1].𝑤=0.
# elements["Ω_"][2].𝓖[2].𝑤=0.
# elements["Ω_"][2].𝓖[4].𝑤=0.
# elements["Ω_"][2].𝓖[5].𝑤=0.
# elements["Ω_"][2].𝓖[7].𝑤=0.
# elements["Ω_"][2].𝓖[8].𝑤=0.
# elements["Ω_"][2].𝓖[10].𝑤=0.
# elements["Ω_"][2].𝓖[11].𝑤=0.
# elements["Ω_"][2].𝓖[13].𝑤=0.
# elements["Ω_"][2].𝓖[14].𝑤=0.

nₚ = getnₚ(elements["Ω"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])

# prescribing
r = 5
u(x,y,z) = (x+y)^r
∂u∂x(x,y,z) = r*(x+y)^abs(r-1)
∂u∂y(x,y,z) = r*(x+y)^abs(r-1)
∂²u∂x²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
∂²u∂y²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
t(x,y,z,n₁,n₂) = ∂u∂x(x,y,z)*n₁+∂u∂y(x,y,z)*n₂
b(x,y,z) = -(∂²u∂x²(x,y,z)+∂²u∂y²(x,y,z))
ū(x,y,z,n₁,n₂) = sign(n₁+n₂)*(x+y)^r

prescribe!(elements["Ω"],:b=>b)
prescribe!(elements["Γᵍ"],:g=>u)
prescribe!(elements["Γᵍ"],:u=>ū)
prescribe!(elements["∂Ω"],:u=>u)
prescribe!(elements["Ω_"],:g=>u)

ops = [
    Operator{:∫∫∇v∇udxdy}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫vgdΓ}(:α=>1e0),
    Operator{:H₁}(),
    Operator{:∫udΓ}(),
    Operator{:∫vudΓ}(:α=>1e9),
    Operator{:𝑓𝑣}()
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

# ops[1](elements["Ω"],k)
# ops[2](elements["Ω"],f)
ops[4](elements["Γᵍ"],k,f)
# ops[7](elements["Ω_"],k,f)
# ops[8](elements["Γᵍ"],f)

d = k\f
# push!(getfield(nodes[1],:data),:d=>(2,d))

dex = ops[6](elements["∂Ω"])

d-dex
# f
# prescribe!(elements["Ω"],:u=>u)
# prescribe!(elements["Ω"],:∂u∂x=>∂u∂x)
# prescribe!(elements["Ω"],:∂u∂y=>∂u∂y)
# H₁,L₂ = ops[5](elements["Ω"])
