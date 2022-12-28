using Revise, TOML, ApproxOperator

config = TOML.parsefile("./toml/ds_lm.toml")
elements,nodes = importmsh("./msh/square_1.msh",config)

ApproxOperator.set𝓖_DB!(elements["Γ"],:SegGI5)
set_memory_𝝭!(elements["Γ"],:𝝭,:𝝭̄,:∂𝝭∂x,:∂𝝭∂y)
set𝝭!(elements["Γ"])
set∇𝝭!(elements["Γ"])
set𝝭̄!(elements["Γ"])

# elements["Γ"][1].𝓖[1].𝑤=0.
# elements["Γ"][1].𝓖[4].𝑤=0.
# elements["Γ"][1].𝓖[7].𝑤=0.
# elements["Γ"][1].𝓖[10].𝑤=0.
# elements["Γ"][1].𝓖[13].𝑤=0.
# elements["Γ"][2].𝓖[3].𝑤=0.
# elements["Γ"][2].𝓖[6].𝑤=0.
# elements["Γ"][2].𝓖[9].𝑤=0.
# elements["Γ"][2].𝓖[12].𝑤=0.
# elements["Γ"][2].𝓖[15].𝑤=0.

# elements["Γ"][1].𝓖[2].𝑤=0.
# elements["Γ"][1].𝓖[3].𝑤=0.
# elements["Γ"][1].𝓖[5].𝑤=0.
# elements["Γ"][1].𝓖[6].𝑤=0.
# elements["Γ"][1].𝓖[8].𝑤=0.
# elements["Γ"][1].𝓖[9].𝑤=0.
# elements["Γ"][1].𝓖[11].𝑤=0.
# elements["Γ"][1].𝓖[12].𝑤=0.
# elements["Γ"][1].𝓖[14].𝑤=0.
# elements["Γ"][1].𝓖[15].𝑤=0.
# elements["Γ"][2].𝓖[1].𝑤=0.
# elements["Γ"][2].𝓖[2].𝑤=0.
# elements["Γ"][2].𝓖[4].𝑤=0.
# elements["Γ"][2].𝓖[5].𝑤=0.
# elements["Γ"][2].𝓖[7].𝑤=0.
# elements["Γ"][2].𝓖[8].𝑤=0.
# elements["Γ"][2].𝓖[10].𝑤=0.
# elements["Γ"][2].𝓖[11].𝑤=0.
# elements["Γ"][2].𝓖[13].𝑤=0.
# elements["Γ"][2].𝓖[14].𝑤=0.

nₚ = getnₚ(elements["Ω"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])

# prescribing
r = 3
u(x,y,z) = (x+y)^r
∂u∂x(x,y,z) = r*(x+y)^abs(r-1)
∂u∂y(x,y,z) = r*(x+y)^abs(r-1)
∂²u∂x²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
∂²u∂y²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
t(x,y,z,n₁,n₂) = ∂u∂x(x,y,z)*n₁+∂u∂y(x,y,z)*n₂
b(x,y,z) = -(∂²u∂x²(x,y,z)+∂²u∂y²(x,y,z))
ū(x,y,z,n₁,n₂) = sign(n₁+n₂)*(x+y)^r

prescribe!(elements["Ω"],:b=>b)
prescribe!(elements["Γᵍ"],:u=>ū)
prescribe!(elements["Γᵍ"],:g=>u)
prescribe!(elements["∂Ω"],:u=>u)

ops = [
    Operator{:∫∫∇v∇udxdy}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫∇𝑛uvdΓ}(),
    Operator{:𝑓𝑣}(),
    Operator{:H₁}(),
    Operator{:∫udΓ}(),
    Operator{:∫vgdΓ}(:α=>1e9)
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
g = zeros(nₚ,4)
q = zeros(4)

# ops[1](elements["Ω"],k)
# ops[2](elements["Ω"],f)
ops[3](elements["Γ"],g)
# ops[4](elements["Γᵍ"],q)
# ops[4](elements["∂Ω"][2],k,f)
# ops[7](elements["Γᵍ"],k,f)
# ops[8](elements["Γ"],g)

# d = [k g;g' zeros(4,4)]\[f;.-q]
# d = k\f
# d_ = [k g[:,[2,4]];g[:,[2,4]]' zeros(2,2)]\[f;zeros(2)]
# push!(getfield(nodes[1],:data),:d=>(2,d))

dex = ops[6](elements["∂Ω"])

# d[1:5] - dex
# prescribe!(elements["Ω"],:u=>u)
# prescribe!(elements["Ω"],:∂u∂x=>∂u∂x)
# prescribe!(elements["Ω"],:∂u∂y=>∂u∂y)
# H₁,L₂ = ops[5](elements["Ω"])
