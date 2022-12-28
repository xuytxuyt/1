using Revise, TOML, ApproxOperator

config = TOML.parsefile("./toml/ds_penalty.toml")
# elements,nodes = importmsh("./msh/square_1_naturalbc.msh",config)
elements,nodes = importmsh("./msh/square_1.msh",config)

ApproxOperator.set𝓖_DB!(elements["Γ"],:SegGI5)
set_memory_𝝭!(elements["Γ"],:𝝭,:𝝭̄)
set𝝭!(elements["Γ"])
set𝝭̄!(elements["Γ"])

# nₚ = getnₚ(elements["Ω"])

# set𝝭!(elements["Ω"])
# set∇𝝭!(elements["Ω"])
# set𝝭!(elements["Γᵍ"])
# set𝝭!(elements["Γᵗ"])

# prescribing
r = 4
u(x,y,z) = (x+y)^r
∂u∂x(x,y,z) = r*(x+y)^abs(r-1)
∂u∂y(x,y,z) = r*(x+y)^abs(r-1)
∂²u∂x²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
∂²u∂y²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
t(x,y,z,n₁,n₂) = ∂u∂x(x,y,z)*n₁+∂u∂y(x,y,z)*n₂
b(x,y,z) = -(∂²u∂x²(x,y,z)+∂²u∂y²(x,y,z))

# prescribe!(elements["Ω"],:b=>b)
# prescribe!(elements["Γᵗ"],:t=>t)
# prescribe!(elements["Γᵍ"],:g=>u)
prescribe!(elements["Γ"],:u=>u)
prescribe!(elements["∂Ω"],:u=>u)

# ops = [
#     Operator{:∫∫∇v∇udxdy}(:k=>1.0),
#     Operator{:∫vbdΩ}(),
#     Operator{:∫vtdΓ}(),
#     Operator{:∫vgdΓ}(:α=>1e0),
#     Operator{:H₁}(),
#     Operator{:∫udΓ}()
# ]

# k = zeros(nₚ,nₚ)
# f = zeros(nₚ)

# ops[1](elements["Ω"],k)
# ops[2](elements["Ω"],f)
# ops[3](elements["Γᵗ"],f)
# ops[4](elements["Γᵍ"][1],k,f)
# ops[4](elements["∂Ω"][2],k,f)

# d = k\f
# push!(getfield(nodes[1],:data),:d=>(2,d))

op_ex = Operator{:∫udΓ}()
dex = op_ex(elements["∂Ω"])

# d - dex
# prescribe!(elements["Ω"],:u=>u)
# prescribe!(elements["Ω"],:∂u∂x=>∂u∂x)
# prescribe!(elements["Ω"],:∂u∂y=>∂u∂y)
# H₁,L₂ = ops[5](elements["Ω"])

g = zeros(5,4)
q = zeros(4)
op_Γ = Operator{:∫vu𝑛dΓ}()
op_Γ(elements["Γ"],g,q)