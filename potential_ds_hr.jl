using Revise, TOML, ApproxOperator

config = TOML.parsefile("./toml/ds_hr.toml")
elements,nodes = importmsh("./msh/square_1.msh",config)

set_memory_𝝭!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:𝝭̄,:∂𝝭̄∂x,:∂𝝭̄∂y)
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭̄!(elements["Ω"])
set∇̄𝝭!(elements["Ω"])

g_ = zeros(5,4)
op_Ω = Operator{:∫∇v̄∇udΩ}()
op_Ω(elements["Ω"],g_)

ApproxOperator.set𝓖_DB!(elements["Γ"],:SegGI5)
set_memory_𝝭!(elements["Γ"],:𝝭,:𝝭̄,:∂𝝭∂x,:∂𝝭∂y)
set𝝭!(elements["Γ"])
set∇𝝭!(elements["Γ"])
set𝝭̄!(elements["Γ"])
g = zeros(5,4)
op_Γ = Operator{:∫∇𝑛uvdΓ}()
op_Γ(elements["Γ"],g) 

g + g_

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

# nₚ = getnₚ(elements["Ω"])

# set𝝭!(elements["Ω̄"])
# set∇𝝭!(elements["Ω"])
# set𝝭̄!(elements["Γᵍ"])
# set∇𝝭!(elements["Γᵍ"])
# set𝝭!(elements["Γᵗ"])

# prescribing
# r = 2
# u(x,y,z) = (x+y)^r
# ∂u∂x(x,y,z) = r*(x+y)^abs(r-1)
# ∂u∂y(x,y,z) = r*(x+y)^abs(r-1)
# ∂²u∂x²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
# ∂²u∂y²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
# t(x,y,z,n₁,n₂) = ∂u∂x(x,y,z)*n₁+∂u∂y(x,y,z)*n₂
# b(x,y,z) = -(∂²u∂x²(x,y,z)+∂²u∂y²(x,y,z))

# prescribe!(elements["Ω̄"],:b=>b)
# prescribe!(elements["Γᵍ"],:g=>u)
# prescribe!(elements["Γᵗ"],:t=>t)
# prescribe!(elements["∂Ω"],:u=>u)

# ops = [
#     Operator{:∫∫∇v∇udxdy}(:k=>1.0),
#     Operator{:∫vbdΩ}(),
#     Operator{:∫vtdΓ}(),
#     Operator{:∫∇𝑛uvdΓ}(),
#     Operator{:∫∇𝑛ugdΓ}(),
#     Operator{:H₁}(),
#     Operator{:∫udΓ}(),
# ]

# k = zeros(nₚ,nₚ)
# f = zeros(nₚ)
# g = zeros(nₚ,4)
# q = zeros(4)

# ops[1](elements["Ω"],k)
# ops[2](elements["Ω̄"],q)
# ops[3](elements["Γᵗ"],q)
# ops[4](elements["Γ"],g)
# ops[4](elements["Γᵍ"],g,f)
# ops[4](elements["∂Ω"][2],k,f)
# ops[7](elements["Γᵍ"],k,f)
# ops[8](elements["Γ"],g)

# d = [k g;g' zeros(4,4)]\[f;.-q]
# d = k\f
# d_ = [k g[:,[2,4]];g[:,[2,4]]' zeros(2,2)]\[f;zeros(2)]
# push!(getfield(nodes[1],:data),:d=>(2,d))

# dex = ops[7](elements["∂Ω"])

# d[1:5] - dex
# prescribe!(elements["Ω"],:u=>u)
# prescribe!(elements["Ω"],:∂u∂x=>∂u∂x)
# prescribe!(elements["Ω"],:∂u∂y=>∂u∂y)
# H₁,L₂ = ops[5](elements["Ω"])
