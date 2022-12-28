using Revise, TOML, ApproxOperator

config = TOML.parsefile("./toml/fem_lm.toml")
elements,nodes = importmsh("./msh/square_1.msh",config)

nₑ = length(elements["Ω"])
nₚ = length(nodes)

ApproxOperator.set𝓖_DB!(elements["Γ"],:SegGI5)
set_memory_𝝭!(elements["Γ"],:𝝭)
set𝝭!(elements["Γ"])

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

elements["Γ"][1].𝓖[2].𝑤=0.
elements["Γ"][1].𝓖[3].𝑤=0.
elements["Γ"][1].𝓖[5].𝑤=0.
elements["Γ"][1].𝓖[6].𝑤=0.
elements["Γ"][1].𝓖[8].𝑤=0.
elements["Γ"][1].𝓖[9].𝑤=0.
elements["Γ"][1].𝓖[11].𝑤=0.
elements["Γ"][1].𝓖[12].𝑤=0.
elements["Γ"][1].𝓖[14].𝑤=0.
elements["Γ"][1].𝓖[15].𝑤=0.
elements["Γ"][2].𝓖[1].𝑤=0.
elements["Γ"][2].𝓖[2].𝑤=0.
elements["Γ"][2].𝓖[4].𝑤=0.
elements["Γ"][2].𝓖[5].𝑤=0.
elements["Γ"][2].𝓖[7].𝑤=0.
elements["Γ"][2].𝓖[8].𝑤=0.
elements["Γ"][2].𝓖[10].𝑤=0.
elements["Γ"][2].𝓖[11].𝑤=0.
elements["Γ"][2].𝓖[13].𝑤=0.
elements["Γ"][2].𝓖[14].𝑤=0.


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
# prescribe!(elements["Γᵗ"],:t=>t)
prescribe!(elements["∂Ω"],:u=>u)
prescribe!(elements["Ω"],:∂u∂x=>∂u∂x)
prescribe!(elements["Ω"],:∂u∂y=>∂u∂y)

k = zeros(nₚ,nₚ)
op_Ω = Operator{:∫∫∇v∇udxdy}(:k=>1.0)
op_Ω(elements["Ω"],k)
f = zeros(nₚ)
op_b = Operator{:∫vbdΩ}()
op_b(elements["Ω"],f)

g = zeros(nₚ,2*nₑ)
op_Γ = Operator{:∫usᵢnᵢdΓ}()
op_Γ(elements["Γ"],g)
q = zeros(2*nₑ)
op_Γᵍ = Operator{:∫gsᵢnᵢdΓ}()
op_Γᵍ(elements["Γᵍ"],q)

d = [k g;g' zeros(2*nₑ,2*nₑ)]\[f;q]

op_∇u = Operator{:∫∇udΩ}()
dx,dy = op_∇u(elements["Ω"])
# ops = [
#     Operator{:∫∫∇v∇udxdy}(:k=>1.0),
#     Operator{:∫vbdΩ}(),
#     Operator{:∫vtdΓ}(),
#     Operator{:∫∇𝑛uvdΓ}(),
#     Operator{:∫∇𝑛ugdΓ}(),
#     Operator{:H₁}(),
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
