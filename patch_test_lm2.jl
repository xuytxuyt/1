
using ApproxOperator, CairoMakie
include("input.jl")

elements, 𝓑 = import_trc("./msh/square_1.msh")
# elements, 𝓑 = import_tr("./msh/square_2_irregular.msh")
np = 5
# nₚ₁ = 5
# nₚ₂ = 4
# nₚ₁ = 16
# nₚ₂ = 9
# nₚ₁ = 56
# nₚ₂ = 25

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ω"])
# set𝝭!(elements["Ω̄"])
# set∇𝝭!(elements["Γ"])
set𝝭!(elements["Γ"])

# ū(x,y,z) = 1.0
ū(x,y,z) = 1.0+2x+3y
∂ū∂x(x,y,z) = 2.
∂ū∂y(x,y,z) = 3.
b̄(x,y,z) = 0.0
# ū(x,y,z) = 1.0+2x+3y+4x^2+5x*y+6y^2
# ∂ū∂x(x,y,z) = 2.0+8x+5y
# ∂ū∂y(x,y,z) = 3.0+5x+12y
# b̄(x,y,z) = -20.0
# prescribe!(elements["Γᵍ"],:g=>ū)
prescribe!(elements["Γ̄ᵍ"],:g=>ū)
prescribe!(elements["∂Ω"],:u=>ū)
prescribe!(elements["Ω"],:b=>b̄)
# prescribe!(elements["Γ̄"],:t=>(x,y,z,n₁,n₂)->∂ū∂x(x,y,z)*n₁ + ∂ū∂y(x,y,z)*n₂)
# prescribe!(elements["Γ"],:g=>ū)

ops = [
    Operator{:∫∫∇v∇udxdy}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫vᵢnᵢuds}(:k=>-1.0),
    Operator{:∫vᵢnᵢuds}(:k=>1.0),
    Operator{:∫vᵢnᵢgds}(:k=>-1.0),
    Operator{:∫vgdΓ}(:α=>1e7),
    Operator{:∫uds}(),
]

k = zeros(np,np)
k_Γ = zeros(np,np)
f₂ = zeros(np)
f = zeros(np)

ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)
ops[4](elements["Γ̄"],elements["Γ"],k_Γ)
# ops[5](elements["Γᵍ"],elements["Γ̄ᵍ"],k_Γ)
ops[6](elements["Γ̄ᵍ"][1:1],f₂)

k_Γ = k_Γ[:,2:end]
k = [k k_Γ;k_Γ' zeros(np-1,np-1)]
f = [f;f₂[2:end]]


# ops[7](elements["Γᵍ"],k,f)

u = ops[8](elements["∂Ω"])
# x = [0.0,1.0,1.0,0.0,0.5,1.0,0.5,0.0,0.5]
# y = [0.0,0.0,1.0,1.0,0.0,0.5,1.0,0.5,0.5]
# x = [1.0,0.0,0.5,1.0,0.5,0.0,0.6]
# y = [0.0,1.0,0.0,0.5,1.0,0.5,0.6]
# v = [ū(x_,y_,0.0) for (x_,y_) in zip(x,y)]

index = [2,3,4,5]
# index = [2,3,7,9,11,12,15,16]
# index = [2,3,7,11,15,17,19,30,32,43,45,46,49,52,55,56]
for i in index
    f .-= k[:,i]*u[i]
    k[:,i] .= 0.
    k[i,:] .= 0.
    k[i,i] = 1.
    f[i] = u[i]
end

d = k\f

# println(k*[u;v]-f)
# println(d-[u;v])
# println(norm(d₀[1:nₚ₁]-u))
# r = k_Γ'*u