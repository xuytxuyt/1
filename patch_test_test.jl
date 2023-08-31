
using ApproxOperator, CairoMakie
include("input.jl")

elements = import_ct("./msh/square_2.msh")
# nₚ₁ = 4
# nₚ₂ = 4
nₚ₁ = 16
nₚ₂ = 9

# plot and check elements
# f = Figure()
# Axis(f[1,1], autolimitaspect = 1)

# xs = Float64[]
# ys = Float64[]
# xg = Float64[]
# yg = Float64[]
# for a in elements["Ω"]
#     𝓒 = collect(a.𝓒)[[1,2,3,1]]
#     𝓖 = a.𝓖
#     for xᵢ in 𝓒
#         push!(xs,xᵢ.x)
#         push!(ys,xᵢ.y)
#     end
#     for ξ in 𝓖
#         push!(xg,ξ.x)
#         push!(yg,ξ.y)
#     end
# end
# scatterlines!(xs,ys,width=2,color=:black,markercolor=:black,marksize=20)
# scatter!(xg,yg,marker=:xcross,marksize=20,color=:red)
# xg = Float64[]
# yg = Float64[]
# for a in elements["Γᵍ"]
#     𝓖 = a.𝓖
#     for ξ in 𝓖
#         push!(xg,ξ.x)
#         push!(yg,ξ.y)
#     end
# end
# scatter!(xg,yg,marker=:cross,marksize=20,color=:blue)
# f

set𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γ"])

# ū(x,y,z) = 1.0
# ū(x,y,z) = 1.0+2x+3y
ū(x,y,z) = 1.0+2x+3y+4x^2+5x*y+6y^2
b̄(x,y,z) = -20.0
# ū(x,y,z) = y^2
# b̄(x,y,z) = -2.0
# ū(x,y,z) = x
# b̄(x,y,z) = 0.0
prescribe!(elements["Γ̄ᵍ"],:g=>ū)
# prescribe!(elements["∂Ω"],:u=>ū)
# prescribe!(elements["Γᵍ"],:g=>(x,y,z)->1.0+2x+3y)
prescribe!(elements["Ω"],:b=>b̄)

ops = [
    Operator{:∫∫∇v∇udxdy}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫vᵢnᵢuds}(:k=>-1.0),
    Operator{:∫vᵢnᵢuds}(:k=>1.0),
    Operator{:∫vᵢnᵢgds}(:k=>1.0),
    Operator{:∫vgdΓ}(:α=>1e7),
    Operator{:∫uds}(),
]

k = zeros(nₚ₁,nₚ₁)
k_Γ = zeros(nₚ₁,nₚ₂)
# k_Γ₁ = zeros(nₚ₁,nₚ₂)
# k_Γ₂ = zeros(nₚ₁,nₚ₂)
f₂ = zeros(nₚ₂)
f = zeros(nₚ₁)

ops[1](elements["Ω̄"],k)
ops[2](elements["Ω"],f₂)
# ops[4](elements["Γ"][1:1],elements["Γ̄"][1:1],k_Γ₁)
# ops[4](elements["Γ"][2:2],elements["Γ̄"][2:2],k_Γ₂)
ops[4](elements["Γ̄"],elements["Γ"],k_Γ)
ops[5](elements["Γ̄ᵍ"],elements["Γᵍ"],k_Γ)
# ops[6](elements["Γ̄ᵍ"][6:6],f)
ops[6](elements["Γ̄ᵍ"],f)

k_Γ = k_Γ[:,[2,4:nₚ₂...]]
k = [k k_Γ;k_Γ' zeros(nₚ₂-2,nₚ₂-2)]
f = [f;-f₂[[2,4:nₚ₂...]]]

# k = [k k_Γ;k_Γ' zeros(nₚ₂,nₚ₂)]
# f = [f;-f₂]

d = k\f
