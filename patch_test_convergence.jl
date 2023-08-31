
using ApproxOperator, CairoMakie, LinearAlgebra
include("input.jl")

elements, 𝓑 = import_tr("./msh/square_2.msh")
# elements, 𝓑 = import_tr("./msh/square_2_irregular.msh")
# nₚ₁ = 5
# nₚ₂ = 4
nₚ₁ = 16
nₚ₂ = 9
# nₚ₁ = 56
# nₚ₂ = 25

index = [2,3,4,5]
# index = [2,3,7,9,11,12,15,16]
# index = [2,3,7,11,15,17,19,30,32,43,45,46,49,52,55,56]

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

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ω"])
set𝝭!(elements["Ω̄"])
set∇𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γ̄ᵍ"])
set∇𝝭!(elements["Γ"])
set𝝭!(elements["Γ"])
set𝝭!(elements["Γ̄"])

# ū(x,y,z) = 1.0
# ū(x,y,z) = 1.0+2x+3y
# ∂ū∂x(x,y,z) = 2.
# ∂ū∂y(x,y,z) = 3.
# b̄(x,y,z) = 0.0
ū(x,y,z) = 1.0+2x+3y+4x^2+5x*y+6y^2+7x^3+8x^2*y+9x*y^2+10y^3
∂ū∂x(x,y,z) = 2.0+8x+5y+21x^2+16x*y+9y^2
∂ū∂y(x,y,z) = 3.0+5x+12y+8x^2+18x*y+30y^2
b̄(x,y,z) = -20.0-60x-76y
# ū(x,y,z) = 10x^2+9x*y+8y^2
# ∂ū∂x(x,y,z) = 20x+9y
# ∂ū∂y(x,y,z) = 9x+16y
# b̄(x,y,z) = -36.0
# ū(x,y,z) = 0.5*x^2 + x*y
# ∂ū∂x(x,y,z) = x + y
# ∂ū∂y(x,y,z) = x
# b̄(x,y,z) = -1.
prescribe!(elements["Γᵍ"],:g=>ū)
prescribe!(elements["∂Ω"],:u=>ū)
# prescribe!(elements["Γᵍ"],:g=>(x,y,z)->1.0+2x+3y)
prescribe!(elements["Ω̄"],:b=>b̄)
prescribe!(elements["Ω"],:b=>b̄)
prescribe!(elements["Γ̄"],:t=>(x,y,z,n₁,n₂)->∂ū∂x(x,y,z)*n₁ + ∂ū∂y(x,y,z)*n₂)
prescribe!(elements["Γ"],:g=>ū)
prescribe!(elements["Γ"],:t=>(x,y,z,n₁,n₂)->∂ū∂x(x,y,z)*n₁ + ∂ū∂y(x,y,z)*n₂)

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
k_Γ₁ = zeros(nₚ₁,nₚ₂)
k_Γ₂ = zeros(nₚ₁,nₚ₂)
f₂ = zeros(nₚ₂)
f = zeros(nₚ₁)

ops[1](elements["Ω"],k)
ops[2](elements["Ω̄"],f₂)
# ops[4](elements["Γ"][1:1],elements["Γ̄"][1:1],k_Γ₁)
# ops[4](elements["Γ"][2:2],elements["Γ̄"][2:2],k_Γ₂)
ops[4](elements["Γ"],elements["Γ̄"],k_Γ)
ops[5](elements["Γᵍ"],elements["Γ̄ᵍ"],k_Γ)
ops[6](elements["Γᵍ"],f)

k_Γ = k_Γ[:,[2,4:nₚ₂...]]
k = [k k_Γ;k_Γ' zeros(nₚ₂-2,nₚ₂-2)]
f = [f;-f₂[[2,4:nₚ₂...]]]


# ops[7](elements["Γᵍ"],k,f)

u = ops[8](elements["∂Ω"])
# x = [0.0,1.0,1.0,0.0,0.5,1.0,0.5,0.0,0.5]
# y = [0.0,0.0,1.0,1.0,0.0,0.5,1.0,0.5,0.5]
# x = [1.0,0.0,0.5,1.0,0.5,0.0,0.6]
# y = [0.0,1.0,0.0,0.5,1.0,0.5,0.6]
# v = [ū(x_,y_,0.0) for (x_,y_) in zip(x,y)]

for i in index
    f .-= k[:,i]*u[i]
    k[:,i] .= 0.
    k[i,:] .= 0.
    k[i,i] = 1.
    f[i] = u[i]
end

# d₀ = k\f

# println(k*[u;v]-f)
# println(d-[u;v])
# println(norm(d₀[1:nₚ₁]-u))
# r = k_Γ'*u


# k = zeros(nₚ₁,nₚ₁)
# k_Γ = zeros(nₚ₁,nₚ₂)
# f₁ = zeros(nₚ₁)
# f₂ = zeros(nₚ₂)
# x = [0.0,1.0,1.0,0.0,0.5,1.0,0.5,0.0,0.6]
# y = [0.0,0.0,1.0,1.0,0.0,0.5,1.0,0.5,0.6]
# v = [ū(x_,y_,0.0) for (x_,y_) in zip(x,y)]
# ops[1](elements["Ω"],k)
# ops[2](elements["Ω̄"],f₂)
# ops[3](elements["Γ̄"],f₂)
# ops[4](elements["Γ"],elements["Γ̄"],k_Γ)
# ops[6](elements["Γ"],f₁)
# println(k_Γ'*u)
# println(f₂)
# println(k_Γ'*u+f₂)
# println(k_Γ*v)
# println(f₁)
# println(k_Γ*v+f₁)
# println(k*u+k_Γ*v)


# ic = 9
# k_Γ = zeros(nₚ₁,nₚ₂)
# f₁ = zeros(nₚ₁)
# ops[4](elements["Γ"][ic:ic],elements["Γ̄"][ic:ic],k_Γ)
# ops[6](elements["Γ"][ic:ic],f₁)
# println(k_Γ)
# println(v)
# println(k_Γ*v)
# println(f₁)
# println(k_Γ*v+f₁)

k = zeros(nₚ₁,nₚ₁)
f = zeros(nₚ₁)
fₜ = zeros(nₚ₁)
ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)
ops[3](elements["Γ"],fₜ)
# index = [2,3,4,5]
# index = [2,3,7,9,11,12,15,16]
# index = [2,3,7,11,15,17,19,30,32,43,45,46,49,52,55,56]
for i in index
    f .-= k[:,i]*u[i]
    k[:,i] .= 0.
    k[i,:] .= 0.
    k[i,i] = 1.
    f[i] = u[i]
end

d₁ = k\f
println(norm(d₁-u)/nₚ₁)

k = zeros(nₚ₁,nₚ₁)
f = zeros(nₚ₁)
ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)
ops[3](elements["Γ"],f)
for i in index
    f .-= k[:,i]*u[i]
    k[:,i] .= 0.
    k[i,:] .= 0.
    k[i,i] = 1.
    f[i] = u[i]
end

d₂ = k\f
println(norm(d₂-u)/nₚ₁)

# f = zeros(nₚ₁)
# ops[3](elements["Γ"][7:7],f)
# ops[3](elements["Γ"],f)
# ops[2](elements["Ω"],f)