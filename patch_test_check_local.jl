
using ApproxOperator, CairoMakie, LinearAlgebra
include("input.jl")

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
# ū(x,y,z) = 0.5*x^2 + 2x*y + 3y^2
# ∂ū∂x(x,y,z) = x + 2y
# ∂ū∂y(x,y,z) = 2x + 6y
# b̄(x,y,z) = -7.
# ū(x,y,z) = 7x^3+8x^2*y+9x*y^2+10y^3
# ∂ū∂x(x,y,z) = 21x^2+16x*y+9y^2
# ∂ū∂y(x,y,z) = 8x^2+18x*y+30y^2
# b̄(x,y,z) = -60x-76y

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
op_debug = Operator{:∫vtdΓ_debug}()

elements, 𝓑 = import_tr("./msh/square_2.msh")
nₚ = 16
nₚ₁ = nₚ

# index = [2,3,4,5]
index = [2,3,7,9,11,12,15,16]
# index = [2,3,7,11,15,17,19,30,32,43,45,46,49,52,55,56]

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵍ"])
set∇𝝭!(elements["Γ"])
set𝝭!(elements["Γ"])

prescribe!(elements["Γᵍ"],:g=>ū)
prescribe!(elements["∂Ω"],:u=>ū)
prescribe!(elements["Ω"],:b=>b̄)
prescribe!(elements["Γ"],:g=>ū)
prescribe!(elements["Γ"],:t=>(x,y,z,n₁,n₂)->∂ū∂x(x,y,z)*n₁ + ∂ū∂y(x,y,z)*n₂)

u = ops[8](elements["∂Ω"])

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
fₜ = zeros(nₚ)
ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)
ops[3](elements["Γ"],fₜ)

for i in index
    f .-= k[:,i]*u[i]
    k[:,i] .= 0.
    k[i,:] .= 0.
    k[i,i] = 1.
    f[i] = u[i]
end

d = k\f
# println(norm(d-u)/nₚ)

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
fₜ = zeros(nₚ)
# ops[1](elements["Ω"][1:2],k)
# ops[2](elements["Ω"][1:2],f)
# ops[3](elements["Γ"][[1,2,3,4,5,7]],fₜ)
ops[1](elements["Ω"][1:1],k)
ops[2](elements["Ω"][1:1],f)
ops[3](elements["Γ"][[1,3,4]],fₜ)

println(norm(k*u-f-fₜ))
println(k*d-f)
u₁ = u
d₁ = d
k₁ = k
f₁ = f
fₜ₁ = fₜ


elements, 𝓑 = import_tr("./msh/square_4.msh")
nₚ = 56
nₚ₂ = nₚ

# index = [2,3,4,5]
# index = [2,3,7,9,11,12,15,16]
index = [2,3,7,11,15,17,19,30,32,43,45,46,49,52,55,56]

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵍ"])
set∇𝝭!(elements["Γ"])
set𝝭!(elements["Γ"])

prescribe!(elements["Γᵍ"],:g=>ū)
prescribe!(elements["∂Ω"],:u=>ū)
prescribe!(elements["Ω"],:b=>b̄)
prescribe!(elements["Γ"],:g=>ū)
prescribe!(elements["Γ"],:t=>(x,y,z,n₁,n₂)->∂ū∂x(x,y,z)*n₁ + ∂ū∂y(x,y,z)*n₂)

u = ops[8](elements["∂Ω"])

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
fₜ = zeros(nₚ)
ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)
ops[3](elements["Γ"],fₜ)

for i in index
    f .-= k[:,i]*u[i]
    k[:,i] .= 0.
    k[i,:] .= 0.
    k[i,i] = 1.
    f[i] = u[i]
end

d = k\f
# println(norm(d-u)/nₚ)

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
fₜ = zeros(nₚ)
# ops[1](elements["Ω"][[1:4...,9:12...]],k)
# ops[2](elements["Ω"][[1:4...,9:12...]],f)
# ops[3](elements["Γ"][[1:14...,29:32...,34:38...,40]],fₜ)
ops[1](elements["Ω"][[1,2,3,9]],k)
ops[2](elements["Ω"][[1,2,3,9]],f)
ops[3](elements["Γ"][[1,2,3,4,5,6,7,8,9,11,29,31]],fₜ)

println(norm(k*u-f-fₜ))
println(k*d-f)
u₂ = u
d₂ = d
k₂ = k
f₂ = f
fₜ₂ = fₜ

println((u₂[6]+u₂[18])/2 - u₁[1])

# k₁ = k₁[1,:]
# f₁ = f₁[1]
# fₜ₁ = fₜ₁[1]
# k₂ = (k₂[6,:]+k₂[18,:])/2
# f₂ = (f₂[6]+f₂[18])/2
# fₜ₂ = (fₜ₂[6]+fₜ₂[18])/2
# println(k₁'*u₁-f₁-fₜ₁)
# println(k₁'*d₁-f₁)
# println(k₂'*u₂-f₂-fₜ₂)
# println(k₂'*d₂-f₂)
# c₁_ = fₜ₁
# c₂_ = fₜ₂
# c₁ = c₁_/(c₁_-c₂_)
# c₂ = c₂_/(c₁_-c₂_)
# println(c₁)
# println((k₂'*u₂-f₂)*c₁ - (k₁'*u₁-f₁)*c₂)
# println((k₂'*d₂-f₂)*c₁ - (k₁'*d₁-f₁)*c₂)
# k = zeros(nₚ₂)
# k .= k₂*c₁
# k[2] -= k₁[2]/2*c₂
# k[3] -= k₁[3]/2*c₂
# k[6] -= k₁[1]/2*c₂
# k[7] -= k₁[2]/2*c₂
# k[9] -= k₁[5]/2*c₂
# k[18] -= k₁[1]/2*c₂
# k[19] -= k₁[3]/2*c₂
# k[20] -= k₁[4]/2*c₂
# k[23] -= k₁[4]/2*c₂
# k[24] -= k₁[5]/2*c₂
# f = f₂*c₁-f₁*c₂
# println(k'*u₂-f)
# d = zeros(nₚ₂)
# d .= d₂
# d[2] -= 
# println(k'*d₂-f)

# fₜ₁ = fₜ₁[1]
# fₜ₂ = (fₜ₂[6]+fₜ₂[18])/2
# fₜ₁ = fₜ₁[2]
# fₜ₂ = (fₜ₂[2]+fₜ₂[7])/2
# fₜ₁ = fₜ₁[3]
# fₜ₂ = (fₜ₂[3]+fₜ₂[19])/2
# c₁_ = fₜ₁
# c₂_ = fₜ₂
# c₁ = c₁_/(c₁_-c₂_)
# c₂ = c₂_/(c₁_-c₂_)
# println(c₁)

println(sum(fₜ₁)-sum(fₜ₂))
# println((fₜ₂[1]+fₜ₂[4]+fₜ₂[5]+fₜ₂[6]+fₜ₂[18])/fₜ₁[1])
# println((fₜ₂[2]+fₜ₂[7])/fₜ₁[2])
# println((fₜ₂[3]+fₜ₂[19])/fₜ₁[3])
# println(sum(fₜ₂[[5,6,7]])/sum(fₜ₁))
# println(sum(fₜ₂[[1,5,4]])/sum(fₜ₁))
# println(sum(fₜ₂[[1,2,3]])/sum(fₜ₁))
# println(sum(fₜ₂[[4,18,19]])/sum(fₜ₁))

# println((fₜ₂[1]+fₜ₂[6]+fₜ₂[18]-f₂[6]-f₂[8])/(fₜ₁[1]-f₁[1]))
println(fₜ₁[1]-fₜ₂[1]-fₜ₂[6]-fₜ₂[18])
println(fₜ₁[2]-fₜ₂[2]-fₜ₂[4]-fₜ₂[7])
println(fₜ₁[3]-fₜ₂[3]-fₜ₂[5]-fₜ₂[19])