
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
# op_debug(elements["Γ"][3:3],fₜ)

for i in index
    f .-= k[:,i]*u[i]
    k[:,i] .= 0.
    k[i,:] .= 0.
    k[i,i] = 1.
    f[i] = u[i]
end

d₁ = k\f
# println(norm(d₁-u)/nₚ)
k₁ = k
invk₁ = inv(k)
fₜ₁ = fₜ
fₜ₁[index] .= 0.0
u₁ = u
a = elements["Γ"][3]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
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

d₀ = k\f
# println(norm(d₀-u)/nₚ₁)


elements, 𝓑 = import_tr("./msh/square_4.msh")
# nₚ = 5
# nₚ = 16
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

d₂ = k\f
# println(norm(d₂-u)/nₚ)
k₂ = k
invk₂ = inv(k)
fₜ₂ = fₜ
fₜ₂[index] .= 0.
u₂ = u

# invk₂₁ = zeros(nₚ₁,nₚ₂)
# invk₂₂ = zeros(nₚ₁,nₚ₁)
# fc = zeros(nₚ₁)
checkf = zeros(nₚ₁)

# invk₂₁[1,:] .= (invk₂[6,:] + invk₂[10,:])/2
# invk₂₁[2,:] .= (invk₂[7,:] + invk₂[2,:])/2
# invk₂₁[3,:] .= (invk₂[3,:] + invk₂[11,:])/2
# invk₂₁[4,:] .= (invk₂[12,:] + invk₂[15,:])/2
# invk₂₁[5,:] .= (invk₂[9,:] + invk₂[16,:])/2

# invk₂₂[:,1] .= (invk₂₁[:,1] + invk₂₁[:,6] + invk₂₁[:,10] + invk₂₁[:,14])
# invk₂₂[:,2] .= (invk₂₁[:,2] + invk₂₁[:,4] + invk₂₁[:,7])
# invk₂₂[:,3] .= (invk₂₁[:,3] + invk₂₁[:,11] + invk₂₁[:,5])
# invk₂₂[:,4] .= (invk₂₁[:,8] + invk₂₁[:,12] + invk₂₁[:,15])
# invk₂₂[:,5] .= (invk₂₁[:,9] + invk₂₁[:,13] + invk₂₁[:,16])

# c₁ = f₁[1]
# c₂ = (f₂[6]+f₂[10])/2
# d = (d₂[6] + d₂[10])/2*c₁/(c₁-c₂) - d₁[1]*c₂/(c₁-c₂)
# println(d)

checkf[1] = fₜ₁[1] - fₜ₂[1] - fₜ₂[6] - fₜ₂[18] - fₜ₂[22]
checkf[2] = fₜ₁[2] - fₜ₂[2] - fₜ₂[4] - fₜ₂[7]
checkf[3] = fₜ₁[3] - fₜ₂[3] - fₜ₂[19] - fₜ₂[5]
checkf[4] = fₜ₁[4] - fₜ₂[8] - fₜ₂[20] - fₜ₂[23] - fₜ₂[33]
checkf[5] = fₜ₁[5] - fₜ₂[9] - fₜ₂[13] - fₜ₂[21] - fₜ₂[24]
checkf[6] = fₜ₁[6] - fₜ₂[10] - fₜ₂[14] - fₜ₂[25] - fₜ₂[28]
checkf[7] = fₜ₁[7] - fₜ₂[11] - fₜ₂[12] - fₜ₂[15]
checkf[8] = fₜ₁[8] - fₜ₂[16] - fₜ₂[26] - fₜ₂[29] - fₜ₂[39]
checkf[9] = fₜ₁[9] - fₜ₂[17] - fₜ₂[30] - fₜ₂[27]
checkf[10] = fₜ₁[10] - fₜ₂[31] - fₜ₂[35] - fₜ₂[44] - fₜ₂[48]
checkf[11] = fₜ₁[11] - fₜ₂[32] - fₜ₂[34] - fₜ₂[45]
checkf[12] = fₜ₁[12] - fₜ₂[36] - fₜ₂[46] - fₜ₂[49]
checkf[13] = fₜ₁[13] - fₜ₂[37] - fₜ₂[40] - fₜ₂[47] - fₜ₂[50]
checkf[14] = fₜ₁[14] - fₜ₂[38] - fₜ₂[41] - fₜ₂[51] - fₜ₂[54]
checkf[15] = fₜ₁[15] - fₜ₂[42] - fₜ₂[52] - fₜ₂[55]
checkf[16] = fₜ₁[16] - fₜ₂[43] - fₜ₂[53] - fₜ₂[56]

f₁ = invk₁*fₜ₁
f₂ = invk₂*fₜ₂
println(norm(u₁-(d₁+f₁)))
println(norm(u₂-(d₂+f₂)))

c₁ = f₁[1]
c₂ = (f₂[6]+f₂[18])/2
d = (d₂[6] + d₂[18])/2*c₁/(c₁-c₂) - d₁[1]*c₂/(c₁-c₂)
println(u₁[1]-d)
println(c₁/(c₁-c₂))

c₁ = f₁[4]
c₂ = (f₂[20]+f₂[23])/2
d = (d₂[20] + d₂[23])/2*c₁/(c₁-c₂) - d₁[4]*c₂/(c₁-c₂)
println(u₁[4]-d)
println(c₁/(c₁-c₂))

c₁ = f₁[5]
c₂ = (f₂[9]+f₂[24])/2
d = (d₂[9] + d₂[24])/2*c₁/(c₁-c₂) - d₁[5]*c₂/(c₁-c₂)
println(u₁[5]-d)
println(c₁/(c₁-c₂))

c₁ = f₁[6]
c₂ = (f₂[14]+f₂[25])/2
d = (d₂[14] + d₂[25])/2*c₁/(c₁-c₂) - d₁[6]*c₂/(c₁-c₂)
println(u₁[6]-d)
println(c₁/(c₁-c₂))

c₁ = f₁[8]
c₂ = (f₂[26]+f₂[29])/2
d = (d₂[26] + d₂[29])/2*c₁/(c₁-c₂) - d₁[8]*c₂/(c₁-c₂)
println(u₁[8]-d)
println(c₁/(c₁-c₂))

c₁ = f₁[10]
c₂ = (f₂[35]+f₂[44])/2
d = (d₂[35] + d₂[44])/2*c₁/(c₁-c₂) - d₁[10]*c₂/(c₁-c₂)
println(u₁[10]-d)
println(c₁/(c₁-c₂))

c₁ = f₁[13]
c₂ = (f₂[37]+f₂[50])/2
d = (d₂[37] + d₂[50])/2*c₁/(c₁-c₂) - d₁[13]*c₂/(c₁-c₂)
println(u₁[13]-d)
println(c₁/(c₁-c₂))

c₁ = f₁[14]
c₂ = (f₂[41]+f₂[51])/2
d = (d₂[41] + d₂[51])/2*c₁/(c₁-c₂) - d₁[14]*c₂/(c₁-c₂)
println(u₁[14]-d)
println(c₁/(c₁-c₂))