
using ApproxOperator, CairoMakie, LinearAlgebra
include("input.jl")

elements, 𝓑 = import_tr("./msh/square_1.msh")
nₚ₁ = 5
nₚ₂ = 4
# nₚ₁ = 16
# nₚ₂ = 9

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
set∇𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γ̄ᵍ"])
set∇𝝭!(elements["Γ"])
set𝝭!(elements["Γ"])
set𝝭!(elements["Γ̄"])

# ū(x,y,z) = 1.0
# ū(x,y,z) = 1.0+2x+3y
# ∂ū∂x(x,y,z) = 2.0
# ∂ū∂y(x,y,z) = 3.0
# b̄(x,y,z) = 0.0

# ū(x,y,z) = 1.0+2x+3y+4x^2+5x*y+6y^2
# ∂ū∂x(x,y,z) = 2.0+8x+5y
# ∂ū∂y(x,y,z) = 3.0+5x+12y
# b̄(x,y,z) = -20.0

# ū(x,y,z) = 4x^2+5x*y+6y^2
# ∂ū∂x(x,y,z) = 8x+5y
# ∂ū∂y(x,y,z) = 5x+12y
# b̄(x,y,z) = -20.0

ū(x,y,z) = 2x^2+1x*y+1y^2
∂ū∂x(x,y,z) = 4x+1y
∂ū∂y(x,y,z) = 1x+2y
b̄(x,y,z) = -6.0

prescribe!(elements["Γᵍ"],:g=>ū)
prescribe!(elements["∂Ω"],:u=>ū)
# prescribe!(elements["Γᵍ"],:g=>(x,y,z)->1.0+2x+3y)
prescribe!(elements["Ω"],:b=>b̄)
prescribe!(elements["Γ"],:t=>(x,y,z,n₁,n₂)->∂ū∂x(x,y,z)*n₁ + ∂ū∂y(x,y,z)*n₂)
prescribe!(elements["Γᵍ"],:t=>(x,y,z,n₁,n₂)->∂ū∂x(x,y,z)*n₁ + ∂ū∂y(x,y,z)*n₂)

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
f = zeros(nₚ₁)
d = zeros(nₚ₁)
k₀ = zeros(nₚ₁,nₚ₁)
fint = zeros(nₚ₁)
f_Γ = zeros(nₚ₁)
f_Ω = zeros(nₚ₁)
Δd = zeros(nₚ₁)
fext = zeros(nₚ₁)
ftest = zeros(nₚ₁)

ops[1](elements["Ω"],k₀)
ops[2](elements["Ω"],f_Ω)
ops[2](elements["Ω"],fext)
ops[3](elements["Γ"],fext)
ops[3](elements["Γᵍ"],ftest)

u = ops[8](elements["∂Ω"])

index = [2,3,4,5]
# index = [2,3,7,9,11,12,15,16]
d[index] .= u[index]

# println(k₀*u-fext)

# println("iter 1")
# println(Δd)
# println(f_Ω - k₀*d)
# println(d - u)

# println(d - u)

# tol = 1e-10
# normΔd = 1.0
# maxiter = 10

# iter = 0
# while normΔd≥tol && iter≤maxiter
#     global iter += 1
    
#     k .= k₀
#     fint .= k₀*d
#     f .= f_Ω - fint
#     if iter ≠ 1
#         f[1] -= (f[2]+f[3]+f[4]+f[5])
#     end

#     for i in index
#         k[:,i] .= 0.
#         k[i,:] .= 0.
#         k[i,i] = 1.
#         f[i] = 0.0
#     end

#     Δd .= k\f
#     d .+= Δd
#     global normΔd = norm(Δd)
#     println("iter = $iter, norm = $normΔd")
#     println(Δd)
#     println(f_Ω - k₀*d)
#     println(d - u)
# end


# r = k_Γ'*u


k .= k₀
fint .= k₀*d
f .= f_Ω - fint
f[1] = -3.5 - fint[1]
# f[1] = -0.541666666666643 - fint[1]
# # f[1] = -4.04166666666666666666 - fint[1]
# f[1] = -10 - fint[1]
# f[1] = -11.874999999999995 - fint[1]

for i in index
    # f .-= k[:,i]*u[i]
    k[:,i] .= 0.
    k[i,:] .= 0.
    k[i,i] = 1.
    f[i] = 0.0
end

Δd = k\f
# println(k₀*d-fext)
d .+= Δd
# println(f)
# println(d)
println(f_Ω - k₀*d)
println(sum((f_Ω - k₀*d)[2:5])/2)
println(d-u)

# k .= k₀
# fint .= k₀*d
# f .= f_Ω - fint
# f[1] = sum((f_Ω - k₀*d)[2:5])/2 - fint[1]

# for i in index
#     k[:,i] .= 0.
#     k[i,:] .= 0.
#     k[i,i] = 1.
#     f[i] = 0.0
# end

# Δd = k\f
# d .+= Δd
# println(d-u)