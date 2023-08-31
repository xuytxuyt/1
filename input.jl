
function import_tri3(filename::String)
    elms,~ = ApproxOperator.importmsh(filename)
    nₚ = length(elms["Ω"][1].x)

    d = zeros(nₚ)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Tri3},:TriGI3)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])
    push!(f_Ω,
        :d=>(:𝐼,d),
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
    )
    push!(f_Γᵍ,
        :𝝭=>:𝑠,
    )
    if haskey(elms,"Γᵗ")
        f_Γᵗ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2)
        elements["Γᵗ"] = f_Γᵗ(elms["Γᵗ"])
        push!(f_Γᵗ,
            :𝝭=>:𝑠,
        )
    end
    return elements, d
end

function import_tr(filename::String)
    elms,~ = ApproxOperator.importmsh(filename)
    elms_ΓᵍΩ = elms["Γᵍ"]∩elms["Ω"]
    
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    𝓑 = ApproxOperator.getboundaries(elms["Ω"])
    ndof = length(𝓑)
    d = zeros(ndof)
    d̄ = zeros(length(elms["Ω"][1].x))
    𝐴 = [ApproxOperator.get𝐴(a) for a in elms["Ω"]]
    𝐴_Γ = [ApproxOperator.get𝐴(b) for (a,b) in elms_ΓᵍΩ]

    f_Ω = ApproxOperator.Field{(:𝐼,:𝐽),2,(:𝑔,:𝐺,:𝐶,:𝑠),4}(TRElement{:Tri3},:TriGI13)
    f_Ω̄ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Tri3},:TriGI13)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,:𝐽),2,(:𝑔,:𝐺,:𝐶,:𝑠),4}(TRElement{:Tri3},:SegGI5)
    elements["Ω"] = f_Ω(elms["Ω"],𝓑)
    elements["Ω̄"] = f_Ω̄(elms["Ω"])
    elements["Γᵍ"] = f_Γᵍ(elms_ΓᵍΩ,𝓑)
    push!(f_Ω,
        :d=>(:𝐼,d),
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝐴=>(:𝐶,𝐴),
    )
    push!(f_Ω̄,
        :d=>(:𝐼,d̄),
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝐴=>(:𝐶,𝐴),
    )
    push!(f_Γᵍ,
        :d=>(:𝐼,d),
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝐴=>(:𝐶,𝐴_Γ),
    )

    elms_Γ̄ᵍΩ = [elm[1] for elm in elms_ΓᵍΩ]
    ndof = length(elms["Ω"][1].x)
    d = zeros(ndof)
    f_Ω̄ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Tri3},:TriGI13)
    f_Γ̄ᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI5)
    elements["Ω̄"] = f_Ω̄(elms["Ω"])
    elements["Γ̄ᵍ"] = f_Γ̄ᵍ(elms_Γ̄ᵍΩ)
    push!(f_Ω̄,
        :d=>(:𝐼,d),
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝐴=>(:𝐶,𝐴),
    )
    push!(f_Γ̄ᵍ,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝐴=>(:𝐶,𝐴_Γ),
    )

    x = elms["Ω"][1].x
    y = elms["Ω"][1].y
    z = elms["Ω"][1].z
    f_∂Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI5)
    elms["∂Ω"] = [ApproxOperator.Seg2(Tuple(s),x,y,z) for s in 𝓑]
    elements["∂Ω"] = f_∂Ω(elms["∂Ω"])
    𝐿 = [ApproxOperator.get𝐿(elm) for elm in elms["∂Ω"]]
    push!(f_∂Ω,
        :𝐿=>(:𝐶,𝐿),
    )
    f_Γ = ApproxOperator.Field{(:𝐼,:𝐽),2,(:𝑔,:𝐺,:𝐶,:𝑠),4}(TRElement{:Tri3},:SegGI5)
    elms_∂ΩΩ = elms["∂Ω"]∩elms["Ω"]
    𝐴 = [ApproxOperator.get𝐴(elm[2]) for elm in elms_∂ΩΩ]
    𝐿 = [ApproxOperator.get𝐿(elm[1]) for elm in elms_∂ΩΩ]
    elements["Γ"] = f_Γ(elms_∂ΩΩ,𝓑)
    elms_Γ = [elm[1] for elm in elms_∂ΩΩ]
    n₁ = zeros(length(elms_Γ))
    n₂ = zeros(length(elms_Γ))
    for (i,a) in enumerate(elms_Γ)
        x₁ = a.x[a.i[1]]
        y₁ = a.y[a.i[1]]
        x₂ = a.x[a.i[2]]
        y₂ = a.y[a.i[2]]
        n₁[i] = (y₂-y₁)/𝐿[i]
        n₂[i] = (x₁-x₂)/𝐿[i]
    end
    f_Γ̄ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI5)
    elements["Γ̄"] = f_Γ̄(elms_Γ)
    push!(f_Γ,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝐴=>(:𝐶,𝐴),
        :n₁=>(:𝐶,n₁),
        :n₂=>(:𝐶,n₂),
    )
    push!(f_Γ̄,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝐿=>(:𝐶,𝐿),
        :n₁=>(:𝐶,n₁),
        :n₂=>(:𝐶,n₂),
    )
    return elements, 𝓑
end

function import_ct(filename::String)
    elms,~ = ApproxOperator.importmsh(filename)
    nₚ = length(elms["Ω"][1].x)

    d = zeros(nₚ)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Tri3},:TriGI3)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])
    push!(f_Ω,
        :d=>(:𝐼,d),
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
    )
    push!(f_Γᵍ,
        :𝝭=>:𝑠,
    )
    if haskey(elms,"Γᵗ")
        f_Γᵗ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2)
        elements["Γᵗ"] = f_Γᵗ(elms["Γᵗ"])
        push!(f_Γᵗ,
            :𝝭=>:𝑠,
        )
    end

    elms["Γ"] = ApproxOperator.Seg2[]
    for elm in elms["Ω"]
        for i in zip(elm.i[[1,2,3]],elm.i[[2,3,1]])
            push!(elms["Γ"],ApproxOperator.Seg2(i,elm.x,elm.y,elm.z))
        end
    end
    f_Γ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2)
    elements["Γ"] = f_Γ(elms["Γ"])
    push!(f_Γ,
        :𝝭=>:𝑠,
    )

    elements["Ω̄"] = Element{:Tri3}[]
    elements["Γ̄"] = Element{:Seg2}[]
    nₑ = length(elms["Ω"])
    B₁ = zeros(2*nₑ)
    # B₁[1:2:2*nₑ] .= 1.0
    B₂ = zeros(2*nₑ)
    # B₂[2:2:2*nₑ] .= 1.0
    n₁ = zeros(6*nₑ)
    n₂ = zeros(6*nₑ)
    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    data𝓒 = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    data𝓖 = Dict{Symbol,Tuple{Int,Vector{Float64}}}([
        :∂𝝭∂x=>(4,B₁),
        :∂𝝭∂y=>(4,B₂),
        :𝑤=>(3,[ApproxOperator.get𝐴(elm) for elm in elms["Ω"]])
    ])
    𝓖_Γ = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    data𝓖_Γ = Dict{Symbol,Tuple{Int,Vector{Float64}}}([
        :∂𝝭∂x=>(4,B₁),
        :∂𝝭∂y=>(4,B₂),
        :𝑤=>(3,[ApproxOperator.get𝐴(elm) for elm in elms["Ω"]]),
        :x=>(2, getfield(getfield(elements["Γ"][1],:𝓖)[3][1],:data)[:x][2]),
        :y=>(2, getfield(getfield(elements["Γ"][1],:𝓖)[3][1],:data)[:y][2]),
        :z=>(2, getfield(getfield(elements["Γ"][1],:𝓖)[3][1],:data)[:z][2]),
        :𝑤=>(2, getfield(getfield(elements["Γ"][1],:𝓖)[3][1],:data)[:𝑤][2]),
        :n₁=>(2, n₁),
        :n₂=>(2, n₂),
    ])
    c = 0
    g = 0
    g_Γ = 0
    𝐼 = 0
    𝐺 = 0
    𝑠 = 0
    𝐺_Γ = 0
    for (𝐶,a) in enumerate(elms["Ω"])
        𝐼 += 1
        push!(𝓒,Node{(:𝐼,),1}((𝐼,),data𝓒))
        𝐼 += 1
        push!(𝓒,Node{(:𝐼,),1}((𝐼,),data𝓒))

        𝐺 += 1
        push!(𝓖,Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((1,𝐺,𝐶,𝑠),data𝓖))

        push!(elements["Ω̄"],Element{:Tri3}((c,2,𝓒),(g,1,𝓖)))

        x₁ = a.x[a.i[1]]
        x₂ = a.x[a.i[2]]
        x₃ = a.x[a.i[3]]
        y₁ = a.y[a.i[1]]
        y₂ = a.y[a.i[2]]
        y₃ = a.y[a.i[3]]

        𝐴 = ApproxOperator.get𝐴(a)
        B₁[2*𝐶-1] = (y₂-y₃)/2/𝐴
        B₁[2*𝐶] = (y₃-y₁)/2/𝐴
        B₂[2*𝐶-1] = (x₃-x₂)/2/𝐴
        B₂[2*𝐶] = (x₁-x₃)/2/𝐴

        𝐿 = ((x₁-x₂)^2+(y₁-y₂)^2)^0.5
        𝐺_Γ += 1
        n₁[𝐺_Γ] = (y₂-y₁)/𝐿
        n₂[𝐺_Γ] = (x₁-x₂)/𝐿
        push!(𝓖_Γ,Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((1,𝐺_Γ,𝐶,𝑠),data𝓖_Γ))
        𝐺_Γ += 1
        n₁[𝐺_Γ] = (y₂-y₁)/𝐿
        n₂[𝐺_Γ] = (x₁-x₂)/𝐿
        push!(𝓖_Γ,Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((1,𝐺_Γ,𝐶,𝑠),data𝓖_Γ))
        push!(elements["Γ̄"],Element{:Seg2}((c,2,𝓒),(g_Γ,2,𝓖_Γ)))
        g_Γ += 2

        𝐿 = ((x₂-x₃)^2+(y₂-y₃)^2)^0.5
        𝐺_Γ += 1
        n₁[𝐺_Γ] = (y₃-y₂)/𝐿
        n₂[𝐺_Γ] = (x₂-x₃)/𝐿
        push!(𝓖_Γ,Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((1,𝐺_Γ,𝐶,𝑠),data𝓖_Γ))
        𝐺_Γ += 1
        n₁[𝐺_Γ] = (y₃-y₂)/𝐿
        n₂[𝐺_Γ] = (x₂-x₃)/𝐿
        push!(𝓖_Γ,Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((1,𝐺_Γ,𝐶,𝑠),data𝓖_Γ))
        push!(elements["Γ̄"],Element{:Seg2}((c,2,𝓒),(g_Γ,2,𝓖_Γ)))
        g_Γ += 2

        𝐿 = ((x₃-x₁)^2+(y₃-y₁)^2)^0.5
        𝐺_Γ += 1
        n₁[𝐺_Γ] = (y₁-y₃)/𝐿
        n₂[𝐺_Γ] = (x₃-x₁)/𝐿
        push!(𝓖_Γ,Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((1,𝐺_Γ,𝐶,𝑠),data𝓖_Γ))
        𝐺_Γ += 1
        n₁[𝐺_Γ] = (y₁-y₃)/𝐿
        n₂[𝐺_Γ] = (x₃-x₁)/𝐿
        push!(𝓖_Γ,Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((1,𝐺_Γ,𝐶,𝑠),data𝓖_Γ))
        push!(elements["Γ̄"],Element{:Seg2}((c,2,𝓒),(g_Γ,2,𝓖_Γ)))
        g_Γ += 2

        𝑠 += 2
        c += 2
        g += 1
    end

    elements["Γ̄ᵍ"] = Element{:Seg2}[]
    for a in elms["Γᵍ"]
        i = findfirst(x->x==a,elms["Γ"])
        push!(elements["Γ̄ᵍ"],elements["Γ̄"][i])
    end

    return elements
end

function import_trc(filename::String)
    elms,~ = ApproxOperator.importmsh(filename)
    
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    𝓑 = ApproxOperator.getboundaries(elms["Ω"])
    ndof = length(𝓑)
    d = zeros(ndof)
    𝐴 = [ApproxOperator.get𝐴(a) for a in elms["Ω"]]

    f_Ω = ApproxOperator.Field{(:𝐼,:𝐽),2,(:𝑔,:𝐺,:𝐶,:𝑠),4}(TRElement{:Tri3},:TriGI3)
    elements["Ω"] = f_Ω(elms["Ω"],𝓑)
    push!(f_Ω,
        :d=>(:𝐼,d),
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝐴=>(:𝐶,𝐴),
    )

    x = elms["Ω"][1].x
    y = elms["Ω"][1].y
    z = elms["Ω"][1].z
    f_∂Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2)
    elms["∂Ω"] = [ApproxOperator.Seg2(Tuple(s),x,y,z) for s in 𝓑]
    elements["∂Ω"] = f_∂Ω(elms["∂Ω"])
    𝐿 = [ApproxOperator.get𝐿(elm) for elm in elms["∂Ω"]]
    push!(f_∂Ω,
        :𝐿=>(:𝐶,𝐿),
    )
    f_Γ = ApproxOperator.Field{(:𝐼,:𝐽),2,(:𝑔,:𝐺,:𝐶,:𝑠),4}(TRElement{:Tri3},:SegGI2)
    elms_∂ΩΩ = elms["∂Ω"]∩elms["Ω"]
    𝐴 = [ApproxOperator.get𝐴(elm[2]) for elm in elms_∂ΩΩ]
    𝐿 = [ApproxOperator.get𝐿(elm[1]) for elm in elms_∂ΩΩ]
    elements["Γ"] = f_Γ(elms_∂ΩΩ,𝓑)
    elms_Γ = [elm[1] for elm in elms_∂ΩΩ]
    n₁ = zeros(length(elms_Γ))
    n₂ = zeros(length(elms_Γ))
    for (i,a) in enumerate(elms_Γ)
        x₁ = a.x[a.i[1]]
        y₁ = a.y[a.i[1]]
        x₂ = a.x[a.i[2]]
        y₂ = a.y[a.i[2]]
        n₁[i] = (y₂-y₁)/𝐿[i]
        n₂[i] = (x₁-x₂)/𝐿[i]
    end
    push!(f_Γ,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝐴=>(:𝐶,𝐴),
        :n₁=>(:𝐶,n₁),
        :n₂=>(:𝐶,n₂),
    )

    scheme = ApproxOperator.quadraturerule(:SegGI2)
    nₑ = length(elms_Γ)
    elements["Γ̄"] = TRElement{:Seg2}[]
    𝓒 = Node{(:𝐼,:𝐽),2}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    data𝓒 = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    n₁ = zeros(nₑ)
    n₂ = zeros(nₑ)
    ∂𝝭∂x = zeros(2*nₑ)
    ∂𝝭∂y = zeros(2*nₑ)
    data𝓖 = Dict([
        :w=>(1,scheme[:w]),
        :ξ=>(1,scheme[:ξ]),
        :𝝭=>(1,scheme[:ξ]),
        :∂𝝭∂x=>(4,∂𝝭∂x),
        :∂𝝭∂y=>(4,∂𝝭∂y),
        :n₁=>(3,n₁),
        :n₂=>(3,n₂),
        :x=>(2,zeros(2*nₑ)),
        :y=>(2,zeros(2*nₑ)),
        :z=>(2,zeros(2*nₑ)),
        :𝑤=>(2,zeros(2*nₑ))
    ])
    G = 0
    for (C,a) in enumerate(elms_Γ)
        push!(elements["Γ̄"],TRElement{:Seg2}((C-1,1,𝓒),(G,2,𝓖)))
        x₁ = a.x[a.i[1]]
        x₂ = a.x[a.i[2]]
        y₁ = a.y[a.i[1]]
        y₂ = a.y[a.i[2]]
        𝐿 = ((x₁-x₂)^2+(y₁-y₂)^2)^0.5
        n₁[C] = (y₂-y₁)/𝐿
        n₂[C] = (x₁-x₂)/𝐿
        I = findfirst(x->x==Set(a.i), 𝓑)
        push!(𝓒,Node{(:𝐼,:𝐽),2}((I,0),data𝓒))
        G += 1
        xg = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((1,G,C,G-1),data𝓖)
        ξ = xg.ξ
        ∂𝝭∂x[G] = n₁[C]*ξ
        ∂𝝭∂y[G] = n₂[C]*ξ
        N₁ = (1-ξ)/2
        N₂ = (1+ξ)/2
        xg.x = N₁*x₁+N₂*x₂
        xg.y = N₁*y₁+N₂*y₂
        xg.𝑤 = xg.w/2*𝐿
        push!(𝓖,xg)
        G += 1
        xg = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((2,G,C,G-1),data𝓖)
        ξ = xg.ξ
        ∂𝝭∂x[G] = n₁[C]*ξ
        ∂𝝭∂y[G] = n₂[C]*ξ
        N₁ = (1-ξ)/2
        N₂ = (1+ξ)/2
        xg.x = N₁*x₁+N₂*x₂
        xg.y = N₁*y₁+N₂*y₂
        xg.𝑤 = xg.w/2*𝐿
        push!(𝓖,xg)
    end

    nₑ = length(elms["Γᵍ"])
    elements["Γ̄ᵍ"] = TRElement{:Seg2}[]
    𝓒 = Node{(:𝐼,:𝐽),2}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    data𝓒 = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    n₁ = zeros(nₑ)
    n₂ = zeros(nₑ)
    ∂𝝭∂x = zeros(2*nₑ)
    ∂𝝭∂y = zeros(2*nₑ)
    data𝓖 = Dict([
        :w=>(1,scheme[:w]),
        :ξ=>(1,scheme[:ξ]),
        :𝝭=>(1,scheme[:ξ]),
        :∂𝝭∂x=>(4,∂𝝭∂x),
        :∂𝝭∂y=>(4,∂𝝭∂y),
        :n₁=>(3,n₁),
        :n₂=>(3,n₂),
        :x=>(2,zeros(2*nₑ)),
        :y=>(2,zeros(2*nₑ)),
        :z=>(2,zeros(2*nₑ)),
        :𝑤=>(2,zeros(2*nₑ))
    ])
    G = 0
    for (C,a) in enumerate(elms["Γᵍ"])
        push!(elements["Γ̄ᵍ"],TRElement{:Seg2}((C-1,1,𝓒),(G,2,𝓖)))
        x₁ = a.x[a.i[1]]
        x₂ = a.x[a.i[2]]
        y₁ = a.y[a.i[1]]
        y₂ = a.y[a.i[2]]
        𝐿 = ((x₁-x₂)^2+(y₁-y₂)^2)^0.5
        n₁[C] = (y₂-y₁)/𝐿
        n₂[C] = (x₁-x₂)/𝐿
        I = findfirst(x->x==Set(a.i), 𝓑)
        push!(𝓒,Node{(:𝐼,:𝐽),2}((I,0),data𝓒))
        G += 1
        xg = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((1,G,C,G-1),data𝓖)
        ξ = xg.ξ
        ∂𝝭∂x[G] = n₁[C]*ξ
        ∂𝝭∂y[G] = n₂[C]*ξ
        N₁ = (1-ξ)/2
        N₂ = (1+ξ)/2
        xg.x = N₁*x₁+N₂*x₂
        xg.y = N₁*y₁+N₂*y₂
        xg.𝑤 = xg.w/2*𝐿
        push!(𝓖,xg)
        G += 1
        xg = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((2,G,C,G-1),data𝓖)
        ξ = xg.ξ
        ∂𝝭∂x[G] = n₁[C]*ξ
        ∂𝝭∂y[G] = n₂[C]*ξ
        N₁ = (1-ξ)/2
        N₂ = (1+ξ)/2
        xg.x = N₁*x₁+N₂*x₂
        xg.y = N₁*y₁+N₂*y₂
        xg.𝑤 = xg.w/2*𝐿
        push!(𝓖,xg)
    end
    return elements, 𝓑
end