function rozwiazanie_8(;
    z::Vector{ComplexF64} = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im],
    p::Vector{ComplexF64} = ComplexF64[0.879940252482802 - 0.1703384931933833im, 0.879940252482802 + 0.1703384931933833im, 0.8011510705587512 - 0.0im],
    k::Float64 = 0.8022305603808205,
) 
    #0.0
    #semistable
    #missing

    l=abs.(p)
    b=1
    for i in eachindex(l)
       if l[i]>1
        b=-1
        break
       elseif l[i]==1
        b=0 
       end
    end
    return b
end

rozwiazanie_8()