function rozwiazanie_7(;
    zz::Vector{ComplexF64} = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im],
    pp::Vector{ComplexF64} = ComplexF64[0.6701939367751139 - 0.6489278834033337im, 0.4615448869440356 + 0.7927714869863219im, 0.6701939367751139 + 0.6489278834033337im, 0.4615448869440356 - 0.7927714869863219im, 0.5357004078881282 - 0.6629989449417352im, 0.5357004078881282 + 0.6629989449417352im],
    k::Float64 = 0.00289819463372143,
    F::Vector{Float64} = [0.12, 0.24, 0.41, 0.44, 0.48],
)
    #0.14602268573773702
    #missing

    H=zeros(ComplexF64, length(F))

    omega = exp.(2*pi*im*F)
    for i in eachindex(H)
        B=1
        for j in eachindex(zz)
            B*=omega[i]-zz[j]
        end

        A=1
        for j in eachindex(pp)
            A*=omega[i]-pp[j]
        end

        H[i]=k*B/A
    end

    return sum(abs.(H)/length(H))
end

rozwiazanie_7()