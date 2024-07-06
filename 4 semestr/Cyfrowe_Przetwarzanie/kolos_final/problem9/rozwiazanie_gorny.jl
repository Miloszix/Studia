function rozwiazanie_9(;
    order::Int = 86,
    fp::Float64 = 109.0,
    f0::Float64 = 30.52,
    z::Vector{Int} = [14, 52, 72, 23],
)
    #-0.006889899137402053
    #missing

    M=order/2
    m=range(-M,M)
    h = -2*(f0/fp)*sinc.(2*(f0/fp).*m)
    #hamming = 0.54 .+ 0.46*cos.(2*pi.*m/(2*M+1))

    hann = 0.5 .+0.5*cos.(2*pi.*m/(2*M+1)) #robimy okno
    h.*= hann #mnozymy okno razy pasmo

    return sum(h[z])
end

rozwiazanie_9()

#hanning = 0.5 .+ 0.5*cos.(2*pi.*m/(2*M+1))
#blackman = 0.42 .+ 0.5*cos.(2*pi.*m/(2*M+1)) + 0.08*cos.(4*pi.*m/(2*M+1))
#hamming = 0.54 .+ 0.46*cos.(2*pi.*m/(2*M+1))
#triang = 1 .- abs.(m)/(M+1)