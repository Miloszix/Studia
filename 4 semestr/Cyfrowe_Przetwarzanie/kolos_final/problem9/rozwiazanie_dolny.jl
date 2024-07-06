function rozwiazanie_9(;
    order::Int = 60,
    fp::Float64 = 149.0,
    f0::Float64 = 67.05,
    z::Vector{Int} = [34, 32, 29],
)
    #0.08934907604478452
    #missing

    M=order/2

    m=range(-M,M)
    h=2*(f0/fp)*sinc.(2*(f0/fp).*m) #robimy pasmo
    hann = 0.5 .+0.5*cos.(2*pi.*m/(2*M+1)) #robimy okno
    h.*= hann #mnozymy okno razy pasmo

    return sum(h[z])
end

rozwiazanie_9()

#hanning = 0.5 .+ 0.5*cos.(2*pi.*m/(2*M+1))
#blackman = 0.42 .+ 0.5*cos.(2*pi.*m/(2*M+1)) + 0.08*cos.(4*pi.*m/(2*M+1))
#hamming = 0.54 .+ 0.46*cos.(2*pi.*m/(2*M+1))
#triang = 1 .- abs.(m)/(M+1)