function rozwiazanie_9(;
    order::Int = 100,
    fp::Float64 = 153.0,
    f1::Float64 = 43.6,
    f2::Float64 = 74.97,
    z::Vector{Int} = [17, 87, 23, 8, 41],
)
    #-0.06435056717791442
    #missing

    M=order/2
    m=range(-M,M)
    h=-2*(f2/fp)*sinc.(2*(f2/fp).*m) + 2*(f1/fp)*sinc.(2*(f1/fp).*m)
    hanning = 0.5 .+ 0.5*cos.(2*pi.*m/(2*M+1))
    

    h.*=hanning
    return sum(h[z])
end

rozwiazanie_9()

#hanning = 0.5 .+ 0.5*cos.(2*pi.*m/(2*M+1))
#blackman = 0.42 .+ 0.5*cos.(2*pi.*m/(2*M+1)) + 0.08*cos.(4*pi.*m/(2*M+1))
#hamming = 0.54 .+ 0.46*cos.(2*pi.*m/(2*M+1))
#triang = 1 .- abs.(m)/(M+1)