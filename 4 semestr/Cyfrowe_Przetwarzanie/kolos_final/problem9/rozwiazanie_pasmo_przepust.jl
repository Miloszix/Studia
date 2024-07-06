function rozwiazanie_9(;
    order::Int = 74,
    fp::Float64 = 123.0,
    f1::Float64 = 23.37,
    f2::Float64 = 60.27,
    z::Vector{Int} = [12, 17, 3, 32, 11, 68],
)
    #-0.06435056717791442
    #missing

    M=order/2
    m=range(-M,M)
    h=2*(f2/fp)*sinc.(2*(f2/fp).*m) - 2*(f1/fp)*sinc.(2*(f1/fp).*m)
    #hanning = 0.5 .+ 0.5*cos.(2*pi.*m/(2*M+1))
    hamming = 0.54 .+ 0.46*cos.(2*pi.*m/(2*M+1))

    h.*=hamming
    return sum(h[z])
end

rozwiazanie_9()

#hanning = 0.5 .+ 0.5*cos.(2*pi.*m/(2*M+1))
#blackman = 0.42 .+ 0.5*cos.(2*pi.*m/(2*M+1)) + 0.08*cos.(4*pi.*m/(2*M+1))
#hamming = 0.54 .+ 0.46*cos.(2*pi.*m/(2*M+1))
#triang = 1 .- abs.(m)/(M+1)


    