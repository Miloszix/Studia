function rozwiazanie_1(;
    fp::Float64 = 422.28,
    t1::Float64 = -3.78,
    N::Int = 449,
)

    #2.411809372446129
    #missing
    dt=1/fp

    t=t1:dt:t1+(N-1)*dt

    triangular_wave(t::Real)::Real = ifelse(mod(t + 1 / 4, 1.0) < 1 / 2, 4mod(t + 1 / 4, 1.0) - 1, -4mod(t + 1 / 4, 1.0) + 3)

    y=4.1*triangular_wave.(1.8.*t.-1.6)

    return sqrt(sum(y.^2)/length(y))
end

rozwiazanie_1()