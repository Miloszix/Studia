function rozwiazanie_1(;
    fp::Float64 = 403.47,
    t1::Float64 = -2.47,
    N::Int = 933,
)
    #0.04930332261521972
    #missing

    dt=1/fp
    g=t->ifelse(mod(t,1)<0.5 , 1, -1)
    t=t1:dt:t1+(N-1)*dt
    y=2.7*g.(2.4.*t.-1.1)
    return sum(y.^2)/length(y)

end

rozwiazanie_1()