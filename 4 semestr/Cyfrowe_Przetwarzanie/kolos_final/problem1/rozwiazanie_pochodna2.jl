function rozwiazanie_1(;
    fp::Float64 = 404.29,
    t1::Float64 = -1.31,
    N::Int = 150,
)  #0.04930332261521972
    #missing

    dt=1/fp
    
    g=t->-2*rem(t,1,RoundNearest)
    t=t1:dt:t1+(N-1)*dt
    y=2.7*g.(2.4.*t.-1.1)
    return sum(y.^2)/length(y)

end

rozwiazanie_1()