function rozwiazanie_4(;
    fp::Int = 175,
    x::Vector{ComplexF64} = ComplexF64[0.31 + 0.03im, -0.38 - 0.15im, 1.04 - 0.04im, 1.19 + 0.71im, 0.11 + 1.51im, -0.06 - 0.75im, 0.07 - 0.92im, 0.09 - 0.03im, -0.51 + 0.04im, -0.34 - 0.04im, 0.32 - 0.53im, -0.56 + 0.98im, -1.41 - 0.31im, -0.19 + 0.39im, -0.03 - 0.64im, -1.17 + 0.41im, -0.17 + 0.33im, -0.55 - 0.48im, 0.11 + 1.09im, -0.83 + 0.63im, 0.42 + 0.47im, -0.29 + 1.21im, -1.65 - 0.09im, 0.99 + 0.66im, -0.33 - 1.13im],
    f::Vector{Int} = [56, 28, -49, 84],
)
    #1.0990464172144918
    #missing

    x_dft=dft(x)

    f_dft=f./(fp/length(x))

    result=zeros(Float64, length(f_dft))

    

    for i in eachindex(f_dft)
        index::Int64=f_dft[i]
        if index<0
            index=length(x)+index
        end
        #result[i]=angle(x_dft[index+1])
        result[i]=abs(x_dft[index+1])/length(x)

    end
return sum(result)
    
end

function dft(x)
    N=length(x)
    zeta=exp(-2*pi*im/N)

    [
        sum((x[n+1]*zeta^(n*k) 
        for n in 0:N-1)
        ) for k in 0:N-1

    ]
end

rozwiazanie_4()