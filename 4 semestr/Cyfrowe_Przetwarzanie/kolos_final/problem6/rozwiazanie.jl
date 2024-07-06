function rozwiazanie_6(;
    b::Vector{Float64} = [0.033194692256649484, -0.07172135788881488, 0.09489407529227799, -0.0717213578888149, 0.03319469225664949],
    a::Vector{Float64} = [1.0, 1.5921708442156168, 1.3135098485645686, 0.49616221058458576, 0.07954938181884069],
    x::Vector{Float64} = [-0.56, -0.74, -0.59, -0.37, 0.43, -0.96, -0.25, -0.96, 0.0, -0.38, 0.17, 0.74, 0.16, 0.07, 0.47, 0.46, -0.53, -0.95],
    L::Int = 57,
)
    #-0.0002646657612752812
    #missing
    N=length(x)
    M=length(b)
    K=length(a)

    y=zeros(Float64,L)

    for n in range(0,L-1)
        for m in range(0, M-1)
            if n-m>=0&&n-m<N
                y[n+1]+=b[m+1]*x[n-m+1]
            end

        end
        
        for k in range(1,K-1)
            if n-k>=0&&n-k<L
                y[n+1]-=a[k+1]*y[n-k+1]
            end
        end
    end

    return sum(y)/length(y)


end

rozwiazanie_6()