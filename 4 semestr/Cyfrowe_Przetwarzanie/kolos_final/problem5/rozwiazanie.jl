function rozwiazanie_5(;
    x::Vector{Float64} = [2.33, 2.24, -0.41, -4.5, -4.96, 1.57, 1.46, 1.73, 1.98, -1.39, -1.38, 4.9, -4.06, -3.88, 1.44, -3.99, -2.69, -4.39, 2.31, 3.83, 1.17, -2.52, -1.59, 2.38, 2.3, -2.61, -0.85, 1.14, 1.43, 2.68, 1.69, 0.32, -0.28, -1.71, 1.63, 2.69, -1.86, 3.28, -3.14, 0.26, -4.08, 1.46, 0.94, 3.47, 3.47, 0.25, -4.52, -2.68, -1.62, 3.26, 0.21, -3.92, 4.36, -3.67, 4.46, 1.43, 2.42, 1.69, 3.1, -4.09, -2.43, -3.84, -1.61, 4.77, -3.8, 2.76, 2.72, 3.97, -4.9, 2.95, -1.11],
    h::Vector{Float64} = [-2.31, -1.07, -4.13, -0.6, 1.47, -0.56, 2.73, -1.33, 3.91, -3.99, 3.34, -4.51, 4.82, 0.11, -0.19, 0.32, 0.85, -4.97, 1.5, 2.36],
)
    #114026.70335463
    #missing

    n=length(x)
    m=length(h)
    y=zeros(eltype(x), n+m-1)
    
    for i in 1:n
        for j in 1:m
            y[i+j-1]+=x[i]*h[j]
        end
    end

    return sum(y.^2)
end

rozwiazanie_5()