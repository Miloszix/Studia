using Polynomials

function rozwiazanie_8(;
    b::Vector{Float64} = [4.953522354214143e-7, 0.0, -2.9721134125284857e-6, 0.0, 7.430283531321214e-6, 0.0, -9.907044708428285e-6, 0.0, 7.430283531321214e-6, 0.0, -2.9721134125284857e-6, 0.0, 4.953522354214143e-7],
    a::Vector{Float64} = [1.0, -9.021258384695642, 39.886213550216425, -112.94469483938911, 226.88388601435287, -339.5114511868411, 387.42822073864284, -339.5114511868411, 226.8838860143528, -112.94469483938909, 39.88621355021641, -9.021258384695642, 0.9999999999999996],
)
z=Polynomial(reverse(a))
@show typeof(z)
    poles = roots(Polynomial(reverse(a))) # tutaj potrzebny moduł Polynomials, nie wiem czy można go użyć na kolosie

    if any(trunc(abs(p), digits = 6) > 1 for p in poles)
        return -1.0
    elseif any(trunc(abs(p), digits = 6) == 1 for p in poles)
        return 0.0
    else
        return 1.0
    end
    # 0.0
end

rozwiazanie_8()