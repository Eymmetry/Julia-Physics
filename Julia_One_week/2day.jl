f(x) = cos(x)+2sin(2x^2)
g(x) = exp(im*x)f(x)
g(3)
g(2im+2)
f(x,y) = cos(x) + 2sin(2y^2)

function f(θ)
    x = cos(θ)
    y = sin(θ)
    R = y / sqrt(x^2 + y^2)
    exp(R)

end

function f(x, a = 2)
    a * x
end
f(3)
f(0.1)
f(3, 4)

function g(x, y)
    x + y, x - y
end
a, b = g(2, 3)

a = g(2, 3)
typeof(a)
typeof(b)

a[1]

T(n, x) = n * acos(x) |> cos #cos(n * acos(x))
G(n, x) = n * acos(x) |> cos |> exp #exp(cos(n * acos(x)))

function ReLU(x)
    if x < 0
        zero(x)
    else
        x
    end
end

ReLU(x) = ifelse(x < 0, zero(x), x) #予め返り値を評価に注意

ix = 10
jx = ix + 1
jx += ifelse(jx > 10, -10, 0) #周期的境界条件とかで使う

x = -3
x > 0 ? sqrt(3) : x #三項演算子

xs = range(0, 2pi, length = 10)
for i in 1:10
    println("cos($(xs[i])) = ", cos(xs[i]))
end
collect(xs) #Vector(Array)にする

function g2(x, T; τ = 0.01, nmax = 1000000, eps = 1e-8)
    n = 0
    ωn = pi * T * (2n + 1)
    a = exp(im * ωn * τ) / (im * ωn -x)
    aold = a
    for n in 1:nmax
        ωn = pi * T * (2n + 1)
        a += exp(im * ωn * τ) / (im * ωn - x)
        ωn = pi * T * (2(-n) + 1)
        a += exp(im * ωn * τ) / (im * ωn - x)
        if abs(a - aold) / abs(aold) < eps
            println("converged at $n step")
            return T * a
        end
        aold = a
    end
    println("not converged in $nmax step")
    return T * a
end

function g2(x, T; τ = 0.01, eps = 1e-8)
    n = 0
    ωn = pi * T * (2n + 1)
    a = exp(im * ωn * τ) / (im * ωn -x)
    aold = 10 * a
    while abs(a - aold) / abs(a) > eps
        aold = a
        n += 1
        ωn = pi * T * (2n + 1)
        a += exp(im * ωn * τ) / (im * ωn - x)
        ωn = pi * T * (2(-n) + 1)
        a += exp(im * ωn * τ) / (im * ωn - x)
    end
    println("converged at $n step")
    return T * a
end

g2(0, 0.1; τ = 0.001, eps = 1e-15)