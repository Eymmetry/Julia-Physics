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

#2.4配列
a = [1, 2]
A = [1 2
     3 4]

b = A * a

a = (1/sqrt(2)) * [1, 1]
b = (1/sqrt(2)) * [1, -1]
A = [0 -im; im 0]
c = b' * A * a
c = (A' * b)' * a

B = zeros(3, 4)
B[1, 2] = 3
B
function make_matrix(n)
    H = zeros(n,n)
    for i in 1:n
        j = i + 1
        j += ifelse(j > n, -n, 0)
        H[i,j] = -1
        j = i - 1
        j += ifelse(j < 1, n, 0)
        H[i, j] = -1
        H[i, i] = 2
    end
    return H
end
make_matrix(10)

using LinearAlgebra
a = [3, 2im]
a ⋅ a

A = [1 -2
        1 1]

e, v = eigen(A)
typeof(e)
typeof(v)

A * v[:,1] - e[1] * v[:,1]

v ./= (Real(v[1,1]))
v

B = Matrix{ComplexF64}(undef, 2, 3)
C = Array{Float64}(undef,2,2,2)

diagm(-3=>[1,2])

a =  [1 2;3 4]
a[:]
a[1,1]

a
reshape(a,(4,1))
reshape(a,(1,4))

function y2(x, W1, W2, b1, b2)
    f(x) = tanh(x)
    return W2*f.(W1*x+b1) .+ b2
end
function y3(x, W1, W2, b1, b2)
    f(x) = tanh(x)
    return W2*f.(W1*x+b1) + b2 #これはエラー
end

W2 = [1 2]
W1 = [1 2 3;4 5 6]
x = [4, 5, 6]
b1 = [3, 5]
b2 =  8
y2(x,W1,W2,b1,b2)

P = ones(3,3)
b = 1
P .+ b

A = [im*pi 1; 1 im*pi]
(x->x^2)(A) #これは行列の2乗
(x->x^2).(A) #これは成分ごとの2乗 ,どちらも無名関数を用いている

#2.5型と多重ディスパッチ
typemin(Int64)
a = -9223372036854775808
typeof(a)
a-1

f(x, y) = x * y
function f(x :: Int64, y :: String)
    z = ""
    for i in 1:x
        z *= y
    end
    return z
end
methods(f)

*
methods(*)

function Base.:*(x :: Int64, y :: String)
    z = ""
    for i in 1:x
        z *= y
    end
    return z
end
*
3 * "dog"

function Base.:*(x :: Integer, y :: String)
    z = ""
    for i in 1:x
        z *= y
    end
    return z
end
*

supertype(Int32)
supertype(Signed)
supertype(Integer)
supertype(Real)
supertype(Number)
supertype(Any)

subtypes(Real)

function get_subtypes(type, num) #型のヒエラルキーを列挙する関数,(num * " "の演算が上のように定義されている必要がある)
    types = subtypes(type)
    num += 1
    if length(types) > 1
        for subtypes in types
            println(num * " ", subtypes)
            types = get_subtypes(subtypes, num)
        end
    end
    return types
end

get_subtypes(Integer, 1)
get_subtypes(Number, 1)

mutable struct Atom
    r
    v
    mass
end

m1 = 0.3
r1 = [0.2, 0.5, 0.1]
v1 = [0.3, 2, -1]

atom1 = Atom(r1, v1, m1)
atom1.mass

mutable struct Atom_new
    r::Array{Float64,1}
    v::Array{Float64,1}
    mass::Float64
end

display(atom1)

function Base.display(a :: Atom)
    println("r = ", a.r)
    println("v = ", a.v)
    println("mass = ", a.mass)
end

display(atom1)

