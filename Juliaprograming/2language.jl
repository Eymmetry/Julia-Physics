x=1
2x+1
2(x-3)^2 - 3(x-2)
θ = pi / 4
sin(θ)+cos(θ)

#2.1.2プリミティブ型
typeof(1)
typeof(0.5)
Sys.WORD_SIZE
Int
Inf64
Inf32
pi
VERSION

s = "Hello Julia"
s[1]
s[end]
s[1:5]
s = "こんにちは, Julia"
chars = Vector{Char}(s)
typeof(chars)

#2.2制御構文

#2.3関数
function add(x, y)
    return x + y
end

#2.4型
struct Point
    x :: Int
    y :: Int
end

function distance(p :: Point)
    sqrt(p.x^2 + p.y^2)
end

p = Point(2, 3)
distance(p)

mutable struct Point2{T, U}
    x :: T
    y :: U
end

function distance(p :: Point2{T, U}) where {T, U}
    sqrt(p.x^2 + p.y^2)
    @__MODULE__
end

#2.5コレクション
list = []
list = [1]
list = [1, 2]
push!(list, 3)
pop!(list)
list
insert!(list, 2, 4)
deleteat!(list, 1)
list2 = [1 2]

d = Dict{String, Int}()
d["a"] = 1; d["b"] = 2;
d

#2.6多次元配列
A = rand(3, 2)
eltype(A)
length(A)
ndims(A)
size(A)
size(A, 1)

A = collect(reshape(1:9, 3, 3))

#2.7モジュール
@__MODULE__

1:3 .|> (x -> x^2) |> sum |> sqrt

1:3 .|> x -> x^2 |> sum 


a = rand(2,1)
A = rand(2,3)

@time broadcast(+, a, A)
@time a .+ A

b = rand(1,2)

@time broadcast(+, a, b)
@time a .+ b