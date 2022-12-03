include("original.jl")
using .Original
function tasu(a,b)
    println("tasu in main.jl")
    return a + 2b
end
println(Original.supercoolfunction(4, 10), " in main.jl")