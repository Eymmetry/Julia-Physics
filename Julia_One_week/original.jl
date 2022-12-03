module Original
    function tasu(a, b)
        println("tasu in original.jl")
        return a + b
    end
    function supercoolfunction(a, b)
        return tasu(a, b)/b
    end
end