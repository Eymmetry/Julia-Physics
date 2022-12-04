function Basel(nc)
    a = 0
    for n in 1:nc
        a+=1/n^2
    end
    return sqrt(6a)
end

function test()
    ncs = [1,10,100,1000,10000]
    for nc in ncs
        b = Basel(nc)
        println("バーゼル級数和(nc = $nc): ",b,"\t",abs(π-b)/π)
    end
end

test()

function test()
    ncs = [10^n for n = 0:9]
    bs = []
    for nc in ncs
        b = Basel(nc)
        push!(bs, abs(pi-b)/pi)
    end
    println(bs)
    return ncs,bs
end

ncs,bs = test()

using Plots

plot(ncs,bs)
savefig("basel.pdf")
plot(ncs,bs,xscale=:log10, yscale=:log10,markershape = :circle,label="Basel", xlabel="cutoff num", ylabel = "relative error")

savefig("basel.pdf")


function Leibniz(nc)
    a = 0
    for n in 0:nc
        a+=(-1)^n/(2n+1)
    end
    return 4a
end


function test()
    ncs = [10^n for n = 0:9]
    bs = []
    ls = []
    for nc in ncs
        b = Basel(nc)
        l = Leibniz(nc)
        push!(bs, abs(pi-b)/pi)
        push!(ls, abs(pi-l)/pi)
    end
    println(bs)
    return ncs,bs,ls
end

ncs,bs,ls = test()

plot(ncs,bs,xscale=:log10, yscale=:log10,markershape = :circle,label="Basel", xlabel="cutoff num", ylabel = "relative error")
plot!(ncs,ls,xscale=:log10, yscale=:log10,markershape = :star5,label="Leibniz", xlabel="cutoff num", ylabel = "relative error")

function Ramanujan(nc)
    a = 0
    for n in 0:nc
        n = big(n)
        a += factorial(4n)*(1103+26390n)/(4^n*99^n*factorial(n))^4
    end
    return 99^2/(2*sqrt(big(2))*a)
end

