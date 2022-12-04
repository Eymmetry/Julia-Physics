function MC(N)
    cnt=0
    for n in 1:N
        x = rand()
        y = rand()
        r=x^2+y^2
        cnt+=ifelse(r>1,0,1)
    end
    return 4*cnt/N
end

N = 10000
println(MC(N))

using Random


function MC(N;seed = 11)
    Random.seed!(seed)
    cnt=0
    for n in 1:N
        x = rand()
        y = rand()
        r=x^2+y^2
        cnt+=ifelse(r>1,0,1)
    end
    return 4*cnt/N
end
println(MC(N))

using Plots
function test()
    ncs = [10^n for n=0:9]
    ms = []
    for nc in ncs
        m = MC(nc)
        push!(ms,abs(pi-m)/pi)
    end
    println(ms)
    plot(ncs, ms, xscale=:log10, yscale=:log10, markershape=:circle, label="Monte Celo",xlabel="num",ylabel="relative error")
    savefig("mc.png")
    return 
end

test()

function Ransu(N)
    xn = 0:0.01:1
    yn = sqrt.(1 .- xn.^2)
    plot(xn,yn,label = nothing, aspect_ratio=:equal,xlims=(0,1),ylims=(0,1))
    xin = []
    xout = []
    yin = []
    yout = []
    for n in 1:N
        x = rand()
        y = rand()
        r = x^2+y^2
        if r>1
            push!(xout,x)
            push!(yout,y)
        else
            push!(xin,x)
            push!(yin,y)
        end
    end
    plot!(xin,yin,label="inside",seriestype=:scatter,markercolor=:blue)
    plot!(xout,yout,label="outside",seriestype=:scatter,markercolor=:yellow)
    savefig("random.png")
    return length(xin)/N
end

Ransu(1000)