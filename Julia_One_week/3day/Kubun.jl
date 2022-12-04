function Kubun(N,x0,x1,f)
    dx=(x1-x0)/N
    a = 0.0
    xn = range(x0,step = dx, length=N)
    for x in xn
        a += f(x)
    end
    return a*dx
end
f(x)=4/(1+x^2)
N = 10000
p = Kubun(N,0,1,f)

function test()
    ncs = [10^n for n in 0:9]
    ks = abs.(pi .- Kubun.(ncs,0,1,f)) ./ pi
    return ncs,ks
end

ncs,ks = test()

using Plots
plot(ncs,ks,xscale=:log10,yscale=:log10,markershape=:circle,label="Kubun",xlabel="sep num",ylabel="relatibe error")


function Daikei(N,x0,x1,f)
    dx=(x1-x0)/N
    a = (f(x0)+f(x1))/2
    xn = range(x0,step = dx, length=N)
    for n in 2:N
        x = xn[n]
        a += f(x)
    end
    return a*dx
end

function test()
    ncs = [10^n for n in 0:9]
    ks = abs.(pi .- Kubun.(ncs,0,1,f)) ./ pi
    ds = abs.(pi .- Daikei.(ncs,0,1,f)) ./ pi
    return ncs,ks,ds
end

@time ncs,ks,ds = test()
plot(ncs,ks,xscale=:log10,yscale=:log10,markershape=:circle,label="Kubun",xlabel="sep num",ylabel="relatibe error")
plot!(ncs,ds,xscale=:log10,yscale=:log10,markershape=:star5,label="Daikei",xlabel="sep num",ylabel="relatibe error")0