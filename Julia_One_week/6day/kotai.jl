using Plots

function calc_energy(M,μ,ε)
    E = 0
    ks = range(-pi,pi,length = M)
    filling = 0
    for kx in ks
        εkx = ε(kx) - μ
        if εkx < 0
          E += εkx  
          filling += 1
        end
    end
    return E/M, filling/M
end

function test2()
    μ = 0
    Ms = [10,50,100,500,1000,5000,10000]
    Es = []
    ε(kx) = -2cos(kx)
    for M in Ms
        @time E,filling = calc_energy(M,μ,ε)
        println("$M $E $filling")
        push!(Es,E)
    end
    plot(Ms,Es,xscale=:log10)
    savefig("TB_edep.png")
end
test2()

function test3()
    ε(kx) = -2cos(kx)
    M = 1000
    μs = range(-3,3,length=100)
    fillings=[]
    for μ in μs 
        E,filling = calc_energy(M,μ,ε)
        push!(fillings,filling)
    end
    plot(μs,fillings,label="filling")
    savefig("TB_mudep.png")
end

test3()

function bisection_method(xmin,xmax,f,eps;itamax = 20)
    fmin = f(xmin)
    fmax = f(xmax)
    @assert fmin*fmax < 0 "f(min)*f(max) should be less than 0!"
    for i in 1:itamax
        xmid = (xmin + xmax)/2
        fmid = f(xmid)
        if abs(fmid) < eps
            return xmid,fmid
        end
        if fmid < 0
            xmin = xmid
        else
            xmax = xmid
        end 
        println("$i $xmid $fmid")
    end
end

function test4()
    filling = 0.25
    M = 1000
    ε(kx)=-2cos(kx)
    f(μ)=calc_energy(M,μ,ε)[2] - filling
    eps = 1e-10
    μ_ans,err = bisection_method(-3,3,f,eps)
    println("μ = $μ_ans, $(calc_energy(M,μ_ans,ε)[2])")
end

test4()

function make_H(N,μ,V)
    H = zeros(Float64,N,N)
    for i in 1:N
        j = i + 1
        j += ifelse(j > N, -N, 0)
        H[i,j] = -1

        j = i-1
        j += ifelse(j < 1, N ,0)

        j=i
        H[i,j] = - μ + V(i)
    end
    return H
end

function make_H(Lx,Ly,μ,V)
    N = Lx*Ly
    H = zeros(Float64,N,N)
    ds = [(1,0),(-1,0),(0,1),(0,-1)]
    for ix in 1:Lx
        for iy in 1:Ly
            i = (iy-1)*Lx+ix
            for d in ds
                jx = ix + d[1]
                jx += ifelse(jx > Lx,-Lx,0)
                jx += ifelse(jx < 1,Lx,0)

                jy = iy + d[2]
                jy += ifelse(jy > Ly,-Ly,0)
                jy += ifelse(jy < 1,Ly,0)

                j = (jy-1)*Lx+jx
                H[i,j]+=-1
            end
            H[i,i] += -μ + V(ix,iy)
        end
    end
    return H
end

function calc_ldos(E, i, ene, ψ, η)
    ldos = 0.0
    for n in (1:length(ene))
        ldos += abs(ψ[i,n])^2*η/((E-ene[n])^2 + η^2)
    end
    return ldos
end

using Plots
using FFTW
using LinearAlgebra
function ldos_plot()
    Lx = 91
    Ly = 91
    μ = -0.2
    V0 = 1
    ix1 = 22
    iy1 = 38
    function V(ix,iy)
        v = ifelse(ix == ix1 && iy == iy1, V0,0)
        return v
    end
    H = make_H(Lx,Ly,μ,V)
    @time ene,ψ = eigen(H)
    ldos = zeros(Float64,Lx,Ly)
    η = 0.01
    E = 0
    for ix in 1:Lx
       for iy in 1:Ly
        i = (iy -1)*Lx+ix
        ldos[ix,iy] = calc_ldos(E,i,ene,ψ,η)
       end 
    end
    heatmap(1:Lx,1:Ly,ldos[:,:],aspect_ratio=:equal,xlims=(1,Lx),ylims=(1,Ly))
    savefig("TB_ldosplot_$μ.png")
    
    ldosfft = fft(ldos)
    ldosfft[1,1]=0
    ldosfftshift = fftshift(ldosfft)
    freq = fftshift(fftfreq(Lx,2pi))
    heatmap(freq,freq,abs.(ldosfftshift),aspect_ratio=:equal,xlims=(-pi,pi),ylims=(-pi,pi))
    savefig("TB_ldosfft_$μ.png")
    
end

ldos_plot()

