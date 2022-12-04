using DifferentialEquations

f(u,p,t) = 1.01u

u0 = 1/2
tspan = (0.0, 1.0)
prob = ODEProblem(f,u0,tspan)

sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

nt = 50

t = range(0.0,1.0,length=nt)
for i in 1:nt
    println("t = $(t[i]),solution: $(sol(t[i])), exact solution $(0.5*exp(1.01t[i]))")
end