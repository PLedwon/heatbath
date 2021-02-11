#using ForwardDiff
using Plots
using DifferentialEquations
include("functions.jl")

N=1
tspan=(0.0,10.0)
beta=1.0 #1/(kB*T)
γ=1.2 #expected diffusion exponent
Ω = 1.0
M=0.2::Float64 #mass of distinguished particle
oscMass=1.0::Float64 #mass of heaviest bath oscillator
ωMin=N^(-0.8323) #eigenfrequency of slowest bath oscillator
ωMax=ωMin*N^(1.05764) #eigenfrequency of fastest bath oscillator
ω = range(ωMin,stop=ωMax, length=N)
m=computeMasses(oscMass,ω,ωMin,γ,Ω)
k=computeSpringConstants(m,ω)
mInv=vcat(1/M, 1 ./m) #!!!mInv is array with dim=N+1

# initial values for position and momentum
Q0=0.0
P0=beta^(-0.5)*M^0.5*randn(1)
q0=beta^(-0.5)*k.^(-0.5) .* randn(N) .+Q0#::Float64
p0=beta^(-0.5)*m.^(0.5)  .* randn(N)#::Float64
u0=vcat(Q0, q0)
v0=vcat(P0, p0)
v0.-=1.0/(N+1)*sum(v0)
nSaved=5
ds=(tspan[2]-tspan[1])/nSaved

# equations of motion
function f1(dv,v,u,p,t) #f1 corresponds to -dH/dq=dp/dt
    dv[1]=0.0
    for i in N
        dv[1]+=k[i]*(u[i+1]-u[1])
    end
    dv[2:end]=-k.*(u[2:end].-u[1])
    return dv
end

function f2(du,v,u,p,t) #f2 corresponds to dH/dp=dq/dt
     du=mInv.*v
     return du
end

prob=DynamicalODEProblem(f1,f2,u0,v0,tspan)

@time begin
    sol=solve(prob,SymplecticEuler(),dt=0.001,saveat=ds)
end

checkEnergy(q0,p0)
println(sol.u)


plt=plot(sol,vars=(1))#,N+2))
savefig(plt,"phasespace.pdf")
