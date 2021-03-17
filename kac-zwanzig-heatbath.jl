using LinearAlgebra, Plots, NPZ#, ForwardDiff
include("functionLib.jl")

N=4000::Int64
tspan=(0.0,10^(1.0))::Tuple{Float64,Float64}
dt = 1.0*10^(-5.0)::Float64
beta=1.0::Float64 #1/(kB*T)
γ=1.2::Float64 #expected diffusion exponent
Ω = 1.0::Float64
M=0.0001::Float64 #mass of distinguished particle
oscMass=10.0::Float64 #mass of heaviest bath oscillator
ωMin=N^(-0.7988)::Float64 #eigenfrequency of slowest bath oscillator
ωMax=ωMin*N^(1.0688)::Float64 #eigenfrequency of fastest bath oscillator
ω = Array{Float64}(range(ωMin,stop=ωMax, length=N))
masses=Array{Float64}(computeMasses(oscMass,ω,ωMin,γ,Ω))
k=Array{Float64}(vcat(0.0,computeSpringConstants(masses,ω)))
m=Array{Float64}(vcat(M,masses))
mInv=Array{Float64}(1 ./m) #!!!mInv is array with dim=N+1


Q0=0.0::Float64
q0= Array{Float64}(vcat(Q0,beta^(-0.5)*k[2:end].^(-0.5) .* randn(N) .+Q0))
p0=Array{Float64}(beta^(-0.5)*m.^(0.5)  .* randn(N+1))
p0.-=1/(length(p0))*sum(p0)
nSaved = 5000::Int64
dts=(tspan[2]-tspan[1])*1.0/nSaved::Int64
saveIndex = max(1,floor(Int,dts/dt)) ::Int64


kl=Array{Float64}(k[2:end])
mInvl=Array{Float64}(mInv[2:end])
QOld=q0[1]::Float64
POld=p0[1]::Float64
qOld=Array{Float64}(q0[2:end])
pOld=Array{Float64}(p0[2:end])
QNew=0.0::Float64
PNew=0.0::Float64
qNew=Array{Float64}(zeros(N))
pNew=Array{Float64}(zeros(N))
qSave=Array{Float64}(zeros(nSaved+1))
pSave=Array{Float64}(zeros(nSaved+1))

energyError = zeros(nSaved+1);
momentumError = zeros(nSaved+1);

@time begin
solveEOM!(QNew,qNew,PNew,pNew,QOld,qOld,POld,pOld,kl,M,mInvl,dt,tspan,qSave,pSave,energyError,momentumError)
end

Q=qSave
P=pSave


avgEnergyError =sum(energyError)/length(energyError)
maxEnergyError =maximum(energyError)
avgMomentumError =sum(momentumError)/length(momentumError)
maxMomentumError =maximum(momentumError)

name=rand(1)*10^5.0
name=string( "../npzFiles/",floor(Int,name[1]),".npz"   )

print(avgEnergyError, maxEnergyError)
print(avgMomentumError, maxMomentumError)
plt=plot(q[1,:])
savefig(plt,"traj.pdf")


################npzwrite(name, Dict("t"=>tspan[1]:dt:(tspan[2]+saveIndex*dt), "Q"=>q[1,:], "P"=>p[1,:], "energyError"=>energyError, "dts"=>dts ,"momentumError"=>momentumError, "avgEnergyError"=>avgEnergyError, "maxEnergyError"=>maxEnergyError,"avgMomentumError"=>avgMomentumError ,"maxMomentumError"=>maxMomentumError, "gamma"=> γ, "omegaMin"=>ωMin, "omegaMax"=>ωMax, "N"=>N ))
