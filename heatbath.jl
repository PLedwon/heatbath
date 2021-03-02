using DifferentialEquations, LinearAlgebra, ForwardDiff, Plots
include("functions.jl")

N=4000::Int64
tspan=(0.0,10.0)::Tuple{Float64,Float64}
beta=1.0::Float64 #1/(kB*T)
γ=1.2::Float64 #expected diffusion exponent
Ω = 1.0::Float64
M=1.0::Float64 #mass of distinguished particle
oscMass=100.0::Float64 #mass of heaviest bath oscillator
ωMin=N^(-0.8)::Float64 #eigenfrequency of slowest bath oscillator
ωMax=ωMin*N^(1.076)::Float64 #eigenfrequency of fastest bath oscillator
ω = Vector{Float64}(range(ωMin,stop=ωMax, length=N))
masses=Array{Float64}(computeMasses(oscMass,ω,ωMin,γ,Ω))
params=Array{Float64}(zeros((N+1,2)))
k=Array{Float64}(vcat(0.0,computeSpringConstants(masses,ω)))
m=Array{Float64}(vcat(M,masses))
mInv=Array{Float64}(1 ./m) #!!! mInv

const params[:,1]=mInv;
const params[:,2]=k;

# initial values for position and momentum
Q0=0.0::Float64
q0= Array{Float64}(vcat(Q0,beta^(-0.5)*k[2:end].^(-0.5) .* randn(N) .+Q0))
p0=Array{Float64}(beta^(-0.5)*m.^(0.5)  .* randn(N+1))
p0.-=1/(N+1)*sum(p0)
nSaved = 2000::Int64


ds=(tspan[2]-tspan[1])*1.0/nSaved::Int64;
initialEnergy=H(q0,p0, params)
initialMomentum=sum(p0);


#@fastmath function loopfunc(param::Float64,qi::Float64,q1::Float64)
#       -1.0*param*(qi-q1)
#    end


#@fastmath function qdot(dq,p,q, params)
#         dq = p.*params[:,1]
#         return dq
#    end

#@fastmath function pdot(dp::Array{Float64,1},p::Array{Float64,1},q::Array{Float64,1}, params::Array{Float64,2})
#    dp[1]= params[1,2]*(q[2]-q[1])
#    dp[2]=-dp[1]
#    @inbounds @simd for i in 1:length(q)-1
#            dp[i+1]=loopfunc(params[i,2],q[i+1],q[1])
#            dp[1]-=dp[i+1]
#       end
#        return dp
#    dp.= -k.*(q.-q[1])
#    dp[1] += sum(@view dp[2:end])
#    end
pdot(dp,p,q,params) = ForwardDiff.gradient!(dp, q->-H(q, p, params), q)
qdot(dq,p,q,params) = ForwardDiff.gradient!(dq, p-> H(q, p, params), p)

#dq = similar(q0)
#@time qdot(dq,p0,q0,params,tspan[1])


prob = DynamicalODEProblem(pdot, qdot, p0, q0, tspan, params);

@time begin
sol = solve(prob,SymplecticEuler(),dt=10^(-5), saveat=ds);
end


q=zeros(Float64,(N+1,nSaved));
p=zeros(Float64,(N+1,nSaved));


for i in 1:nSaved
    q[:,i] = sol.u[i][2,:]
    p[:,i] = sol.u[i][1,:]
end


energyError = zeros(nSaved);
momentumError = zeros(nSaved);
@simd for i in 1:nSaved
    energyError[i]= abs( (initialEnergy-H(q[:,i],p[:,i],params))/initialEnergy )
    momentumError[i]=abs( sum(p[:,i])-sum(p[:,1]) )
end

avgEnergyError=sum(energyError)/length(energyError)
maxEnergyError=maximum(energyError)
avgMomentumError=sum(momentumError)/length(momentumError)
maxMomentumError=maximum(momentumError)
println(avgEnergyError)
println(maxEnergyError)
println(avgMomentumError)
println(maxMomentumError)


plot_trajectory(sol) = plot(sol,vars=(1), lab="Trajectory", title="Distinguished particle")
plot_energyError(sol) = plot(energyError, lab="Energy variation", title="First Integrals")
plot_momentumError(sol) = plot(momentumError, lab="momentum variation", title="First Integrals")

plt_trajectory=plot(plot_trajectory(sol))
savefig(plt_trajectory,"Trajectory.pdf")

plt_energy=plot(plot_energyError(sol))
savefig(plt_energy,"energyError.pdf")

plt_momentum=plot(plot_momentumError(sol))
savefig(plt_momentum,"momentumError.pdf")
