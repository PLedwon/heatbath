using QuadGK,FFTW, InverseLaplace

gamma=1.2::Float64
M=1.0::Float64
tspan=(0.0,10^(1.0))::Tuple{Float64,Float64}
dt=0.1::Float64
times=Array{Float64}(tspan[1]:dt:tspan[2])
R=Array{Float64}(zeros(size(times)))

function K(t)
    cos(gamma*atan(t))*(t^2+1)^(-gamma/2)
end
    
function KLap(times)
    fftw(K.(-im*times))
end
    
function H(times)
    Talbot(s -> 1/(s*(M*s+KLap(times))) ,100)
end

function computeR!(R,tspan,dt)
    for i in eachindex(times)
        R[i]= quadgk(t -> H(t), times[1], tspan[i], rtol=1e-5)
        end
end

computeR!(R,tspan,dt)
print(R)
