function computeMasses(oscMass::Float64,ω::Array{Float64,1},ωMin::Float64,γ::Float64,Ω::Float64)
    oscMass*(ω/ωMin).^(γ-3).*exp.(-ω/Ω)
end

function computeSpringConstants(m::Array{Float64,1},ω::Array{Float64,1})
    m.*ω.^2
end

@fastmath function makeTimestep(QOld,qOld,POld,pOld,k,M,mInv, dt)           
        PNew = POld + (dot(kl,(qOld.-QOld)))* dt #-V'(Q)*dt  for potential
        pNew = pOld - kl .*  (qOld.-QOld)  * dt
        QNew = QOld + PNew./M              * dt 
        qNew = qOld + pNew.*mInvl          * dt 
    
    return PNew,pNew,QNew,qNew
end


@fastmath function solveEOM!(QNew,qNew,PNew,pNew,QOld,qOld,POld,pOld,k,M,mInv,dt,tspan,qSave,PSave,energyError,momentumError)
           
    kH=Array{Float64}(vcat(0.0,k))       
    mInvH=Array{Float64}(vcat(1/M,mInv))       
    initialEnergy=H(vcat(QOld,qOld),vcat(POld,pOld),kH,mInvH)::Float64
    initialMomentum=sum(pOld)+POld::Float64
           
    for i in eachindex(tspan[1]:dt:(tspan[2]+saveIndex*dt)) 
     P,p,Q,q = makeTimestep(QOld,qOld,POld,pOld,k,M,mInv,dt)
      
      if i % saveIndex == 0
           j= Int(i/saveIndex) 
           qSave[j] = Q
           pSave[j] = P           
           energyError[j]=computeE(vcat(Q,q),vcat(P,p),mInvH,kH,initialEnergy)
           momentumError[j]=computeM(vcat(P,p),initialMomentum)
           
       end 
        
       QOld=Q
       qOld=q
       POld=P
       pOld=p
        
    end 
    
end
    
@fastmath function interact(q::Array{Float64,1})
    (q.-q[1]).^2

end

@fastmath function HInt(q::Array{Float64,1}, k::Array{Float64,1})
        0.5*dot(k,interact(q))
end

@fastmath function HMom(p::Array{Float64,1},mInv::Array{Float64,1})
        0.5*dot(mInv,p.^2)
end   

@fastmath function H(q::Array{Float64,1},p::Array{Float64,1},mInv::Array{Float64,1},k::Array{Float64,1}) 
     HMom(p,mInv)+ HInt(q,k)     
end

function computeE(q::Array{Float64,1},p::Array{Float64,1},mInv::Array{Float64,1},k::Array{Float64,1})
    (H(q,p,mInv,k))#-initialEnergy)/initialEnergy
end

function computeM(p::Array{Float64,1})
        sum(p)
end





