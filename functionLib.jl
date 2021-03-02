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


@fastmath function solveEOM!(QNew,qNew,PNew,pNew,QOld,qOld,POld,pOld,k,M,mInv,dt,tspan,qSave,PSave)
           
   @inbounds for i in eachindex(tspan[1]:dt:(tspan[2]+saveIndex*dt)) 
     P,p,Q,q = makeTimestep(QOld,qOld,POld,pOld,k,M,mInv,dt)
      
      if i % saveIndex == 0
           j= Int(i/saveIndex) 
           qSave[:,j] = vcat(Q,q)
           pSave[:,j] = vcat(P,p)          
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





