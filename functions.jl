#masses
function computeMasses(oscMass,ω,ωMin,γ,Ω)
    oscMass*(ω/ωMin).^(γ-3).*exp.(-ω/Ω)
end

function computeSpringConstants(m,ω)
    m.*ω.^2
end

function checkEnergy(u,v)
    v.^2+u
end





#k
