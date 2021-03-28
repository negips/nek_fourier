function GradTheta(U,r,θ,t)
       
  rp1   = (Sym.diff(er,r),0)
  rp2   = (Sym.diff(er,θ),eθ)
  rp3   = (Sym.diff(eθ,r),0)
  rp4   = (Sym.diff(eθ,θ),-er)

  rp5   = (et*et,1.0)
  rp6   = (et*er,0.0)
  rp7   = (et*eθ,0.0)
  rp8   = (er*er,1.0)
  rp9   = (er*eθ,0.0)
  rp10  = (eθ*eθ,1.0)

  dudθ  = eθ*(1/r)*Sym.diff(U,θ)
  tmp   = gu.subs(Dict([rp1 rp2 rp3 rp4]))

# Gradients
#  tmp = gu.subs(Dict([rp1 rp2 rp3 rp4 rp5 rp6 rp7 rp8 rp9 rp10]))
  dudθ = tmp

  return dudr
end 
