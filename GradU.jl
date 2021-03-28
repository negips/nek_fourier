function GradU(U,r,θ,t)
       
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

  gu    = Sym.expand(er*Sym.diff(U,r) + eθ*Sym.diff(U,θ)*(1. /r) + et*Sym.diff(U,t))
  tmp   = gu.xreplace(Dict([rp1 rp2 rp3 rp4]))

# Gradients
#  tmp = gu.subs(Dict([rp1 rp2 rp3 rp4 rp5 rp6 rp7 rp8 rp9 rp10]))
  gradu = tmp

#  etet  = gradu.collect(et*et)
#  eθet  = gradu.collect(eθ*et)
#  eret  = gradu.collect(er*et)
#
#  eteθ  = gradu.collect(et*eθ)
#  eθeθ  = gradu.collect(eθ*eθ)
#  ereθ  = gradu.collect(er*eθ)
#
#  eter  = gradu.collect(et*er)
#  eθer  = gradu.collect(eθ*er)
#  erer  = gradu.collect(er*er)
#
#  graduM = [etet eter eteθ;
#            eret erer ereθ;
#            eθet eθer eθeθ]

  return gradu
end 
