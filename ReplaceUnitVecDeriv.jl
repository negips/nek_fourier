function GradR(U,r,θ,t,dir)

# Calculate the result as a tensor
# if the input is a vector
#
  
  if dir == nothing
    dir = "left"
  end  

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

  if lower(dir)=="left"
    dudr  = er*Sym.diff(U,r)
  else
    dudr  = Sym.diff(U,r)*er
  end 

  tmp   = gu.subs(Dict([rp1 rp2 rp3 rp4]))

# Gradients
#  tmp = gu.subs(Dict([rp1 rp2 rp3 rp4 rp5 rp6 rp7 rp8 rp9 rp10]))
  dudr = tmp

  return dudr
end 
