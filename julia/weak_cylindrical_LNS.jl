#!/usr/bin/julia

#    Derivation of Cylindrical coordinates using symbolic Julia

println("Symbolic Math in Julia using SymPy (Python Package)")
println("Derivation for Cylindrical Coordinates")

using PyCall
using SymPy

include("ReplaceUnitVecDeriv.jl")
include("ReplaceUnitVecDot.jl")
include("GetComponents.jl")
include("GradR.jl")
include("GradTheta.jl")
include("GradT.jl")
include("GradVector.jl")
include("GradScalar.jl")
include("Div.jl")
include("WeakDiv.jl")
include("Weak2DDiv.jl")
include("Dot.jl")

Sym = pyimport("sympy")


IfFourier = true

# # Cartesian Coordinates
# x = Sym.Symbol("x")
# y = Sym.Symbol("y")
# z = Sym.Symbol("z")

# Cylindrical Coordinates
t = Sym.Symbol("t")          # x in the cylindrical coordinates
r = Sym.Symbol("r", positive=true)
θ = Sym.Symbol("θ")

#x = t
#y = r*cos(θ)
#z = r*sin(θ)
#
#
# # x
# dxdr = Sym.diff(x,r)
# dxdθ = Sym.diff(x,θ)
# dxdt = Sym.diff(x,t)
# 
# # y
# dydr = Sym.diff(y,r)
# dydθ = Sym.diff(y,θ)
# dydt = Sym.diff(y,t)
# 
# # z
# dzdr = Sym.diff(z,r)
# dzdθ = Sym.diff(z,θ)
# dzdt = Sym.diff(z,t)
# 
# row_t = [dxdt  dydt dzdt]
# row_r = [dxdr  dydr dzdr]
# row_θ = [dxdθ  dydθ dzdθ]

# A     = Sym.simplify([row_t; row_r; row_θ])

# A_inv = Sym.simplify(Sym.inv_quick(A))
#
# dtdx = A_inv[0,0] 
# dtdy = A_inv[0,1]  
# dtdz = A_inv[0,2]
# drdx = A_inv[1,0] 
# drdy = A_inv[1,1]  
# drdz = A_inv[1,2]
# dθdx = A_inv[2,0] 
# dθdy = A_inv[2,1]  
# dθdz = A_inv[2,2]

# Unit Vectors in Clindrical coordinates
er = Sym.Function("er", commutative=false)(r,θ)
eθ = Sym.Function("eθ", commutative=false)(r,θ)
et = Sym.Function("et", commutative=false)()

# # Commutatitive symbols for collection
# erc = Sym.Function("erc", commutative=true)(r,θ)
# eθc = Sym.Function("eθc", commutative=true)(r,θ)
# etc = Sym.Function("etc", commutative=true)()

# # Unit Vectors in Cartesian coordinates
# ex = Sym.Symbol("ex")
# ey = Sym.Symbol("ey")
# ez = Sym.Symbol("ez")


# er1 = cos(θ)*ex + sin(θ)*ey
# eθ1 = -sin(θ)*ex + cos(θ)*ey
# et1 = ez

# derdθ = Sym.diff(er,θ)
# deθdθ = Sym.diff(eθ,θ)  
#replace(derdθ,-ex*sin(θ) + ey*cos(θ),eθ)
#replace(deθdθ,-ex*cos(θ) - ey*sin(θ),-er)

# Fourier Coefficients in Clindrical coordinates

if IfFourier
  ur = Sym.Function("ur\'",commutative=false)(r,t)
  uθ = Sym.Function("uθ\'",commutative=false)(r,t)
  ut = Sym.Function("ut\'",commutative=false)(r,t)
  p  = Sym.Function("p\'" ,commutative=false)(r,t)
else
  ur = Sym.Function("ur\'",commutative=false)(r,θ,t)
  uθ = Sym.Function("uθ\'",commutative=false)(r,θ,t)
  ut = Sym.Function("ut\'",commutative=false)(r,θ,t)
  p  = Sym.Function("p\'" ,commutative=false)(r,θ,t)
end  

# Base Flow variables
Ur = Sym.Function("Ur",commutative=false)(r,t)
Uθ = Sym.Function("Uθ",commutative=false)(r,t)
Ut = Sym.Function("Ut",commutative=false)(r,t)
P  = Sym.Function("P" ,commutative=false)(r,t)


# Density
ρ  = Sym.Symbol("ρ",commutative=true)
# Dynamic Viscosity
μ  = Sym.Symbol("μ",commutative=true)
# Kinematic Viscosity
ν  = μ/ρ


# Wavenumber
k  = Sym.Symbol("k",commutative=true)
j  = sqrt(Complex(-1.0))                        # well you know...

expz  = exp(j*k*θ)
expzi = exp(-j*k*θ)

# Test Function
if IfFourier
  vrt      = Sym.Function("v",commutative=false)(r,t)
  v        = vrt*expzi
else
  v        = Sym.Function("v",commutative=false)(r,θ,t)
end

if IfFourier
  u_right  = ur*expz*er + uθ*expz*eθ + ut*expz*et
  u_left   = er*ur*expz + eθ*uθ*expz + et*ut*expz
  p        = p*expz
else
  u_right  = ur*er + uθ*eθ + ut*et
  u_left   = er*ur + eθ*uθ + et*ut
end

U_right  = Ur*er + Uθ*eθ + Ut*et
U_left   = er*Ur + eθ*Uθ + et*Ut


gradu,mGu  = GradVector(u_right,r,θ,t,"left")
mGuT       = mGu.transpose()

# 
GjR   = er*mGu[1,1] + eθ*mGu[2,1] + et*mGu[3,1]
#GRj   = er*mGu[1,1] + eθ*mGu[1,2] + et*mGu[1,3]

Gjθ   = er*mGu[1,2] + eθ*mGu[2,2] + et*mGu[3,2]
#Gθj   = er*mGu[2,1] + eθ*mGu[2,2] + et*mGu[2,3]

Gjt   = er*mGu[1,3] + eθ*mGu[2,3] + et*mGu[3,3]
#Gtj   = er*mGu[3,1] + eθ*mGu[3,2] + et*mGu[3,3]

# Gradient Transpose
GTjR   = er*mGuT[1,1] + eθ*mGuT[2,1] + et*mGuT[3,1]
#GTRj   = er*mGuT[1,1] + eθ*mGuT[1,2] + et*mGuT[1,3]

GTjθ   = er*mGuT[1,2] + eθ*mGuT[2,2] + et*mGuT[3,2]
#GTθj   = er*mGuT[2,1] + eθ*mGuT[2,2] + et*mGuT[2,3]

GTjt   = er*mGuT[1,3] + eθ*mGuT[2,3] + et*mGuT[3,3]
#GTtj   = er*mGuT[3,1] + eθ*mGuT[3,2] + et*mGuT[3,3]

graduT = GTjR*er + GTjθ*eθ + GTjt*et

# Strain-rate
ϵ      = 1/2*(gradu + graduT)
# Stress-rate
σ      = 2*ν*ϵ

#ERj   = (er*mGU[1,1] + er*mGU[1,1]) + (eθ*mGU[2,1] + eθ*mGU[1,2]) + (et*mGU[3,1] + et*mGU[1,3])
ER    = GjR + GTjR
Eθ    = Gjθ + GTjθ
Et    = Gjt + GTjt

#WkdivER = Weak2DDiv(ER*er,v,"left")
#WkdivEθ = Weak2DDiv(Eθ*eθ,v,"left")
#WkdivEt = Weak2DDiv(Et*et,v,"left")

#WkdivER = WeakDiv(ER*er,v,"left")
#WkdivEθ = WeakDiv(Eθ*eθ,v,"left")
#WkdivEt = WeakDiv(Et*et,v,"left")

#Wkdivall = WkdivER + WkdivEθ + WkdivEt
#

#divER    = Div(ER*er,"left")
#divEθ    = Div(Eθ*eθ,"left")
#divEt    = Div(Et*et,"left")

#divall   = divER + divEθ + divEt

divall   = Div(σ,"left")

divcom   = GetComponents(divall)

gradUb,mGUb  = GradVector(U_right,r,θ,t,"left")

conv1    = Dot(U_right,gradu)
conv2    = Dot(u_right,gradUb)

conv     = conv1 + conv2

convcom  = GetComponents(conv)

gradP,mGP = GradScalar(p,r,θ,t,"right")

# This is in the "Stress" formulation
# Some Divergence terms need to be eliminated
# to get the simplified Laplacian term
NS       = convcom .+ (1. /ρ)*mGP .- divcom 


# Weak Formulation

gradV,mGV   = GradScalar(v,r,θ,t,"right")


conv1       = Dot(U_right,gradu)
conv2       = Dot(u_right,gradUb)

Wkconv      = v*(conv1 + conv2)
Wkconvcom   = GetComponents(Wkconv)

WkgradP     = -p*gradV
WkgradPcom  = -p.*mGV

WkLapl      = Dot(-gradV,σ)
WkLaplcom   = GetComponents(WkLapl)

WkNS        = Wkconvcom .+ (1. /ρ)*WkgradPcom .- WkLaplcom 

println("Done")














