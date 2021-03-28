#!/usr/bin/julia

#    Derivation of Cylindrical coordinates using symbolic Julia

println("Symbolic Math in Julia using SymPy (Python Package)")
println("Derivation for Cylindrical Coordinates")

using PyCall
using SymPy

include("GradU.jl")
include("DivU.jl")


Sym = pyimport("sympy")

# Cartesian Coordinates
x = Sym.Symbol("x")
y = Sym.Symbol("y")
z = Sym.Symbol("z")

# Cylindrical Coordinates
t = Sym.Symbol("t")          # x in the cylindrical coordinates
r = Sym.Symbol("r", positive=true)
θ = Sym.Symbol("θ")

x = t
y = r*cos(θ)
z = r*sin(θ)


# x
dxdr = Sym.diff(x,r)
dxdθ = Sym.diff(x,θ)
dxdt = Sym.diff(x,t)

# y
dydr = Sym.diff(y,r)
dydθ = Sym.diff(y,θ)
dydt = Sym.diff(y,t)

# z
dzdr = Sym.diff(z,r)
dzdθ = Sym.diff(z,θ)
dzdt = Sym.diff(z,t)

row_t = [dxdt  dydt dzdt]
row_r = [dxdr  dydr dzdr]
row_θ = [dxdθ  dydθ dzdθ]

A     = Sym.simplify([row_t; row_r; row_θ])

A_inv = Sym.simplify(Sym.inv_quick(A))
#
dtdx = A_inv[0,0] 
dtdy = A_inv[0,1]  
dtdz = A_inv[0,2]
drdx = A_inv[1,0] 
drdy = A_inv[1,1]  
drdz = A_inv[1,2]
dθdx = A_inv[2,0] 
dθdy = A_inv[2,1]  
dθdz = A_inv[2,2]

# Unit Vectors in Clindrical coordinates
er = Sym.Function("er", commutative=false)(r,θ)
eθ = Sym.Function("eθ", commutative=false)(r,θ)
et = Sym.Function("et", commutative=false)()

# Commutatitive symbols for collection
erc = Sym.Function("erc", commutative=true)(r,θ)
eθc = Sym.Function("eθc", commutative=true)(r,θ)
etc = Sym.Function("etc", commutative=true)()

# Unit Vectors in Cartesian coordinates
ex = Sym.Symbol("ex")
ey = Sym.Symbol("ey")
ez = Sym.Symbol("ez")


er1 = cos(θ)*ex + sin(θ)*ey
eθ1 = -sin(θ)*ex + cos(θ)*ey
et1 = ez

derdθ = Sym.diff(er,θ)
deθdθ = Sym.diff(eθ,θ)  
#replace(derdθ,-ex*sin(θ) + ey*cos(θ),eθ)
#replace(deθdθ,-ex*cos(θ) - ey*sin(θ),-er)

# Tangent Vectors in Clindrical coordinates
Ur = Sym.Function("Ur",commutative=false)(r,θ,t)
Uθ = Sym.Function("Uθ",commutative=false)(r,θ,t)
Ut = Sym.Function("Ut",commutative=false)(r,θ,t)

U_right  = Ur*er + Uθ*eθ + Ut*et
U_left   = er*Ur + eθ*Uθ + et*Ut

#gradU = [dUdr.__pyobject__[0] dUdr.__pyobject__[1] dUdr.__pyobject__[2]; 
#         RinvdUdθ.__pyobject__[0] RinvdUdθ.__pyobject__[1] RinvdUdθ.__pyobject__[2];
#         dUdt.__pyobject__[0] dUdt.__pyobject__[1] dUdt.__pyobject__[2]]

#gradU = Sym.simplify(gradU)

gradU  = GradU(U_right,r,θ,t)

#gradUM = Sym.simplify(gradUM)

#parsed_expr = Sym.parse_expr(gradU,evaluate=False)

Err   = Sym.diff(Ur,r)
Eθr   = r/2*Sym.diff(1/r*Uθ,r) + (1/2)*(1/r)*Sym.diff(Ur,θ)
Etr   = 1/2*Sym.diff(Ur,t) + 1/2*Sym.diff(Ut,r)

Erθ   = copy(Eθr)
Eθθ   = (1/r)*Sym.diff(Uθ,θ) + 1/r*Ur
Etθ   = (1/2)*(1/r)*Sym.diff(Ut,θ) + 1/2*Sym.diff(Uθ,t)

Ert   = copy(Etr)
Eθt   = copy(Etθ)
Ett   = Sym.diff(Ut,t)

Erj   = (er*Err + eθ*Eθr + et*Etr)*2
Eθj   = (er*Erθ + eθ*Eθθ + et*Etθ)*2
Etj   = (er*Ert + eθ*Eθt + et*Ett)*2

divEr = DivU(Erj,"left")
divEθ = DivU(Eθj,"left")
divEt = DivU(Etj,"left")
#
#
#divU  = Sym.simplify(DivU(U_left,"left"))

#divEr1  = divEr.subs(divU,0)

#divE2  = Sym.simplify(divE)

#tmp    = divE.subs(divU,0.)


VecdivEr = DivU(Erj*er,"left")
VecdivEθ = DivU(Eθj*eθ,"left")
VecdivEt = DivU(Etj*et,"left")

alldiv   = VecdivEr + VecdivEθ + VecdivEt


## Change to commutative Symbols
#rp1   = (er,erc)
#rp2   = (eθ,eθc)
#rp3   = (et,etc)
#
#tmp    = Sym.expand(alldiv.xreplace(Dict([rp1 rp2 rp3])))
#
#divR   = Sym.collect(tmp,erc)
#divθ   = Sym.collect(eθc)
#divt   = Sym.collect(etc)  

println("Done")








