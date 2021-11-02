println("Check the Growth rates")

using DelimitedFiles
using PyPlot

f1 = "glnorm1"
d1 = readdlm(f1,Float64)

# f2 = "glnorm2"
# d2 = readdlm(f2,Float64)
# 
# f3 = "glnorm3"
# d3 = readdlm(f3,Float64)

t1 = d1[:,1]
y1 = d1[:,2]

# t2 = d2[:,1]
# y2 = d2[:,2]
# 
# t3 = d3[:,1]
# y3 = d3[:,2]

semilogy(t1,y1)
# semilogy(t2,y2)
# semilogy(t3,y3)




