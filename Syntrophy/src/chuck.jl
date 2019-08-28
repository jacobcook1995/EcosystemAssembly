# Throwaway script to test that my PyPlot can annotate correctly
using PyPlot

x = collect(0.0:0.1:10.0)
y = collect(0.0:0.1:10.0)
plot(x,y)
annotate("A",xy=[-2.0;10.0],annotation_clip=false)
savefig("Output/test2.png")
