using Plots
using TradeOff

function wrap_test()
    ωs, λf, af, ϕsf = ω_test()
    # Extract max ribosome fraction
    plot(ωs,λf,label="")
    savefig("Output/Fig8/FixedGrowth.png")
    plot(ωs,af,label="")
    savefig("Output/Fig8/Fixeda.png")
    plot(ωs,ϕsf,label="")
    savefig("Output/Fig8/FixedPhi.png")
    return(nothing)
end

@time wrap_test()

# C, T = ω_test()
# plot(T,C[:,1],label="")
# savefig("Output/Fig8/TestPops.png")
# plot(T,C[:,2],label="")
# savefig("Output/Fig8/Testa.png")
# plot(T,C[:,3],label="")
# savefig("Output/Fig8/Testphi.png")
# plot(T,log.(C[:,1]),label="")
# savefig("Output/Fig8/TestLogPops.png")
