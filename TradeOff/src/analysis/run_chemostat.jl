using Plots
using TradeOff

function wrap_test()
    C, T = ω_test()
    # Extract max ribosome fraction
    plot(T,C[:,1],label="")
    savefig("Output/Fig8/TestPops.png")
    plot(T,C[:,2:3],label="")
    savefig("Output/Fig8/TestConcs.png")
    plot(T,C[:,4],label="")
    savefig("Output/Fig8/Testa.png")
    plot(T,C[:,5],label="")
    savefig("Output/Fig8/Testphi.png")
    plot(T,log.(C[:,1]),label="")
    savefig("Output/Fig8/TestLogPops.png")
    return(nothing)
end

@time wrap_test()
