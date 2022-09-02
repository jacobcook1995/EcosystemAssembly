using Plots
using TradeOff

function wrap_test()
    ωs, λl, ϕsl, χfl, χvl, λh, ϕsh, χfh, χvh = ω_test()
    # Extract max ribosome fraction
    plot(ωs,λl,label="")
    savefig("Output/Fig8/GrowthLow.png")
    plot(ωs,ϕsl,label="")
    savefig("Output/Fig8/PhiLow.png")
    plot(ωs,χfl,label="")
    savefig("Output/Fig8/FixedCostLow.png")
    plot(ωs,χvl,label="")
    savefig("Output/Fig8/VariableCostLow.png")
    plot(ωs,λh,label="")
    savefig("Output/Fig8/GrowthHigh.png")
    plot(ωs,ϕsh,label="")
    savefig("Output/Fig8/PhiHigh.png")
    plot(ωs,χfh,label="")
    savefig("Output/Fig8/FixedCostHigh.png")
    plot(ωs,χvh,label="")
    savefig("Output/Fig8/VariableCostHigh.png")
    return(nothing)
end

@time wrap_test()
