using Plots
using TradeOff

function wrap_test()
    ωs, λlf, λhf, λlv, λhv = ω_test()
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig8")
        mkdir("Output/Fig8")
    end
    # Extract max ribosome fraction
    plot(ωs, λlf, label = "")
    savefig("Output/Fig8/GrowthLowFixed.png")
    plot(ωs, λhf, label = "")
    savefig("Output/Fig8/GrowthHighFixed.png")
    plot(ωs, λlv, label = "")
    savefig("Output/Fig8/GrowthLowVar.png")
    plot(ωs, λhv, label = "")
    savefig("Output/Fig8/GrowthHighVar.png")
    return (nothing)
end

@time wrap_test()
