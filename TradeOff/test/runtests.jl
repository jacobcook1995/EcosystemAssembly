using Test
using TradeOff

function generate_test_reactions()
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(3, 3, 2e5, 293.15)
    # Preallocate vector of reactions
    test_reactions = Array{Reaction, 1}(undef, 3)
    for i in 1:3
        test_reactions[i] = make_Reaction(i, RP[i, 1], RP[i, 2], ΔG[i])
    end
    return test_reactions
end

function generate_test_microbe()
    # Unless otherwise noted these are the values used in the main simulation
    MC = 10^8
    γm = 1260.0 / 60.0
    n = zeros(Int64, 4)
    n[1] = 7459
    n[2:4] .= 300
    χl = 29.0
    Kγ = 5e8
    Pb = 0.7
    d = 6.0e-5
    ϕH = 0.45
    KΩ = 1e9
    fd = log(100) / log(2)
    kc = 10.0
    KS = (1 / 4) * 5.5e-3
    kr = 10.0
    # The value is randomised in the simulation, but going to fix it here
    ω = 0.5
    # Want a species with 2 reactions, relative expressions and η values differ, but kinetic
    # values don't
    R = 2
    η = [1.0, 2.0]
    ϕP = [0.7, 0.3]
    kcs = [kc, kc]
    KSs = [KS, KS]
    krs = [kr, kr]
    Reacs = [1, 2]
    # Make a microbe to help with testing
    test_microbe = make_Microbe(MC, γm, Kγ, χl, Pb, d, ϕH, KΩ, fd, ω, R, Reacs, η, kcs, KSs,
        krs, n, ϕP, 1, "test pool")
    return test_microbe
end

reactions = generate_test_reactions()
microbe = generate_test_microbe()

# Test that the Q function has sensible behaviour for positive and negative concentrations
@testset "Test Q function" begin
    @test Q(10.0, 2.0) == 1.0 / 5.0
    @test Q(2.0, 10.0) == 5.0
    @test Q(-10.0, 2.0) == -1.0 / 5.0
    @test Q(10.0, -2.0) == -1.0 / 5.0
    @test Q(-10.0, -2.0) == 1.0 / 5.0
end

# Want Keq functions gives something reasonable for positive and negative Gibbs free energy
# changes
@testset "Test Keq function" begin
    @test Keq(293.15, 1.0, -1e5) == 28478.23426211885
    @test Keq(293.15, 2.0, -1e5) == 1.2330306823576351e-9
    @test Keq(293.15, 1.0, 1e5) == 6.582768655458186e-32
end

# Test that θ function has sensible behaviour for positive and negative concentrations
@testset "Test θ function" begin
    @test θ(10.0, 2.0, 293.15, 1.0, -1e5) == 7.022907324912199e-6
    @test θ(10.0, 2.0, 293.15, 2.0, -1e5) == 1.62201965337624e8
    @test θ(0.0, 2.0, 293.15, 2.0, -1e5) == 1.0
    @test θ(-10.0, 2.0, 293.15, 2.0, -1e5) == 1.0
    @test θ(10.0, 0.0, 293.15, 2.0, -1e5) == 0.0
    @test θ(10.0, -2.0, 293.15, 2.0, -1e5) == 0.0
end

# Same tests for the smoothed version of the function
@testset "Test θ_smooth function" begin
    @test θ_smooth(10.0, 2.0, 293.15, 1.0, -1e5) == 7.022907324912199e-6
    @test θ_smooth(10.0, 2.0, 293.15, 2.0, -1e5) == 1.0
    @test θ_smooth(0.0, 2.0, 293.15, 2.0, -1e5) == 1.0
    @test θ_smooth(-10.0, 2.0, 293.15, 2.0, -1e5) == 1.0
    @test θ_smooth(10.0, 0.0, 293.15, 2.0, -1e5) == 0.0
    @test θ_smooth(10.0, -2.0, 293.15, 2.0, -1e5) == 0.0
end

# Test that the enzyme function works as expected
@testset "Test Eα function" begin
    @test Eα(0.35, microbe, 1) == 46666.666666666664
    @test Eα(0.35, microbe, 2) == 20000.0
    # invalid input generates bad values, but simulation setup stops them from being reached
    @test Eα(-0.35, microbe, 1) == 210000.00000000003
    @test Eα(1.05, microbe, 1) == -116666.66666666667
end

# Test that simplified q function has sensible behaviour for positive and negative
# concentrations
@testset "Test simplified q function" begin
    @test qs(microbe, 10.0, 2.0, 20000.0, 7.022907324912199e-6) == 199957.0585240782
    @test qs(microbe, 10.0, 2.0, 20000.0, 1.5e8) == 0.0
    @test qs(microbe, 10.0, 2.0, 20000.0, 1.0) == 0.0
    @test qs(microbe, 5.0, 2.0, 20000.0, 7.022907324912199e-6) == 199929.57391701653
    @test qs(microbe, 10.0, -2.0, 20000.0, 0.0) == 199972.50378073016
    @test qs(microbe, -10.0, 2.0, 20000.0, 1.0) == 0.0
    @test qs(microbe, -10.0, -2.0, 20000.0, 1.0) == 0.0
end

# Test that full version of q function has sensible behaviour for positive and negative
# concentrations, and that if a large enough product value is given rates drop to zero
# rather than reversing
@testset "Test full q function" begin
    @test qs(10.0, 2.0, 46666.666666666664, 1, microbe, 293.15, reactions[1]) ==
          466566.4698895157
    @test qs(10.0, 2.0, 20000.0, 2, microbe, 293.15, reactions[2]) == 199972.50323833904
    @test qs(10.0, 200000.0, 20000.0, 2, microbe, 293.15, reactions[1]) == 0.0
    @test qs(5.0, 2.0, 20000.0, 2, microbe, 293.15, reactions[2]) == 199945.01403634422
    @test qs(10.0, -2.0, 20000.0, 2, microbe, 293.15, reactions[2]) == 199972.50378073016
    @test qs(-10.0, 2.0, 20000.0, 2, microbe, 293.15, reactions[2]) == 0.0
    @test qs(-10.0, -2.0, 20000.0, 2, microbe, 293.15, reactions[2]) == 0.0
end
