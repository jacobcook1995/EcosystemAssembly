using Test
using Assembly

# WHEN I ADD NEW FUNCTIONS ADD TESTS FOR THEM IN HERE

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
    # TODO - This value seems suspect, so I should look into it further
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
