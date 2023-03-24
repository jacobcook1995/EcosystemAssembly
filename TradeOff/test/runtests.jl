using Test
using TradeOff

# WHEN I ADD NEW FUNCTIONS ADD TESTS FOR THEM IN HERE

# Test functions defined in TradeOff
@test Q(10.0, 2.0) == 1.0 / 5.0
@test Keq(293.0, 1.0, -1e5) == 28628.16592393438
