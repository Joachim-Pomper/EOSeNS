using Test
using EOSeNS

@testset "EOSeNS" begin
    @test EOSeNS.NEOS isa Module
    @test EOSeNS.TOV isa Module
    @test EOSeNS.Utils isa Module
    @test NumberDensityDim() == DimensionLessDim() / VolumeDim()
end
