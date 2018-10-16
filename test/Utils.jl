@testset "Utils" begin
    @testset "get" begin
        z = rand(Complex128)

        @getfield z (im, re)
        @test im == imag(z)
        @test re == real(z)

        @getfield z (re, im)
        @test im == imag(z)
        @test re == real(z)

        @getfield z (re, im) (r, i)
        @test i == imag(z)
        @test r == real(z)

        @test_throws ArgumentError (@eval @getfield z z)
        @test_throws ArgumentError (@eval @getfield z (re, im) nothing)
        @test_throws ArgumentError (@eval @getfield z (re, im) (r, i) nothing)
        @test_throws ArgumentError (@eval @getfield z (re, im) (im,) (im,))
        @test_throws AssertionError (@eval @getfield z (re, im) (im,))
end

end
