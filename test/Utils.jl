@testset "Utils" begin
    @testset "get" begin
        z = rand(Complex128)

        @get z (im, re)
        @test im == imag(z)
        @test re == real(z)

        @get z (re, im)
        @test im == imag(z)
        @test re == real(z)

        @get z (re, im) (r, i)
        @test i == imag(z)
        @test r == real(z)

        @test_throws ArgumentError (@eval @get z z)
        @test_throws ArgumentError (@eval @get z (re, im) nothing)
        @test_throws ArgumentError (@eval @get z (re, im) (r, i) nothing)
        @test_throws ArgumentError (@eval @get z (re, im) (im,) (im,))
        @test_throws AssertionError (@eval @get z (re, im) (im,))
end

end
