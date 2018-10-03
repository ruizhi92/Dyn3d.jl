@testset "Utils" begin
    @testset "get" begin
        z = rand(Complex128)

        @get_field z (im, re)
        @test im == imag(z)
        @test re == real(z)

        @get_field z (re, im)
        @test im == imag(z)
        @test re == real(z)

        @get_field z (re, im) (r, i)
        @test i == imag(z)
        @test r == real(z)

        @test_throws ArgumentError (@eval @get_field z z)
        @test_throws ArgumentError (@eval @get_field z (re, im) nothing)
        @test_throws ArgumentError (@eval @get_field z (re, im) (r, i) nothing)
        @test_throws ArgumentError (@eval @get_field z (re, im) (im,) (im,))
        @test_throws AssertionError (@eval @get_field z (re, im) (im,))
end

end
