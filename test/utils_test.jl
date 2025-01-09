@testset "utils" begin
    field = 42*ones(512,512)
    x_coords = LinRange(1.0, 512.0, 1000)
    z_coords = LinRange(1.0, 512.0, 1000)

    position = [rand(x_coords), rand(z_coords)]

    interp_val = RayTracing.interpolate_field(field, position[1], position[2])

    @test field - interp_val*ones(512,512) == zeros(size(field))
end