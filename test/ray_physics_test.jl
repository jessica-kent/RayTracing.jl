using Test
using Test, MultipleScattering, LinearAlgebra, Statistics
@testset "Snells law reflection test" begin
    #test snells law
    mat = zeros(300,300)
    θ0 = rand(LinRange(-pi/2, pi/2, 100))
    θ = rand(θ0:0.1:2*θ0)
    x_grad = cos(θ0)
    z_grad = sin(θ0)
    x_0 = rand(1:300);  
    z_0 = rand(1:300);
    for i in 1:300, j in 1:300
        if z_grad*(i-Int(z_0)) + x_grad*(j-Int(x_0)) < 0
            mat[i,j] += 1.0
        end
    end
    position  = [x_0 , z_0];
    amp_0     = 1;
    domain = Domain(size(mat,1), size(mat, 2), mat, get_normals(mat))
    normal = domain.boundary_normals[Int(round(z_0)), Int(round(x_0)),:]
    tangent = cross(vcat(normal, [0.0]), [0.0, 0.0, 1.0])[1:2]
    direction = [cos(θ), sin(θ)]
    ray = Ray(position, direction, amp_0, rand(2))
    ref_dir = snells_law(domain, ray)

    @test isapprox(dot(ref_dir, -dot(direction, normal)*normal+dot(direction, tangent)*tangent), 1.0; atol = 1e-14)
end