using RayTracing
using Test, MultipleScattering, LinearAlgebra, Statistics

@testset "Normal calculation test" begin
    mat = zeros(300,300)
    θ0 = pi/6
    x_grad = cos(θ0)
    z_grad = sin(θ0)
    x_0 = 150;  
    z_0 = 150;
    for i in 1:300, j in 1:300
        if z_grad*(i-Int(z_0)) + x_grad*(j-Int(x_0)) < 0
            mat[i,j] += 1.0
        end
    end

    normals = get_normals(mat)
    normal = normals[x_0, z_0,:]

    reference_normal = -[cos(θ0), sin(θ0)]
    norm(normal - reference_normal)
    

    @test isapprox(norm(normal - reference_normal), 0.0; atol = 1e-2)
end

#not very robust, need a better error that 1e-2!!!
#run some tests to see what improves error the most.