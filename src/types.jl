struct Ray
    position::Vector{Float64}
    direction::Vector{Float64}
    amp::Float64
    slowness::Vector{Float64}
end

abstract type BoundaryType end

struct TransmissionType <: BoundaryType end

struct ReflectionType <: BoundaryType end

struct RayData 
    position::Matrix{Float64}
    indices::Vector{Int64}
    amplitudes::Matrix{Float64}
end

struct InitialConditions
    position::Vector{Float64}
    amplitude::Float64
    theta::Float64
    u::Float64
end

struct Parameters
    time::Vector{Float64}
    angles::Vector{Float64}
end

struct Domain
   z_dims::Int64
   x_dims::Int64
   binary_domain::Matrix{Bool}
   boundary_normals::Array{Float64, 3}
end

struct DomainFields
    slowness_field::Matrix{Float64}
    du_dx::Matrix{Float64}
    du_dz::Matrix{Float64}
end

struct VelocityGradients
    x_grad::Float64
    z_grad::Float64
    base_velocity::Float64
end

struct RaySimulation
    domain::Domain
    ray::Ray
    parameters::Parameters
end