module RayTracing

#types.jl
export Ray, RayData, BoundaryType, ReflectionType, TransmissionType
export InitialConditions, Parameters, Domain
export DomainFields, VelocityGradients, RaySimulation

#utils.jl

export grad_to_field, initialise_fields, image_to_domain
export get_submatrix, interpolate_field

# sources.jl

export PointSource
export point_source

# ray_tracing.jl

export initialise_ray, initialise_ray_data, initialise_parameters
export advance_ray, calculate_ray

# ray_physics

export snells_law, reflect_ray

#plots.jl

export plot_ray

# geometry.jl

export Region, Boundary
export identify_regions
export get_region, get_normals, get_boundary


using Plots
using RecipesBase
using MultipleScattering
using LinearAlgebra
using Images, ImageSegmentation
using Statistics
using FFTW

include("types.jl")
include("utils.jl")
include("geometry.jl")
include("plots.jl")
include("ray_physics.jl")
include("boundary_logic.jl")
include("sources.jl")
include("ray_tracing.jl")

end
