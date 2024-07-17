struct PointSource
    position::Vector{Float64}
    amplitude::Float64
    θs::AbstractVector
    slowness::Vector{Float64}
end

function source(s::PointSource)
    position = s.position
    amp = s.amplitude
    θs = s.θs
    directions = [[cos(θ), sin(θ)] for θ in θs]

    rays = [Ray(position, directions[i], amp, s.slowness) for i in 1:length(θs)]

    return rays
end

