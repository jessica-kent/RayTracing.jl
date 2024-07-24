abstract type Source end

struct PointSource <:Source
    position::Vector{Float64}
    amplitude::Float64
    θs::AbstractVector
    slowness::Float64
end

function source(s::PointSource)
    position = s.position
    amp = s.amplitude
    θs = s.θs
    directions = [[cos(θ), sin(θ)] for θ in θs]

    rays = [Ray(position, d, amp, s.slowness*d) for d in directions]

    return rays
end

