# ray physics

function snells_law(domain::Domain, ray::Ray, region1::Region, region2::Region)
    position = Int.(round.(ray.position))
    normal = domain.boundary_normals[position[2], position[1],:]
    
    slowness_field1 = region1.slowness_field
    slowness_field2 = region2.slowness_field
    
    cs = [real(m.c) for m in media]

    #calculate angle from normal
    cosθ1 = dot(normal, ray.direction)
    sinθ1 = sqrt(1.0 - cosθ1^2 + 0.0im) |> real

    # snells law

    sinθ2 = cs[2]*sinθ1/cs[1]
    cosθ2 = sqrt(1.0 - sinθ2^2 + 0.0im) |> real

    trans_dir = [cosθ2, sinθ2]

    return trans_dir
end

function snells_law(domain::Domain, ray::Ray)
    #Only returns reflected direction
    
    position = Int.(round.(ray.position))
    normal = domain.boundary_normals[position[2], position[1],:]
    
    #calculate angle from normal
    cosθ1 = dot(normal, ray.direction)
    sinθ1 = sqrt(1.0 - cosθ1^2 + 0.0im) |> real


    ref_dir = [-cosθ1, sinθ1]
    return ref_dir
end

function scattered_amplitudes(domain::Domain, ray::Ray, region1::Region, region2::Region)
    position = Int.(round.(ray.position))
    normal = domain.boundary_normals[position[2], position[1],:]
    us = [1/real(m.c) for m in media]
    ρs = [m.ρ for m in media]
    incident_amplitude = ray.amp
    
    inc_direction = ray.direction
    trans_direction = snells_law(domain, ray, media)

    scat_matrix_element1 = (us[2]*ρs[1]*dot(trans_direction, normal))/(us[1]*ρs[2]*dot(inc_direction, normal))

    scat_matrix = [1.0 -1.0; scat_matrix_element1 1.0]

    scat_amplitudes = scat_matrix \ [incident_amplitude, incident_amplitude]

    return scat_amplitudes
end

function reflect_ray(domain::Domain,params::Parameters,fields::DomainFields,ray::Ray, boundary_type::BoundaryType)
    dt = params.time[2] - params.time[1]

    x = ray.position;
    v2 = snells_law(domain, ray)

    dx = v2[1]; dz = v2[2];
    x += v2;
    u = interpolate_field(fields.slowness_field, x)
    ρ = interpolate_field(fields.ρ_field, x)
    sx = (dx/dt)*u^2; sz = (dz/dt)*u^2;
    
    position  = x
    direction = v2
    slowness  = [sx, sz]
    
    if boundary_type == TransmissionType()
        media = [Acoustic(1/u, ρ, 2), Acoustic(1/u, ρ, 2)] #placeholder...
        ref_amp = scattered_amplitudes(domain, ray, media)[2]
    else
        ref_amp = -ray.amp
    end

    return Ray(position, direction, ref_amp, slowness)

end

function transmit_ray(domain::Domain, ray::Ray, params::Parameters, media::Vector{Acoustic{Float64, 2}})
    cs = [real(m.c) for m in media]
    pos = ray.position
    trans_dir = snells_law(domain, ray, media)

    dt = params.time[2] - params.time[1]
    u2 = 1/cs[2]
    slowness_vec = (u2^2/dt)*trans_dir
    trans_amp = scattered_amplitudes(domain, ray, media)[1]
    pos += trans_dir
    transmitted_ray = Ray(pos, trans_dir, trans_amp, slowness_vec)

    return transmitted_ray

end

function update_amplitude(ray::Ray, parameters::Parameters)
    #think this function is kinda wrong...
    current_amp = ray.amplitude
    times = parameters.time
    slowness = norm(ray.slowness)

    dt = times[2]-times[1]
    ds = dt/slowness

    f(x) = 1.0 - 0.5*x + 0.75*x^2
    
    scale = f(ds)
    new_amp = scale*current_amp
    
    return new_amp
end

