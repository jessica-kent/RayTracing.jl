# ray physics

function snells_law(domain::Domain, ray::Ray, media::Vector{Acoustic{Float64, 2}})
    
    normal = [interpolate_field(domain.boundary_normals[:,:,1], ray.position[1], ray.position[2]), interpolate_field(domain.boundary_normals[:,:,2], ray.position[1], ray.position[2])]
    tangent = [normal[2], -normal[1]]
    cs = [real(m.c) for m in media]

    #calculate angle from normal
    sinθ1 = dot(tangent, ray.direction)

    # snells law
    
    sinθ2 = maximum(cs)*sinθ1/minimum(cs)
    cosθ2 = sign(dot(normal, ray.direction))*sqrt(1.0 - sinθ2^2 + 0.0im) |> real

    trans_dir = cosθ2*normal + sinθ2*tangent
        
    trans_dir = normalize(trans_dir)
    return trans_dir
end

function snells_law(domain::Domain, ray::Ray)
    #Only returns reflected direction
    
    normal = [interpolate_field(domain.boundary_normals[:,:,1], ray.position[1], ray.position[2]), interpolate_field(domain.boundary_normals[:,:,2], ray.position[1], ray.position[2])]
    tangent = [normal[2], -normal[1]]
    
    #calculate angle from normal
    cosθ1 = dot(normal, ray.direction)
    sinθ1 = dot(tangent, ray.direction)


    ref_dir = normalize(-cosθ1*normal + sinθ1*tangent)
    return ref_dir
end

function scattered_amplitudes(domain::Domain, ray::Ray, media::Vector{Acoustic{Float64, 2}})
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

function reflect_ray(domain::Domain,params::Parameters,fields::DomainFields,ray::Ray)
    dt = params.time[2] - params.time[1]
    position = ray.position
    x = ray.position[1];
    z = ray.position[2]


    
    v2 = snells_law(domain, ray)

    dx = v2[1]; dz = v2[2];
    u = interpolate_field(fields.slowness_field, x, z)
    sx = (dx/dt)*u^2; sz = (dz/dt)*u^2;
    
    direction = v2
    slowness  = [sx, sz]
    ref_amp = -ray.amp


    return Ray(position, direction, ref_amp, slowness)

end

function transmit_ray(domain::Domain, ray::Ray, params::Parameters, media::Vector{Acoustic{Float64, 2}})
    cs = [real(m.c) for m in media]
    pos = ray.position
    trans_dir = snells_law(domain, ray, media)

    dt = params.time[2] - params.time[1]
    u2 = 1/cs[2]
    slowness_vec = (u2^2/dt)*trans_dir
    trans_amp = ray.amp
    pos += trans_dir
    transmitted_ray = Ray(pos, trans_dir, trans_amp, slowness_vec)

    return transmitted_ray

end

function get_media(fields::DomainFields, z::Float64, x::Float64, domain::Domain, id::Int)
    
    slowness_field = fields.slowness_field
    ρ_field = fields.ρ_field
    
    ρ_submatrix, nodex, nodez = get_submatrix(ρ_field,z,x)
    slowness_submatrix, nodex, nodez = get_submatrix(slowness_field, z, x)
    dom_submatrix, nodex, nodez = get_submatrix(domain.domain, z, x)

    M = dom_submatrix .- id

    zero_inds = findall(iszero, M)
    one_inds = findall(!iszero, M)

    ρ1s = [ρ_submatrix[i] for i in zero_inds]
    ρ2s = [ρ_submatrix[i] for i in one_inds]
    u1s = [slowness_submatrix[i] for i in zero_inds]
    u2s = [slowness_submatrix[i] for i in one_inds]

    medium1 = Acoustic(mean(ρ1s), 1/mean(u1s), 2)
    medium2 = Acoustic(mean(ρ2s), 1/mean(u2s), 2)

    return [medium1, medium2]
end


