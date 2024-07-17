# ray physics

function snells_law(domain::Domain, ray::Ray)
    #Only returns reflected direction
    
    position = Int.(round.(ray.position))
    normal = domain.boundary_normals[position[2], position[1],:]
    
    #calculate angle from normal
    cosθ1 = dot(normal, ray.direction)


    ref_dir = ray.direction - 2normal*cosθ1
    return ref_dir
end

function reflect_ray(domain::Domain,params::Parameters,fields::DomainFields,ray::Ray)
    dt = params.time[2] - params.time[1]
    position = ray.position
    x = ray.position[1];
    z = ray.position[2]
    v2 = snells_law(domain, ray)

    dx = v2[1]; dz = v2[2];
    position += v2;
    u = interpolate_field(fields.slowness_field, x, z)
    sx = (dx/dt)*u^2; sz = (dz/dt)*u^2;
    
    direction = v2
    slowness  = [sx, sz]
    ref_amp = -ray.amp

    return Ray(position, direction, ref_amp, slowness)

end

