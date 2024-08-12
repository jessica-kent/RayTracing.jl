# boundary logic

function in_domain(ray::Ray, domain::Domain)
    x = ray.position[1];   z = ray.position[2];
    is_in_domain = true
    if 1 <= z <= domain.z_dims && 1 <= x <= domain.x_dims
        is_in_domain = true
    else
        is_in_domain = false
    end

    return is_in_domain
end

function boundary_detect(domain::Domain, ray::Ray; transmitted::Bool=false)
    x = ray.position[1];   z = ray.position[2];
    at_boundary = false
    submatrix,nodex,nodez = get_submatrix(domain.domain,z,x)
    id = domain.domain[nodez,nodex]
    
    if sum(submatrix) != 9*id && transmitted == false
        at_boundary = true
    else
        at_boundary = false
    end

    return at_boundary
end

function id_check(domain::Domain, ray::Ray, id::Int, ICs::InitialConditions, params::Parameters, fields::DomainFields, distance_travelled::Number)
    new_ray = advance_ray(ICs,params,fields,ray, distance_travelled)
    new_position = new_ray.position
    new_nodes = Int.(round.(new_position))
    new_id = domain.domain[new_nodes[2],new_nodes[1]]
    difference = id - new_id
    same_id = true

    if difference != 0
        same_id = false
    else
        same_id = true
    end
    return same_id
end

function is_critical(domain::Domain, ray::Ray, ratio::Float64)
    critical = false

    position = Int.(round.(ray.position))
    normal = domain.boundary_normals[position[2], position[1],:]
    tangent = [-normal[2], normal[1]]
    
    sinθ1 = dot(tangent, ray.direction)
    angle = asin(sinθ1)

    if ratio < 1.0
        if abs(angle) > asin(ratio)
            critical = true
        else
            critical = false
        end
    else
        critical = false
    end

    return critical
end

function scatter_logic(domain::Domain, fields::DomainFields, ray::Ray, params::Parameters, media::Vector{Acoustic{Float64, 2}}, critical::Bool)
    if critical == false
        rays = scatter_ray(domain, fields, ray, params, media)
    else
        ray = reflect_ray(domain,params,fields,ray)
        rays = ScatteredRays(ray, ray)
    end

    return rays
end

function in_boundary(normal::Vector, direction::Vector)
    is_in_boundary = false

    if dot(normal, direction) == 0.0
        is_in_boundary = true
    end

    return is_in_boundary
end
