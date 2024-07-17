function initialise_ray(ICs::InitialConditions,params::Parameters)
    amp = ICs.amplitude
    x = ICs.position[1];    z = ICs.position[2]

    dt = params.time[2] - params.time[1]
    ds = dt/ICs.u

    sx = ICs.u*sind(ICs.theta);     sz = ICs.u*cosd(ICs.theta);
    dx = (sx/ICs.u)*ds;             dz = (sz/ICs.u)*ds; 
    x += dx;                        z += dz;

    position    = [x , z ]
    direction   = [dx, dz]
    slowness    = [sx, sz]
    return Ray(position,direction,amp,slowness)
end

function initialise_ray_data(domain::Domain, params::Parameters, no_samples::Int64)
    
    dom = domain.binary_domain
    
    t_max = maximum(params.time)
    sampling_time = collect(range(0,t_max,no_samples));

    indices = Vector{Int64}(undef,no_samples)
    for i = 1:no_samples
        indices[i] = argmin(abs.(params.time .- sampling_time[i]))
    end

    position = zeros(no_samples,2)
    amplitudes = zeros(size(dom,1), size(dom,2))

    return RayData(position, indices, amplitudes)
end

function initialise_parameters(t_max::Number,no_rays::Int,angle_range::Vector, fields::DomainFields)
    v_max = maximum(1 ./ fields.slowness_field)
    ds = 0.5
    dt = ds/v_max
    time = collect(0:dt:t_max)
    angles = 0
    if no_rays != 1
        angles = collect(range(angle_range[1],angle_range[2],no_rays))
    else
        angles = [angle_range[1]]
    end
    return Parameters(time,angles)
end

function advance_ray(params::Parameters,fields::DomainFields,ray::Ray)
    x  = ray.position[1];   z  = ray.position[2]
    sx = ray.slowness[1];   sz = ray.slowness[2]

    u    = interpolate_field(fields.slowness_field,x,z)
    dudx = interpolate_field(fields.du_dx,x,z)
    dudz = interpolate_field(fields.du_dz,x,z)

    dt = params.time[2] - params.time[1]
    ds = dt/u

    sx += dudx*ds;                  sz += dudz*ds
    dx = (sx/u)*ds;                 dz = (sz/u)*ds; 
    x += dx;                        z += dz;
    
    position  = [x, z]
    direction = [dx, dz]
    slowness  = [sx, sz]

    amp = ray.amp
    
    return Ray(position,direction,amp,slowness)
end

function calculate_ray(ICs::InitialConditions,domain::Domain,
    params::Parameters,fields::DomainFields,ray_data::RayData)
    
    function save_ray_data(i::Int64,ray::Ray,ray_data::RayData)
        if i in ray_data.indices
            j = findfirst(x -> x == i, ray_data.indices)
            ray_data.position[j:end,:] .= ray.position'
        end
        return ray_data
    end

    ray = initialise_ray(ICs,params)
    if size(ray_data.position,1) != 0
        ray_data.position[1,:] = ICs.position
    end
    steps = size(params.time, 1)

    for i = 2:steps-1
        x = ray.position[1];   z = ray.position[2];

        if 1 <= z <= domain.z_dims && 1 <= x <= domain.x_dims
            # Confirms ray is in domain
            submatrix = get_submatrix(domain.binary_domain,z,x)
            if sum(submatrix) > 2
                # Reflective boundary encountered 
                # needs to be updated
                ray = reflect_ray(domain,params,fields,ray)
            else
                # No reflection, continue as normal
                ray = advance_ray(params,fields,ray)
            end
        else
            # Ray has left the Domain
            break
        end

        if size(ray_data.position,1) != 0
            ray_data = save_ray_data(i,ray,ray_data)
        end    

    end
    return ray, ray_data
end