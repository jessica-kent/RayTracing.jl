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
    
    t_max = maximum(params.time)
    sampling_time = collect(range(0,t_max,no_samples));

    indices = Vector{Int64}(undef,no_samples)
    for i = 1:no_samples
        indices[i] = argmin(abs.(params.time .- sampling_time[i]))
    end

    position = zeros(no_samples,2)
    amplitudes = zeros(domain.z_dims, domain.x_dims)
   
    return RayData(position, indices, amplitudes)
end

function initialise_parameters(t_max::Number,no_rays::Int,angle_range::Vector, fields::DomainFields)
    v_max = maximum(1 ./ fields.slowness_field)
    ds = 1
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

function advance_ray(ICs::InitialConditions,params::Parameters,fields::DomainFields,ray::Ray, distance_travelled::Float64)
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
    direction = normalize([dx, dz])
    slowness  = [sx, sz]
    if distance_travelled != 0.0
        amp = ICs.amplitude / sqrt(distance_travelled)
    else
        amp = ICs.amplitude
    end
    
    return Ray(position,direction,amp,slowness)
end


function calculate_ray(sim::RaySimulation, boundary_type::TransmissionType)
    
    domain = sim.domain
    fields = sim.fields
    ICs = sim.ICs
    ray_data = sim.ray_data
    params = sim.params
    ray = initialise_ray(ICs,params)

    if size(ray_data.position,1) != 0
        ray_data.position[1,:] = ICs.position
    end

    steps = size(params.time, 1)
    distance_travelled = 0.0
    amps = ray_data.amplitudes
    transmitted = false

    for i = 2:steps-1
        is_in_domain = in_domain(ray, domain)

        if is_in_domain == true
            at_boundary = boundary_detect(domain, ray; transmitted)
            if at_boundary==true
                nodes = Int.(round.(ray.position))
                nodex = nodes[1]
                nodez = nodes[2]

                amps[nodez,nodex] += ray.amp
                id = domain.domain[nodez,nodex]
                media = get_media(fields, ray.position[2], ray.position[1], domain, id)

                ratio = media[1].c/media[2].c |>real
                
                critical = is_critical(domain, ray, ratio)
                
                if critical == true
                    ray = reflect_ray(domain,params,fields,ray)
                else
                    ray = transmit_ray(domain, ray, params, media)
                end

                transmitted = true
            else
                nodes = Int.(round.(ray.position))
                nodex = nodes[1]
                nodez = nodes[2]
                old_position = ray.position

                ray = advance_ray(ICs,params,fields,ray, distance_travelled)

                distance_travelled += norm(ray.position - old_position)
                transmitted = false
                amps[nodez,nodex] += ray.amp
            end
        else
            # Ray has left the Domain
            break
        end

        if size(ray_data.position,1) != 0
            ray_data = save_ray_data(i,ray,ray_data)
        end    

    end
    return RaySimulationResult(ray, ray_data, amps)
end


function calculate_ray(sim::RaySimulation, boundary_type::ReflectionType)
    
    domain = sim.domain
    fields = sim.fields
    ICs = sim.ICs
    ray_data = sim.ray_data
    params = sim.params

    ray = initialise_ray(ICs,params)

    if size(ray_data.position,1) != 0
        ray_data.position[1,:] = ICs.position
    end

    steps = size(params.time, 1)
    distance_travelled = 0.0
    amps = ray_data.amplitudes
    reflected = false
    for i = 2:steps-1
        is_in_domain = in_domain(ray, domain)
        
        if is_in_domain == true
            nodes = Int.(round.(ray.position))
            nodex = nodes[1]
            nodez = nodes[2]
            amps[nodez,nodex] += ray.amp
            at_boundary = boundary_detect(domain, ray; transmitted = reflected)
            
            nodes = Int.(round.(ray.position))
            nodex = nodes[1]
            nodez = nodes[2]
            if at_boundary==true
                nodes = Int.(round.(ray.position))
                nodex = nodes[1]
                nodez = nodes[2]
                ray = reflect_ray(domain,params,fields,ray)
                reflected = true
                println("position:", ray.position)
                println("direction:", ray.direction)
                println("normal:", domain.boundary_normals[nodez, nodex, :])
            end
            
            old_position = ray.position
            ray = advance_ray(ICs,params,fields,ray, distance_travelled)
            distance_travelled += norm(ray.position - old_position)
            reflected = false
            
        else
            # Ray has left the Domain
            break
        end

        if size(ray_data.position,1) != 0
            ray_data = save_ray_data(i,ray,ray_data)
        end    

    end
    return RaySimulationResult(ray, ray_data, amps)
end
