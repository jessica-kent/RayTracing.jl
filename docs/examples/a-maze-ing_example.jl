using RayTracing
using Plots

# we can use the image_to_domain function to convert an image to a domain
domain = image_to_domain("dev/RayTracing.jl/docs/examples/maze.png")

# this step ensures that the boundaries are in the right place

domain = Domain(size(domain.domain,1), size(domain.domain, 2), iszero.(domain.domain), -domain.boundary_normals)

# plot domain!
heatmap(domain.domain)

dvdx          = 0.0;
dvdz          = 0.0;
base_velocity = 1.0;        

#define a velocity gradient.

velocity_grads = VelocityGradients(dvdx,dvdz,base_velocity); 

# inialise the fields
fields = initialise_fields(domain, grad_to_field(domain, velocity_grads), zeros(512,512))

# slowness map

heatmap(fields.slowness_field)

#define Parameters

t_max       = 1000;
no_rays     = 100;
angle_range = [0, 360];
params      = RayTracing.initialise_parameters(t_max,no_rays,angle_range,fields);

# initial ray parameters

x_0       = 256.0;  
z_0       = 215.0;
position  = [x_0 , z_0];
amp_0     = 1;
theta     = params.angles[1];
u_initial = RayTracing.interpolate_field(fields.slowness_field, x_0,z_0);

# initials raays and ray data
ray = Ray([513.0, 256.0], rand(2), 1.0, rand(2))
RayTracing.in_domain(ray, domain)
results = Vector{RayTracing.RaySimulationResult}(undef, size(params.angles,1))
no_samples  = 1000

for (i,theta) in enumerate(params.angles)
    ICs      = InitialConditions(position,amp_0,theta,u_initial)
    ray_data = initialise_ray_data(domain, params,no_samples)
    sim = RaySimulation(ICs,domain,params,fields,ray_data)
    results[i] = calculate_ray(sim, ReflectionType())
end

source_data = [r.ray_data for r in results]

plot_ray(domain,source_data)
scatter!([x_0],[z_0])
my_palette = palette(:Oranges, 100)
anim = @animate for i in 1:10:1000
    x = [s.position[1:i,1] for s in source_data]
    z = [s.position[1:i,2] for s in source_data]
    heatmap(domain.domain, c=:blues, leg=false)
    for i in 1:length(x)
        plot!(x[i], z[i],yflip=true, xlims=(0,512), ylims=(0,512), palette=my_palette)
    end
    scatter!([x_0],[z_0])
end

gif(anim,"dev/RayTracing.jl/docs/imgs/a-maze-ing.gif", fps = 7)
