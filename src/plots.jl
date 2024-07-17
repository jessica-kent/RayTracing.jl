#add recipes later...

function plot_ray(domain::Domain,ray_data::Vector{RayData})
    colour_gradient = cgrad([:white, :black])
    ray_plot = heatmap(
        domain.binary_domain, 
        yflip=true, 
        color=colour_gradient, 
        c=:blues, 
        legend=false, 
        title="Domain",
        xlabel="Horizontal Distance, x",
        ylabel="Vertical Distance, z"
    ) 
    for data in ray_data 
        x = data.position[:,1]; z = data.position[:,2]
        ray_plot = plot!(x,z,yflip=true, linecolor="red")
    end
    display(ray_plot)
end
