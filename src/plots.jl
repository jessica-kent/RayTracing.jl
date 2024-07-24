#add recipes later...

function plot_ray(domain::Domain,ray_data::Vector{RayData})
    ray_plot = heatmap(
        domain.domain, 
        yflip=true, 
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

@recipe function plot(d::Domain)
    x_max = d.x_dims
    z_max = d.z_dims
    @series begin
        mat = d.domain
        seriestype --> :heatmap
        colour --> :blues

        (1:x_max, 1:z_max, mat)
    end
end

@recipe function plot(ray_data::RayData)

    @series begin
        x = ray_data.position[:,1]
        y = ray_data.position[:,2]
        (x, y)
    end

end