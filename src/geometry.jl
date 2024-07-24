struct Boundary
    indices::Vector{CartesianIndex{2}}
    id::Int
    normals::Array{Float64, 3}
end

function identify_regions(filename::String, threshold::Float64)
    img = load(filename)
    
    segments = unseeded_region_growing(img, threshold)

    domain_ids = segments.segment_labels
    id_map = labels_map(segments)
    return id_map, domain_ids
end

function get_region(domain::Matrix{Float64}, id::Int)
    inds = findall(iszero.(domain .- id), domain)

    return Region(inds)
end

function get_normals(boundaries::Matrix)
    nz,nx = size(boundaries)
    boundary_normals = zeros(Float64,nz,nx,2)
    for i = 2:nx-1, j = 2:nz-1
        grad_z = 0; grad_x = 0
        for ii = (i-1):(i+1)
            grad_z += boundaries[(j-1),ii] - boundaries[(j+1),ii]
        end
        for jj = (j-1):(j+1)
            grad_x += boundaries[jj,(i+1)] - boundaries[jj,(i-1)]
        end
        normals = [grad_x/3, -grad_z/3]
        if normals != [0, 0]
            boundary_normals[j,i,1] = normalize(normals)[1];
            boundary_normals[j,i,2] = normalize(normals)[2];
        end
    end
    boundary_normals[[1,end],:,:] .= boundary_normals[[2,end-1],:,:]
    boundary_normals[:,[1,end],:] .= boundary_normals[:,[2,end-1],:]

    return boundary_normals
end

function get_boundaries(domain::Matrix{Float64})
    
    domain_field = domain.domain
    normals = domain.boundary_normals
    boundaries_field = (abs2.(normals[:,:,1])+abs2.(normals[:,:,2]))/2
    boundary_inds = findall(!iszero, boundaries_field)
    
    for inds in boundary_inds
        boundaries_field[inds] *= domain_field[inds]
    end





end