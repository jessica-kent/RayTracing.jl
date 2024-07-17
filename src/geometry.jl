struct Region
    indices::Vector{CartesianIndex{2}}
    slowness_field::Matrix{Float64}
    ρ_field::Matrix{Float64}
end

struct Boundary
    indices::Vector{CartesianIndex{2}}
    normals::Array{Float64, 3}
    boundary_type::BoundaryType
end


function get_region(domain::Matrix{Float64}, id::Int, fields::DomainFields)
    inds = findall(iszero.(domain .- id), domain)
    full_slowness_field = fields.slowness_field
    full_ρ_field = fields.ρ_field
    slowness_field = zeros(size(full_slowness_field))
    ρ_field = zeros(size(full_ρ_field))

    for i in inds
        slowness_field[i] += full_slowness_field[i]
        ρ_field[i] += full_ρ_field[i]
    end

    return Region(inds, slowness_field, ρ_field)
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
        normals = [-grad_x/3, grad_z/3]
        if normals != [0, 0]
            boundary_normals[j,i,1] = normalize(normals)[1];
            boundary_normals[j,i,2] = normalize(normals)[2];
        end
    end
    boundary_normals[[1,end],:,:] .= boundary_normals[[2,end-1],:,:]
    boundary_normals[:,[1,end],:] .= boundary_normals[:,[2,end-1],:]

    return boundary_normals
end

function get_boundary(domain::Matrix{Float64}, region::Region, boundary_type::BoundaryType)
   
    subdomain = get_subdomain(domain, region)
    nz,nx = size(subdomain)
    boundary = zeros(nz,nx)
    
    for i=2:nz-1,j=2:nx-1
        three_by_three = get_submatrix(subdomain,i,j)
        condition = three_by_three .- subdomain[i,j]
        if iszero(condition)
        else
            boundary[i,j] = 1
        end
    end
    
    indices = findall(!iszero, boundary)
    normals = get_normals(boundary)
    
    return Boundary(indices, normals, boundary_type)
end