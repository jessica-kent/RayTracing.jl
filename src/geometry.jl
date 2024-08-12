struct Boundary
    boundary_points::AbstractArray
    id::Int
end

function identify_regions(filename::String, threshold::Float64)
    img = load(filename)
    
    segments = unseeded_region_growing(img, threshold)

    domain_ids = segments.segment_labels
    id_map = labels_map(segments)
    return id_map, domain_ids
end

function smooth_domain(domain::Matrix)

    x_dims = size(domain,2)
    z_dims = size(domain,1)
    area = x_dims*z_dims

    α = 2pi*0.5^2/(9*(var(domain)*0.01*area))
    sample_rate = Int(round(sqrt(area)/2))

    gaussians1 = [[exp((-(norm([(i-k)/z_dims,j/x_dims]))^2)/2α) for i in 1:z_dims, j in 1:x_dims] for k in 0:sample_rate:z_dims]

    gaussians2 = [[exp((-(norm([(i-k)/z_dims,(j-sample_rate)/x_dims]))^2)/2α) for i in 1:z_dims, j in 1:x_dims] for k in 0:sample_rate:z_dims]

    gaussians3 = [[exp((-(norm([(i-k)/z_dims,(j-2*sample_rate)/x_dims]))^2)/2α) for i in 1:z_dims, j in 1:x_dims] for k in 0:sample_rate:z_dims]

    gaussians = sum(gaussians1)+sum(gaussians2)+sum(gaussians3)

    dom_spectrum = fft(domain)

    pointwise = [dom_spectrum[i,j]*gaussians[i,j] for i in 1:z_dims, j in 1:x_dims]

    smoothed = ifft(pointwise)

    return abs.(smoothed)
end

function get_normals(boundaries::Matrix)
    nz,nx = size(boundaries)
    boundary_normals = zeros(Float64,nz,nx,2)
    boundaries = smooth_domain(boundaries)
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

function center_of_mass(matrix::AbstractArray)
    total_mass = sum(matrix)
    xs = 1:size(matrix,2)
    zs = 1:size(matrix,1)
    weighted_positions = [matrix[i,j]*[j,i] for i in zs, j in xs]

    center = sum(weighted_positions)/total_mass

    return center
end


function is_closed(matrix::AbstractArray)
    closed = false
    center = Int.(round.(center_of_mass(matrix)))
    x_slice = matrix[:,center[1]]
    z_slice = matrix[center[2],:]
    first_zero_x = findfirst(iszero, x_slice)
    last_one_x = findlast(isone, x_slice)
    first_zero_z = findfirst(iszero, z_slice)
    last_one_z = findlast(isone, z_slice)
    
    if first_zero_x + (size(matrix,2)-last_one_x) -1 != size(matrix,2) && first_zero_z + (size(matrix,1)-last_one_z) -1 != size(matrix,1)
        closed = true
    end
    return closed
end

function march(direction::Vector, center::Vector, matrix::AbstractArray)
    position = center + direction
    at_boundary = false

    while at_boundary == false
        submatrix, nodex, nodez = get_submatrix(matrix, position[2], position[1])

        if sum(submatrix) != 9
            at_boundary = true
        else
            position += direction
        end
    end

    return position
end