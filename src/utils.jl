
function grad_to_field(domain::Domain, grad::VelocityGradients)
    nx = domain.x_dims; nz = domain.z_dims
    velocity_field = zeros(nz,nx)

    for i = 1:nz, j = 1:nx
        velocity_field[i,j] = grad.base_velocity + i*grad.z_grad + j*grad.x_grad
    end

    slowness_field = 1 ./ velocity_field

    return slowness_field
end

function initialise_fields(domain::Domain,vel::VelocityGradients)
    nx = domain.x_dims; nz = domain.z_dims
    
    slowness_field = grad_to_field(domain, vel)
    du_dx = zeros(nz,nx)
    du_dz = zeros(nz,nx)

    for i = 1:nx-1
        du_dx[:,i] = (slowness_field[:,i+1] .- slowness_field[:,i]) 
    end
    for j = 1:nz-1
        du_dz[j,:] = (slowness_field[j+1,:] .- slowness_field[j,:])
    end
    # Boundary conditions (Dirichlet)
    du_dz[[end],:] = du_dz[[end-1],:]
    du_dx[:,[end]] = du_dx[:,[end-1]]
    
    return DomainFields(slowness_field,du_dx,du_dz)
end

function image_to_domain(image_name::String)
    grey_img = Gray.(load(image_name))
    min_value = minimum(grey_img)
    max_value = maximum(grey_img)
    threshold = 0.7*(max_value-min_value) + min_value
    binary_domain = Int.((grey_img) .< threshold) 
    
    (z_dims,x_dims) = size(binary_domain)
    boundary_normals = get_normals(binary_domain)
    
    return Domain(z_dims,x_dims,binary_domain,boundary_normals)
end

function get_submatrix(matrix::Matrix, z::Float64, x::Float64)
    # Gets 3x3 matrix around points z,x  
    zz = Int(round(z)); xx = Int(round(x))
    if (xx in [1,size(matrix)[2]]) || (zz in [1,size(matrix)[1]])
        submatrix = matrix[zz,xx]
    else
        submatrix = matrix[(zz-1):(zz+1),(xx-1):(xx+1)]
    end
    return submatrix
end

function interpolate_field(field::Matrix, x::Number, z::Number)
    # name ^
    x1 = Int(floor(x)); x2 = Int(ceil(x))
    z1 = Int(floor(z)); z2 = Int(ceil(z))
    if x1 == x || x2 == x 
        interp_col = field[:,Int(x)] 
    else
        col1 = field[:,x1]; col2 = field[:,x2]
        interp_col = col1 .+ (col2 .- col1).*(x-x1)./(x2-x1)
    end
    if z1 == z || z2 == z
        interp_val = interp_col[Int(z)]
    else
        row1 = interp_col[z1]; row2 = interp_col[z2]
        interp_val = row1 + (row2-row1)*(z-z1)/(z2-z1)
    end 
    if isnan(interp_val)
        interp_val = 0
    end

    return interp_val
end
