


# Macroscopic variable computation
def macroscopic(fin):
    # density
    rho = zeros((nx,ny))
    rho[openPath] = sum(fin[:,openPath], axis=0)
    # velocity
    u = zeros((2, nx, ny))
    for i in range(9):
        u[0,openPath] += v[i,0] * fin[i,openPath]
        u[1,openPath] += v[i,1] * fin[i,openPath]
    u[:,openPath] /= rho[openPath]
    return rho, u

# Equilibrium distribution function (rho = array)
def equilibrium(rho, u):          
    usqr = 3/2 * (u[0,openPath]**2 + u[1,openPath]**2)
    feq = zeros((9,nx,ny))
    for i in range(9):
        cu = 3 * (v[i,0]*u[0,openPath] + v[i,1]*u[1,openPath])
        feq[i,openPath] = rho[openPath]*w[i] * (1 + cu + 0.5*cu**2 - usqr)
    return feq


def dissolveClot(rho, u, K, clot):

    # Determining norm of the velocity
    u_norm = sqrt(u[0]**2+u[1]**2)

    # Computing flow
    flow = rho*u_norm
    
    # Computing amount of dissolution 
    dissolution1 = d*K[0,clot]*flow[clot]
    dissolution2 = d*K[1,clot]*flow[clot]

    # Calculating the resulting K values
    K_tmp1 = K[0,clot] - dissolution1
    K_tmp2 = K[0,clot] - dissolution2

    # Assessing if dissolution is > 0 -> if not set value to 0
    condition1 = K_tmp1 > 0  # Create a boolean mask
    condition2 = K_tmp2 > 0
    filtered_K_tmp1 = where(condition1, K_tmp1, 0) # Check if value validates condition
    filtered_K_tmp2 = where(condition2, K_tmp2, 0) # else replace by 0

    # Updating the dissolved K
    K[0, clot] =  filtered_K_tmp1
    K[1, clot] =  filtered_K_tmp2

    return K