! Initialization routine for homework #2 (Fortran version)

! Arguments:

!    rho(:,:)     (real) array to receive initial density field
!    P(:,:)       (real) array to receive initial pressure field
!    Etot(:,:)    (real) array to receive initial total energy density field
!    ux(:,:)      (real) array to receive initial x-velocity field
!    uy(:,:)      (real) array to receive initial y-velocity field
!    Nx, Ny       (integer) number of interior zones in x- and y-directions
!    Nbz          (integer) number of boundary zones
!    dx, dy       (real) cell sizes in x- and y-directions
!    gamma        (real) gas ratio of specific heats
!    E            (real) total explosion energy
!    rho_ambient  (real) initial ambient gas density
!    p_ambient    (real) initial ambient gas pressure


subroutine init (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, gamma, &
                 E, rho_ambient, p_ambient, B_ambient, theta)

!------------------------------------------------------------------------------
integer, intent(in)   :: Nx, Ny, Nbz
real, intent(in)      :: dx, dy, gamma, E, rho_ambient, p_ambient, B_ambient, theta

real, dimension(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), &
      intent(inout)                         :: rho, P, ux, uy, Etot, Bx, By

!------------------------------------------------------------------------------

! Initialize the grid

!print *, rho_ambient, p_ambient
rho(:,:) = rho_ambient
P(:,:)   = p_ambient
ux(:,:)  = 0.0
uy(:,:)  = 0.0
Etot(:,:) = 1.0 / (gamma - 1.0) * p_ambient + 0.5*B_ambient**2.
Bx(:,:) = B_ambient*cos(theta)
By(:,:) = B_ambient*sin(theta)

! Create peturbation

!Etot(Nx/2,Ny/2) = E / (dx * dy)
!P(Nx/2,Ny/2) = (gamma - 1.0) * E / (dx * dy)
Etot(Nx/2:Nx/2+1, Ny/2:Ny/2+1) = E / (4 * dx * dy)
P(Nx/2:Nx/2+1, Ny/2:Ny/2+1) = (gamma - 1.0) * E / (4 * dx * dy)

!------------------------------------------------------------------------------

return
end
