! 1D Godunov update routine
!
! This routine updates a 2D uniform Cartesian mesh in either the x or the
! y direction.
!
! Rewritten by Yinghe, from the original Ricker test code

subroutine god1d (rho, P, Etot, ux, uy, Nx, Ny, Nbz, dx, dy, dt, gamma, dir)

!------------------------------------------------------------------------------

implicit none

integer, intent(in)   :: Nx, Ny, Nbz
real, intent(in)      :: dx, dy, dt, gamma
character, intent(in) :: dir

real, dimension(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), intent(inout) :: rho, P, ux, uy, &
                                                             Etot

integer :: i, j
real    :: dtdx, dtdy, dxinv, dyinv, s

real, dimension(1-Nbz:max(Nx,Ny)+Nbz) :: rho0, Etot0, mx, my, ekin, a, &
                                         jrho, jmx, jmy, je, rhoedge, &
                                         uxedge, uyedge, Eedge, Pedge, &
                                         rhol, pl, uxl, uyl, &
                                         rhor, pr, uxr, uyr

real, parameter    :: smallp = 1.E-10, smlrho = 1.E-10, smallu = 1.E-10
real, parameter    :: rieman_tol = 1.E-6
integer, parameter :: nriem = 10

!------------------------------------------------------------------------------

! Sweeps in the x direction.

if (dir == "x") then

  dxinv = 1./dx
  dtdx  = dt*dxinv

! loop over y-rows

  do j = 1, Ny

! compute current timestep conserved quantities

    rho0(:) = rho(:,j)
    mx(:)   = rho(:,j) * ux(:,j)
    my(:)   = rho(:,j) * uy(:,j)
    Etot0(:)= Etot(:,j)

! set up left and right input states for riemann problems

    do i = 0, Nx
      rhol(i) = rho(i,j)
      rhor(i) = rho(i+1,j)
      pl(i)   = P(i,j)
      pr(i)   = P(i+1,j)
      uxl(i)  = ux(i,j)
      uxr(i)  = ux(i+1,j)
      uyl(i)  = uy(i,j)
      uyr(i)  = uy(i+1,j)
    enddo

! call riemann solver to obtain time-averaged, cell-edge quantities

    call riemann (Nx, rhol(1-Nbz), pl(1-Nbz), uxl(1-Nbz), uyl(1-Nbz), &
                  rhor(1-Nbz), pr(1-Nbz), uxr(1-Nbz), uyr(1-Nbz), &
                  rhoedge(1-Nbz), uxedge(1-Nbz), uyedge(1-Nbz), Pedge(1-Nbz), &
                  gamma, smallp, smlrho, smallu, rieman_tol, nriem)
    Eedge = 0.5*rhoedge*(uxedge**2 + uyedge**2) + Pedge/(gamma-1.)

! compute time-averaged fluxes

    jrho(:) = rhoedge(:) * uxedge(:)
    jmx(:)  = rhoedge(:) * uxedge(:) * uxedge(:) + Pedge(:)
    jmy(:)  = rhoedge(:) * uxedge(:) * uyedge(:)
    je(:)   = (Eedge(:) + Pedge(:)) * uxedge(:)

! update

    do i = 1, Nx
      rho(i,j) = rho(i,j) - dtdx*(jrho(i) - jrho(i-1))
      rho(i,j) = max( rho(i,j), smlrho )
      ux(i,j)  = (mx(i) - dtdx*(jmx(i) - jmx(i-1))) / rho(i,j)
      uy(i,j)  = (my(i) - dtdx*(jmy(i) - jmy(i-1))) / rho(i,j)
      Etot(i,j)= Etot0(i) - dtdx*(je(i) - je(i-1))
      ekin(i)  = 0.5 * rho(i,j) * (ux(i,j)**2 + uy(i,j)**2)
      P(i,j)   = max( (Etot(i,j)-ekin(i))*(gamma-1.), smallp )
    enddo

! end loop over y-rows

  enddo

!------------------------------------------------------------------------------

! Sweeps in the y direction.

else

  dyinv = 1./dy
  dtdy  = dt*dyinv

! loop over x-rows

  do i = 1, Nx

! compute current timestep conserved quantities

    rho0(:) = rho(i,:)
    mx(:)   = rho(i,:) * ux(i,:)
    my(:)   = rho(i,:) * uy(i,:)
    Etot0(:)= Etot(i,:)

! set up left and right input states for riemann problems

    do j = 0, Ny
      rhol(j) = rho(i,j)
      rhor(j) = rho(i,j+1)
      pl(j)   = P(i,j)
      pr(j)   = P(i,j+1)
      uxl(j)  = ux(i,j)
      uxr(j)  = ux(i,j+1)
      uyl(j)  = uy(i,j)
      uyr(j)  = uy(i,j+1)
    enddo

! call riemann solver to obtain time-averaged, cell-edge quantities

    call riemann (Ny, rhol(1-Nbz), pl(1-Nbz), uyl(1-Nbz), uxl(1-Nbz), &
                  rhor(1-Nbz), pr(1-Nbz), uyr(1-Nbz), uxr(1-Nbz), &
                  rhoedge(1-Nbz), uyedge(1-Nbz), uxedge(1-Nbz), Pedge(1-Nbz), &
                  gamma, smallp, smlrho, smallu, rieman_tol, nriem)
    Eedge = 0.5*rhoedge*(uxedge**2 + uyedge**2) + Pedge/(gamma-1.)

! compute time-averaged fluxes

    jrho(:) = rhoedge(:) * uyedge(:)
    jmx(:)  = rhoedge(:) * uxedge(:) * uyedge(:)
    jmy(:)  = rhoedge(:) * uyedge(:) * uyedge(:) + Pedge(:)
    je(:)   = (Eedge(:) + Pedge(:)) * uyedge(:)

! update

    do j = 1, Ny
      rho(i,j) = rho(i,j) - dtdy*(jrho(j) - jrho(j-1))
      rho(i,j) = max( rho(i,j), smlrho )
      ux(i,j)  = (mx(j) - dtdy*(jmx(j) - jmx(j-1))) / rho(i,j)
      uy(i,j)  = (my(j) - dtdy*(jmy(j) - jmy(j-1))) / rho(i,j)
      Etot(i,j)= Etot0(j) - dtdy*(je(j) - je(j-1))
      ekin(j)  = 0.5 * rho(i,j) * (ux(i,j)**2 + uy(i,j)**2)
      P(i,j)   = max( (Etot(i,j)-ekin(j))*(gamma-1.), smallp )
    enddo

! end loop over x-rows

  enddo

endif

!------------------------------------------------------------------------------

return
end
