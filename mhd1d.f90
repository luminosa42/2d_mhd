! 1D MHD update routine
!
! This routine updates a 2D uniform Cartesian mesh in either the x or the
! y direction.
!
! Yinghe Lu 4/7/16

subroutine mhd1d (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, dt, gamma, dir)

!------------------------------------------------------------------------------

implicit none

integer, intent(in)   :: Nx, Ny, Nbz
real, intent(in)      :: dx, dy, dt, gamma
character, intent(in) :: dir

real, dimension(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), intent(inout) :: rho, P, ux, uy, &
                                                             Etot, Bx, By

integer :: i, j
real    :: dtdx, dtdy, dxinv, dyinv, s

real, dimension(1-Nbz:max(Nx,Ny)+Nbz) :: rho0, Etot0, mx, my, jbx, jby, ekin, a, &
                                         jrho, jmx, jmy, je, jB

real, parameter :: smallp = 1.E-10, smlrho = 1.E-10, smallu = 1.E-10

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
    jby(:)   = By(:,j)

! compute time-averaged fluxes

    jrho(:) = rho(:,j) * ux(:,j)
    jmy(:)  = rho(:,j) * ux(:,j) * uy(:,j) - Bx(:,j) * By(:,j)
    jmx(:)  = rho(:,j) * ux(:,j) * ux(:,j) + P(:,j) + 0.5*By(:,j)*By(:,j)
    jB(:)   = ux(:,j) * By(:,j) - uy(:,j) * Bx(:,j)
    je(:)   = (Etot(:,j) + P(:,j)+ 0.5* By(:,j) * By(:,j) ) * ux(:,j) - uy(:,j) * Bx(:,j) * By(:,j)

! Lax-Friedrichs update

    do i = 1, Nx
      rho(i,j) = 0.5*(rho(i+1,j)+rho(i-1,j))- 0.5*dtdx*(jrho(i+1) - jrho(i-1))
      rho(i,j) = max( rho(i,j), smlrho )
      ux(i,j)  = (0.5*(mx(i+1)+mx(i-1))- 0.5*dtdx*(jmx(i+1) - jmx(i-1))) / rho(i,j)
      uy(i,j)  = (0.5*(my(i+1)+my(i-1))- 0.5*dtdx*(jmy(i+1) - jmy(i-1))) / rho(i,j)
      By(i,j)  = 0.5*(jby(i+1)+jby(i-1)) - 0.5*dtdx*(jB(i+1)-jB(i-1))
      Etot(i,j)= 0.5*(Etot0(i+1)+Etot0(i-1)) - 0.5*dtdx*(je(i+1) - je(i-1))
      ekin(i)  = 0.5 * rho(i,j) * (ux(i,j)**2 + uy(i,j)**2)
      P(i,j)   = max( (Etot(i,j)-ekin(i)-0.5*(Bx(i,j)*Bx(i,j)+By(i,j)*By(i,j)) )*(gamma-1.), smallp )
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
    jbx(:)   = Bx(i,:)

! compute time-averaged fluxes

    jrho(:) = rho(i,:) * uy(i,:)
    jmx(:)  = rho(i,:) * ux(i,:) * uy(i,:) - Bx(i,:) * By(i,:)
    jmy(:)  = rho(i,:) * uy(i,:) * uy(i,:) + P(i,:) + 0.5*Bx(i,:)*Bx(i,:)
    jB(:)   = uy(i,:) * Bx(i,:) - ux(i,:) * By(i,:)
    je(:)   = (Etot(i,:) + P(i,:)+ 0.5* Bx(i,:) * Bx(i,:) ) * uy(i,:) - ux(i,:) * Bx(i,:) * By(i,:)

! update

    do j = 1, Ny
      rho(i,j) = 0.5*(rho(i,j+1)+rho(i,j-1))- 0.5*dtdy*(jrho(j+1) - jrho(j-1))
      rho(i,j) = max( rho(i,j), smlrho )
      ux(i,j)  = (0.5*(mx(j+1)+mx(j-1))- 0.5*dtdy*(jmx(j+1) - jmx(j-1))) / rho(i,j)
      uy(i,j)  = (0.5*(my(j+1)+my(j-1))- 0.5*dtdy*(jmy(j+1) - jmy(j-1))) / rho(i,j)
      Bx(i,j)  = 0.5*(jbx(j+1)+jbx(j-1)) - 0.5*dtdy*(jB(j+1)-jB(j-1))
      Etot(i,j)= 0.5*(Etot0(j+1)+Etot0(j-1)) - 0.5*dtdy*(je(j+1) - je(j-1))

      ekin(j)  = 0.5 * rho(i,j) * (ux(i,j)**2 + uy(i,j)**2)
      P(i,j)   = max( (Etot(i,j)-ekin(j)-0.5*(Bx(i,j)*Bx(i,j)+By(i,j)*By(i,j)) )*(gamma-1.), smallp )
    enddo

! end loop over x-rows

  enddo

endif

!------------------------------------------------------------------------------

return
end

