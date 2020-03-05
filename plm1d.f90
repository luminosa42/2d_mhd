! 1D piecewise linear method (PLM) update routine based on van Leer's
! MUSCL technique -- derived by stripping out quadratic features from PPM.
! Contact steepening is not included.
!
! This routine updates a 2D uniform Cartesian mesh in either the x or the
! y direction.
!

subroutine plm1d (rho, P, Etot, ux, uy, Nx, Ny, Nbz, dx, dy, dt, gamma, dir)

!------------------------------------------------------------------------------

implicit none

integer, intent(in)   :: Nx, Ny, Nbz
real, intent(in)      :: dx, dy, dt, gamma
character, intent(in) :: dir

real, dimension(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), intent(inout) :: rho, P, ux, uy, &
                                                             Etot

integer :: i, j, iplm
real    :: dtdx, dtdy, dxinv, dyinv, s

real, dimension(1-Nbz:max(Nx,Ny)+Nbz) :: rho0, Etot0, mx, my, ekin, a, &
                                         jrho, jmx, jmy, je, rhoedge, &
                                         uxedge, uyedge, Eedge, Pedge, &
                                         rhoslp, pslp, uxslp, uyslp, &
                                         rhol, pl, uxl, uyl, &
                                         rhor, pr, uxr, uyr, &
                                         rholp, rholm, rholz, &
                                         rhorp, rhorm, rhorz, &
                                         plp, plm, plz, prp, prm, prz, &
                                         uxlp, uxlm, uxlz, uxrp, uxrm, uxrz, &
                                         uylp, uylm, uylz, uyrp, uyrm, uyrz, &
                                         rholtw, rhortw, pltw, prtw, &
                                         uxltw, uxrtw, uyltw, uyrtw, &
                                         altw, artw

real :: chl(-1:1), chr(-1:1), chil(-1:1), chir(-1:1)

real, parameter    :: smallp = 1.E-10, smlrho = 1.E-10, smallu = 1.E-10
real, parameter    :: rieman_tol = 1.E-6
integer, parameter :: nriem = 20

real :: minmod

!------------------------------------------------------------------------------

! If iplm == 0, we ignore linear corrections (ie. we get Godunov's method).

iplm = 1

!------------------------------------------------------------------------------

! Sweeps in the x direction.

if (dir == "x") then

  dxinv = 1./dx
  dtdx  = dt*dxinv

! loop over y-rows

  do j = 1, Ny

! compute current timestep conserved quantities

    rho0(1-Nbz:Nx+Nbz) = rho(1-Nbz:Nx+Nbz,j)
    mx(1-Nbz:Nx+Nbz)   = rho(1-Nbz:Nx+Nbz,j) * ux(1-Nbz:Nx+Nbz,j)
    my(1-Nbz:Nx+Nbz)   = rho(1-Nbz:Nx+Nbz,j) * uy(1-Nbz:Nx+Nbz,j)
    Etot0(1-Nbz:Nx+Nbz)= Etot(1-Nbz:Nx+Nbz,j)

! obtain monotonized slopes

    if (iplm /= 0) then
      do i = -1, Nx+2
        rhoslp(i) = dxinv * minmod( 0.5*(rho(i+1,j)-rho(i-1,j)), &
                                    minmod( 2.*(rho(i,j)-rho(i-1,j)), &
                                            2.*(rho(i+1,j)-rho(i,j)) ) )
        rhoslp(i) = sign( rhoslp(i), rho(i+1,j)-rho(i-1,j) )
        pslp(i)   = dxinv * minmod( 0.5*(P(i+1,j)-P(i-1,j)), &
                                    minmod( 2.*(P(i,j)-P(i-1,j)), &
                                            2.*(P(i+1,j)-P(i,j)) ) )
        pslp(i)   = sign( pslp(i), P(i+1,j)-P(i-1,j) )
        uxslp(i)  = dxinv * minmod( 0.5*(ux(i+1,j)-ux(i-1,j)), &
                                    minmod( 2.*(ux(i,j)-ux(i-1,j)), &
                                            2.*(ux(i+1,j)-ux(i,j)) ) )
        uxslp(i)  = sign( uxslp(i), ux(i+1,j)-ux(i-1,j) )
        uyslp(i)  = dxinv * minmod( 0.5*(uy(i+1,j)-uy(i-1,j)), &
                                    minmod( 2.*(uy(i,j)-uy(i-1,j)), &
                                            2.*(uy(i+1,j)-uy(i,j)) ) )
        uyslp(i)  = sign( uyslp(i), uy(i+1,j)-uy(i-1,j) )
      enddo
    else
      rhoslp(:) = 0.
      pslp(:)   = 0.
      uxslp(:)  = 0.
      uyslp(:)  = 0.
    endif

! compute sound speeds

    do i = -1, Nx+2
      a(i) = sqrt(gamma*P(i,j)/rho(i,j))
    enddo

! compute averages over domains of dependence for all characteristics
! that reach interfaces during timestep

    do i = 0, Nx
      chl(-1) = max( (ux(i,j)-a(i))*dt/dx, 0. )
      chl(0)  = max( ux(i,j)*dt/dx, 0. )
      chl(1)  = max( (ux(i,j)+a(i))*dt/dx, 0. )
      chr(-1) = max( -(ux(i+1,j)-a(i+1))*dt/dx, 0. )
      chr(0)  = max( -ux(i+1,j)*dt/dx, 0. )
      chr(1)  = max( -(ux(i+1,j)+a(i+1))*dt/dx, 0. )

      rholp(i) = rho(i,j)   + 0.5*dx*rhoslp(i)   * (1. - chl(1))
      rholm(i) = rho(i,j)   + 0.5*dx*rhoslp(i)   * (1. - chl(-1))
      rholz(i) = rho(i,j)   + 0.5*dx*rhoslp(i)   * (1. - chl(0))
      rhorp(i) = rho(i+1,j) - 0.5*dx*rhoslp(i+1) * (1. - chr(1))
      rhorm(i) = rho(i+1,j) - 0.5*dx*rhoslp(i+1) * (1. - chr(-1))
      rhorz(i) = rho(i+1,j) - 0.5*dx*rhoslp(i+1) * (1. - chr(0))

      plp(i) = P(i,j)   + 0.5*dx*pslp(i)   * (1. - chl(1))
      plm(i) = P(i,j)   + 0.5*dx*pslp(i)   * (1. - chl(-1))
      plz(i) = P(i,j)   + 0.5*dx*pslp(i)   * (1. - chl(0))
      prp(i) = P(i+1,j) - 0.5*dx*pslp(i+1) * (1. - chr(1))
      prm(i) = P(i+1,j) - 0.5*dx*pslp(i+1) * (1. - chr(-1))
      prz(i) = P(i+1,j) - 0.5*dx*pslp(i+1) * (1. - chr(0))

      uxlp(i) = ux(i,j)   + 0.5*dx*uxslp(i)   * (1. - chl(1))
      uxlm(i) = ux(i,j)   + 0.5*dx*uxslp(i)   * (1. - chl(-1))
      uxlz(i) = ux(i,j)   + 0.5*dx*uxslp(i)   * (1. - chl(0))
      uxrp(i) = ux(i+1,j) - 0.5*dx*uxslp(i+1) * (1. - chr(1))
      uxrm(i) = ux(i+1,j) - 0.5*dx*uxslp(i+1) * (1. - chr(-1))
      uxrz(i) = ux(i+1,j) - 0.5*dx*uxslp(i+1) * (1. - chr(0))

      uylp(i) = uy(i,j)   + 0.5*dx*uyslp(i)   * (1. - chl(1))
      uylm(i) = uy(i,j)   + 0.5*dx*uyslp(i)   * (1. - chl(-1))
      uylz(i) = uy(i,j)   + 0.5*dx*uyslp(i)   * (1. - chl(0))
      uyrp(i) = uy(i+1,j) - 0.5*dx*uyslp(i+1) * (1. - chr(1))
      uyrm(i) = uy(i+1,j) - 0.5*dx*uyslp(i+1) * (1. - chr(-1))
      uyrz(i) = uy(i+1,j) - 0.5*dx*uyslp(i+1) * (1. - chr(0))
    enddo

! "twiddle" quantities are first guesses at correct left and right states

    do i = 0, Nx
      rholtw(i) = rholp(i)
      rhortw(i) = rhorm(i)
      pltw(i)   = plp(i)
      prtw(i)   = prm(i)
      uxltw(i)  = uxlp(i)
      uxrtw(i)  = uxrm(i)
      uyltw(i)  = uylz(i)
      uyrtw(i)  = uyrz(i)
      altw(i)   = sqrt(gamma*pltw(i)/rholtw(i))
      artw(i)   = sqrt(gamma*prtw(i)/rhortw(i))
    enddo

! correct the initial guesses by subtracting characteristic information that
! doesn't make it to the interfaces during the timestep

    do i = 0, Nx
      chil(-1) = (uxltw(i) - uxlm(i) - &
                 (pltw(i)-plm(i))/(rholtw(i)*altw(i))) / (2.*rholtw(i)*altw(i))
      chil(0)  = (pltw(i)-plz(i))/(rholtw(i)*altw(i))**2 + &
                 1./rholtw(i) - 1./rholz(i)
      chil(1)  = -(uxltw(i) - uxlp(i) + &
                 (pltw(i)-plp(i))/(rholtw(i)*altw(i))) / (2.*rholtw(i)*altw(i))
      chir(-1) = (uxrtw(i) - uxrm(i) - &
                 (prtw(i)-prm(i))/(rhortw(i)*artw(i))) / (2.*rhortw(i)*artw(i))
      chir(0)  = (prtw(i)-prz(i))/(rhortw(i)*artw(i))**2 + &
                 1./rhortw(i) - 1./rhorz(i)
      chir(1)  = -(uxrtw(i) - uxrp(i) + &
                 (prtw(i)-prp(i))/(rhortw(i)*artw(i))) / (2.*rhortw(i)*artw(i))

      chl(-1) = max( (ux(i,j)-a(i))*dt*dxinv, 0. )
      chl(0)  = max( ux(i,j)*dt*dxinv, 0. )
      chl(1)  = max( (ux(i,j)+a(i))*dt*dxinv, 0. )
      chr(-1) = max( -(ux(i+1,j)-a(i+1))*dt*dxinv, 0. )
      chr(0)  = max( -ux(i+1,j)*dt*dxinv, 0. )
      chr(1)  = max( -(ux(i+1,j)+a(i+1))*dt*dxinv, 0. )

      if (chl(-1) <= 0.) chil(-1) = 0.
      if (chl(0) <= 0.)  chil(0)  = 0.
      if (chl(1) <= 0.)  chil(1)  = 0.
      if (chr(-1) >= 0.) chir(-1) = 0.
      if (chr(0) >= 0.)  chir(0)  = 0.
      if (chr(1) >= 0.)  chir(1)  = 0.

      rhol(i) = 1. / (1./rholtw(i) - chil(-1) - chil(0) - chil(1))
      rhor(i) = 1. / (1./rhortw(i) - chir(-1) - chir(0) - chir(1))
      pl(i)   = pltw(i) + (rholtw(i)*altw(i))**2 * (chil(-1) + chil(1))
      pr(i)   = prtw(i) + (rhortw(i)*artw(i))**2 * (chir(-1) + chir(1))
      uxl(i)  = uxltw(i) + rholtw(i)*altw(i) * (chil(1) - chil(-1))
      uxr(i)  = uxrtw(i) + rhortw(i)*artw(i) * (chir(1) - chir(-1))
      uyl(i)  = uyltw(i)
      uyr(i)  = uyrtw(i)
    enddo

! call riemann solver to obtain time-averaged, cell-edge quantities

    call riemann (Nx, rhol(1-Nbz), pl(1-Nbz), uxl(1-Nbz), uyl(1-Nbz), &
                  rhor(1-Nbz), pr(1-Nbz), uxr(1-Nbz), uyr(1-Nbz), &
                  rhoedge(1-Nbz), uxedge(1-Nbz), uyedge(1-Nbz), Pedge(1-Nbz), &
                  gamma, smallp, smlrho, smallu, rieman_tol, nriem)

    do i = 0, Nx
      Eedge(i) = 0.5*rhoedge(i)*(uxedge(i)**2 + uyedge(i)**2) + &
                 Pedge(i)/(gamma-1.)
    enddo

! compute time-averaged fluxes

    do i = 0, Nx
      jrho(i) = rhoedge(i) * uxedge(i)
      jmx(i)  = rhoedge(i) * uxedge(i) * uxedge(i) + Pedge(i)
      jmy(i)  = rhoedge(i) * uxedge(i) * uyedge(i)
      je(i)   = (Eedge(i) + Pedge(i)) * uxedge(i)
    enddo

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

    rho0(1-Nbz:Ny+Nbz) = rho(i,1-Nbz:Ny+Nbz)
    mx(1-Nbz:Ny+Nbz)   = rho(i,1-Nbz:Ny+Nbz) * ux(i,1-Nbz:Ny+Nbz)
    my(1-Nbz:Ny+Nbz)   = rho(i,1-Nbz:Ny+Nbz) * uy(i,1-Nbz:Ny+Nbz)
    Etot0(1-Nbz:Ny+Nbz)= Etot(i,1-Nbz:Ny+Nbz)

! obtain monotonized slopes

    if (iplm /= 0) then
      do j = -1, Ny+1
        rhoslp(j) = dyinv * minmod( 0.5*(rho(i,j+1)-rho(i,j-1)), &
                                    minmod( 2.*(rho(i,j)-rho(i,j-1)), &
                                            2.*(rho(i,j+1)-rho(i,j)) ) )
        rhoslp(j) = sign( rhoslp(j), rho(i,j+1)-rho(i,j-1) )
        pslp(j)   = dyinv * minmod( 0.5*(P(i,j+1)-P(i,j-1)), &
                                    minmod( 2.*(P(i,j)-P(i,j-1)), &
                                            2.*(P(i,j+1)-P(i,j)) ) )
        pslp(j)   = sign( pslp(j), P(i,j+1)-P(i,j-1) )
        uxslp(j)  = dyinv * minmod( 0.5*(ux(i,j+1)-ux(i,j-1)), &
                                    minmod( 2.*(ux(i,j)-ux(i,j-1)), &
                                            2.*(ux(i,j+1)-ux(i,j)) ) )
        uxslp(j)  = sign( uxslp(j), ux(i,j+1)-ux(i,j-1) )
        uyslp(j)  = dyinv * minmod( 0.5*(uy(i,j+1)-uy(i,j-1)), &
                                    minmod( 2.*(uy(i,j)-uy(i,j-1)), &
                                            2.*(uy(i,j+1)-uy(i,j)) ) )
        uyslp(j)  = sign( uyslp(j), uy(i,j+1)-uy(i,j-1) )
      enddo
    else
      rhoslp(:) = 0.
      pslp(:)   = 0.
      uxslp(:)  = 0.
      uyslp(:)  = 0.
    endif

! compute sound speeds

    do j = -1, Ny+1
      a(j) = sqrt(gamma*P(i,j)/rho(i,j))
    enddo

! compute averages over domains of dependence for all characteristics
! that reach interfaces during timestep

    do j = 0, Ny
      chl(-1) = max( (uy(i,j)-a(j))*dt*dyinv, 0. )
      chl(0)  = max( uy(i,j)*dt*dyinv, 0. )
      chl(1)  = max( (uy(i,j)+a(j))*dt*dyinv, 0. )
      chr(-1) = max( -(uy(i,j+1)-a(j+1))*dt*dyinv, 0. )
      chr(0)  = max( -uy(i,j+1)*dt*dyinv, 0. )
      chr(1)  = max( -(uy(i,j+1)+a(j+1))*dt*dyinv, 0. )

      rholp(j) = rho(i,j)   + 0.5*dy*rhoslp(j)   * (1. - chl(1))
      rholm(j) = rho(i,j)   + 0.5*dy*rhoslp(j)   * (1. - chl(-1))
      rholz(j) = rho(i,j)   + 0.5*dy*rhoslp(j)   * (1. - chl(0))
      rhorp(j) = rho(i,j+1) - 0.5*dy*rhoslp(j+1) * (1. - chr(1))
      rhorm(j) = rho(i,j+1) - 0.5*dy*rhoslp(j+1) * (1. - chr(-1))
      rhorz(j) = rho(i,j+1) - 0.5*dy*rhoslp(j+1) * (1. - chr(0))

      plp(j) = P(i,j)   + 0.5*dy*pslp(j)   * (1. - chl(1))
      plm(j) = P(i,j)   + 0.5*dy*pslp(j)   * (1. - chl(-1))
      plz(j) = P(i,j)   + 0.5*dy*pslp(j)   * (1. - chl(0))
      prp(j) = P(i,j+1) - 0.5*dy*pslp(j+1) * (1. - chr(1))
      prm(j) = P(i,j+1) - 0.5*dy*pslp(j+1) * (1. - chr(-1))
      prz(j) = P(i,j+1) - 0.5*dy*pslp(j+1) * (1. - chr(0))

      uxlp(j) = ux(i,j)   + 0.5*dy*uxslp(j)   * (1. - chl(1))
      uxlm(j) = ux(i,j)   + 0.5*dy*uxslp(j)   * (1. - chl(-1))
      uxlz(j) = ux(i,j)   + 0.5*dy*uxslp(j)   * (1. - chl(0))
      uxrp(j) = ux(i,j+1) - 0.5*dy*uxslp(j+1) * (1. - chr(1))
      uxrm(j) = ux(i,j+1) - 0.5*dy*uxslp(j+1) * (1. - chr(-1))
      uxrz(j) = ux(i,j+1) - 0.5*dy*uxslp(j+1) * (1. - chr(0))

      uylp(j) = uy(i,j)   + 0.5*dy*uyslp(j)   * (1. - chl(1))
      uylm(j) = uy(i,j)   + 0.5*dy*uyslp(j)   * (1. - chl(-1))
      uylz(j) = uy(i,j)   + 0.5*dy*uyslp(j)   * (1. - chl(0))
      uyrp(j) = uy(i,j+1) - 0.5*dy*uyslp(j+1) * (1. - chr(1))
      uyrm(j) = uy(i,j+1) - 0.5*dy*uyslp(j+1) * (1. - chr(-1))
      uyrz(j) = uy(i,j+1) - 0.5*dy*uyslp(j+1) * (1. - chr(0))
    enddo

! "twiddle" quantities are first guesses at correct left and right states

    do j = 0, Ny
      rholtw(j) = rholp(j)
      rhortw(j) = rhorm(j)
      pltw(j)   = plp(j)
      prtw(j)   = prm(j)
      uxltw(j)  = uxlz(j)
      uxrtw(j)  = uxrz(j)
      uyltw(j)  = uylp(j)
      uyrtw(j)  = uyrm(j)
      altw(j)   = sqrt(gamma*pltw(j)/rholtw(j))
      artw(j)   = sqrt(gamma*prtw(j)/rhortw(j))
    enddo

! correct the initial guesses by subtracting characteristic information that
! doesn't make it to the interfaces during the timestep

    do j = 0, Ny

      chil(-1) = (uyltw(j) - uylm(j) - &
                 (pltw(j)-plm(j))/(rholtw(j)*altw(j))) / (2.*rholtw(j)*altw(j))
      chil(0)  = (pltw(j)-plz(j))/(rholtw(j)*altw(j))**2 + &
                 1./rholtw(j) - 1./rholz(j)
      chil(1)  = -(uyltw(j) - uylp(j) + &
                 (pltw(j)-plp(j))/(rholtw(j)*altw(j))) / (2.*rholtw(j)*altw(j))
      chir(-1) = (uyrtw(j) - uyrm(j) - &
                 (prtw(j)-prm(j))/(rhortw(j)*artw(j))) / (2.*rhortw(j)*artw(j))
      chir(0)  = (prtw(j)-prz(j))/(rhortw(j)*artw(j))**2 + &
                 1./rhortw(j) - 1./rhorz(j)
      chir(1)  = -(uyrtw(j) - uyrp(j) + &
                 (prtw(j)-prp(j))/(rhortw(j)*artw(j))) / (2.*rhortw(j)*artw(j))

      chl(-1) = max( (uy(i,j)-a(j))*dt*dyinv, 0. )
      chl(0)  = max( uy(i,j)*dt*dyinv, 0. )
      chl(1)  = max( (uy(i,j)+a(j))*dt*dyinv, 0. )
      chr(-1) = max( -(uy(i,j+1)-a(j+1))*dt*dyinv, 0. )
      chr(0)  = max( -uy(i,j+1)*dt*dyinv, 0. )
      chr(1)  = max( -(uy(i,j+1)+a(j+1))*dt*dyinv, 0. )

      if (chl(-1) <= 0.) chil(-1) = 0.
      if (chl(0) <= 0.)  chil(0)  = 0.
      if (chl(1) <= 0.)  chil(1)  = 0.
      if (chr(-1) >= 0.) chir(-1) = 0.
      if (chr(0) >= 0.)  chir(0)  = 0.
      if (chr(1) >= 0.)  chir(1)  = 0.

      rhol(j) = 1. / (1./rholtw(j) - chil(-1) - chil(0) - chil(1))
      rhor(j) = 1. / (1./rhortw(j) - chir(-1) - chir(0) - chir(1))
      pl(j)   = pltw(j) + (rholtw(j)*altw(j))**2 * (chil(-1) + chil(1))
      pr(j)   = prtw(j) + (rhortw(j)*artw(j))**2 * (chir(-1) + chir(1))
      uyl(j)  = uyltw(j) + rholtw(j)*altw(j) * (chil(1) - chil(-1))
      uyr(j)  = uyrtw(j) + rhortw(j)*artw(j) * (chir(1) - chir(-1))
      uxl(j)  = uxlz(j)
      uxr(j)  = uxrz(j)

    enddo

! call riemann solver to obtain time-averaged, cell-edge quantities

    call riemann (Ny, rhol(1-Nbz), pl(1-Nbz), uyl(1-Nbz), uxl(1-Nbz), &
                  rhor(1-Nbz), pr(1-Nbz), uyr(1-Nbz), uxr(1-Nbz), &
                  rhoedge(1-Nbz), uyedge(1-Nbz), uxedge(1-Nbz), Pedge(1-Nbz), &
                  gamma, smallp, smlrho, smallu, rieman_tol, nriem)
    do j = 0, Ny
      Eedge(j) = 0.5*rhoedge(j)*(uxedge(j)**2 + uyedge(j)**2) + &
                 Pedge(j)/(gamma-1.)
    enddo

! compute time-averaged fluxes

    do j = 0, Ny
      jrho(j) = rhoedge(j) * uyedge(j)
      jmx(j)  = rhoedge(j) * uxedge(j) * uyedge(j)
      jmy(j)  = rhoedge(j) * uyedge(j) * uyedge(j) + Pedge(j)
      je(j)   = (Eedge(j) + Pedge(j)) * uyedge(j)
    enddo

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

!==============================================================================

! MINMOD function for slope limiting

real function minmod (a, b)

real :: a, b

if (a*b > 0) then

  if (abs(a) < abs(b)) then
    minmod = a
  else
    minmod = b
  endif

else

  minmod = 0.

endif

return
end
