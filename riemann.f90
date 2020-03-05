!! NAME
!! 
!!  riemann
!!
!! DESCRIPTION
!!
!!  Default FLASH Riemann solver, modified to be stand-alone and to work
!!  only with constant gamma.
!!  
!!  Solve Riemann shock tube problem for a general equation of state using 
!!  the method of Colella and Glaz.  Use a two shock approximation, and
!!  linearly interpolation between the head and tail of a rarefaction to
!!  treat rarefactions.
!!
!!  Take as input the effective left and right states, obtained by 
!!  integrating the parabolic reconstructions of the data over the domain
!!  of dependence of each of the characteristics on the respective side
!!  of the interface.  This is accomplished by states().  Return the
!!  solution to Riemann's problem on axis -- this is used in computing
!!  the fluxes.
!!
!!  The Riemann problem for the Euler equations produces 4 states, 
!!  separated by the three characteristics (u - cs, u, u + cs):
!!
!!
!!        l_1      t    l_2       l_3
!!         \       ^   .       /
!!          \  *L  |   . *R   /
!!           \     |  .     /
!!            \    |  .    /
!!        L    \   | .   /    R
!!              \  | .  /
!!               \ |. / 
!!                \|./
!!       ----------+----------------> x
!!
!!       l_1 = u - cs   eigenvalue
!!       l_2 = u        eigenvalue (contact)
!!       l_3 = u + cs   eigenvalue
!!
!!       only density jumps across l_2
!!
!!  references:
!!
!!   CG:   Colella & Glaz 1985, JCP, 59, 264.
!!
!!   CW:   Colella & Woodward 1984, JCP, 54, 174.
!!
!!   Fry:  Fryxell et al. 2000, ApJS, 131, 273.
!!
!!   Toro: Toro 1999, ``Riemann Solvers and Numerical Methods for Fluid
!!         Dynamcs: A Practical Introduction, 2nd Ed.'', Springer-Verlag
!!
!!***

subroutine riemann (nzn, rholft, plft, ulft, utlft, &
                    rhorght, prght, urght, utrght, &
                    rhoav, uav, utav, pav, gamma, &
                    smallp, smlrho, smallu, rieman_tol, nriem)

!==============================================================================
  
  implicit none

  integer, intent(in) :: nzn
  real, intent(in)    :: gamma
  
  real, DIMENSION(nzn+8), intent(in)  :: rholft, plft, ulft, utlft, &
                                         rhorght, prght, urght, utrght
  real, DIMENSION(nzn+8), intent(out) :: rhoav, uav, utav, pav

  real, intent(in)    :: smallp, smlrho, smallu, rieman_tol
  integer, intent(in) :: nriem

  real               :: pstor(nriem+5)

  real, DIMENSION(nzn+8) :: clft, crght, vlft, vrght

  real, DIMENSION(nzn+8) :: wlft, wrght, pstar, ustar, vstar, cestar, &
       rhostr, westar, ps, us, uts, vs, rhos, ces, ws, wes

  real, DIMENSION(nzn+8) :: pstar1, pstar2, wlft1, wrght1, aux
  
  real  ::  ustrl1, ustrr1, ustrl2, ustrr2, delu1, delu2, pres_err
  
  integer :: i, j, k, n, nzn4, ierr
  
  real, DIMENSION(nzn+8) :: scrch1, scrch2, scrch3, scrch4

!==============================================================================

  vlft(:)   = 1./max(rholft, smlrho)
  vrght(:)  = 1./max(rhorght, smlrho)

  clft(:)   = sqrt(gamma*plft*vlft)
  crght(:)  = sqrt(gamma*prght*vrght)

  nzn4 = nzn + 4

  do i = 4, nzn4
    aux(i) = sqrt(0.5 * (gamma-1.)/gamma)
  enddo
  
! construct first guess for secant iteration by assuming that the nonlinear 
! wave speed is equal to the sound speed -- the resulting expression is the
! same as Toro, Eq. 9.28 in the Primitive Variable Riemann Solver (PVRS).
! See also Fry Eq. 72.
  
  do i = 4, nzn4
     pstar1(i) = prght(i) - plft(i) - crght(i) * (urght(i) - ulft(i))
     pstar1(i) = plft(i) + pstar1(i) * (clft(i) / (clft(i) + crght(i)))
     pstar1(i) = max (smallp, pstar1(i))
  enddo


! calculate nonlinear wave speeds for the left and right moving waves based
! on the first guess for the pressure jump.  Again, there is a left and a 
! right wave speed.  Compute this using CG Eq. 34.
  
  do i = 4, nzn4
     wlft1(i)  = sqrt( rholft(i)*(pstar1(i) + &
                                  0.5*(gamma-1.)*(pstar1(i)+plft(i))) )
     wrght1(i) = sqrt( rhorght(i)*(pstar1(i) + &
                                  0.5*(gamma-1.)*(pstar1(i)+prght(i))) )
  enddo


! construct second guess for the pressure using the nonlinear wave speeds
! from the first guess.  This is basically the same thing we did to get
! pstar1, except now we are using the better wave speeds instead of the 
! sound speed.
  
  do i = 4, nzn4
     pstar2(i) = prght(i) - plft(i) - wrght1(i) * (urght(i) - ulft(i))
     pstar2(i) = plft(i) + pstar2(i) * wlft1(i) / (wlft1(i) + wrght1(i))
     pstar2(i) = max (smallp, pstar2(i))
  enddo
  

! begin the secant iteration -- see CG Eq. 17 for details.  We will continue to
! interate for convergence until the error falls below tol (in which case, 
! things are good), or we hit nriem iterations (in which case we have a 
! problem, and we spit out an error).

  do i = 4, nzn4

     pstor(1) = pstar1(i)
     pstor(2) = pstar2(i)

     do n = 1, nriem
        
! new nonlinear wave speeds, using CG Eq. 34
        
        wlft(i)  = sqrt( rholft(i)*(pstar2(i) + &
                                    0.5*(gamma-1.)*(pstar2(i)+plft(i))) )
        wrght(i) = sqrt( rhorght(i)*(pstar2(i) + &
                                    0.5*(gamma-1.)*(pstar2(i)+prght(i))) )
        
! compute the velocities in the "star" state -- using CG Eq. 18 -- ustrl2 and
! ustrr2 are the velocities they define there.  ustrl1 and ustrl2 seem to be
! the velocities at the last time, since pstar1 is the old 'star' pressure, and
! wlft1 is the old wave speed.

        ustrl1    =  ulft(i) - (pstar1(i) -  plft(i)) /  wlft1(i)
        ustrr1    = urght(i) + (pstar1(i) - prght(i)) / wrght1(i)
        ustrl2    =  ulft(i) - (pstar2(i) -  plft(i)) /   wlft(i)
        ustrr2    = urght(i) + (pstar2(i) - prght(i)) /  wrght(i)
        
        delu1     = ustrl1 - ustrr1
        delu2     = ustrl2 - ustrr2
        scrch1(i) = delu2  - delu1
        
        if (abs(pstar2(i)-pstar1(i)) .le. smallp) scrch1(i) = 0.

        if (abs(scrch1(i)) .lt. smallu) then
           delu2 = 0.
           scrch1(i) = 1.
        endif

! pressure at the "star" state -- using CG Eq. 18

        pstar(i)  = pstar2(i) - delu2 * (pstar2(i) - pstar1(i)) / scrch1(i)
        pstar(i)  = max (smallp, pstar(i))

! check for convergence of iteration

        pres_err = abs(pstar(i)-pstar2(i)) / pstar(i)
        if (pres_err .lt. rieman_tol) goto 10
        
! reset variables for next iteration
        
        pstar1(i) = pstar2(i)
        pstar2(i) = pstar(i)
        pstor(n+2) = pstar(i)

        wlft1(i)  = wlft(i)
        wrght1(i) = wrght(i)
        
     enddo

     n = n - 1
     
! print error message and stop code if iteration fails to converge
     
     print *, ' '
     print *, 'Nonconvergence in subroutine rieman!'
     print *, ' '
     print *, 'Zone index       = ', i
     print *, 'Iterations tried = ', n+2
     print *, 'Pressure error   = ', pres_err
     print *, ' '
     print *, 'pL       = ', plft(i),   ' pR       =', prght(i)
     print *, 'uL       = ', ulft(i),   ' uR       =', urght(i)
     print *, 'cL       = ', clft(i),   ' cR       =', crght(i)
     print *, ' '
     print *, 'Iteration history:'
     print *, ' '
     print '(A4, 2X, A20)', 'n', 'p*'
     do j = 1, n+2
       print '(I4, 2X, E20.12)', j, pstor(j)
     enddo
     print *, ' '
     print *, 'Terminating execution.'
     stop

! land here if the iterations have converged

10     continue

  enddo

! end of secant iteration


! calculate fluid velocity for the "star" state -- this comes from the shock
! jump equations, Fry Eq. 68 and 69.  The ustar velocity can be computed
! using either the jump eq. for a left moving or right moving shock -- we use
! the average of the two.

  do i = 4, nzn4
     scrch3(i) = ulft (i) - (pstar(i) -  plft(i)) /  wlft(i)
     scrch4(i) = urght(i) + (pstar(i) - prght(i)) / wrght(i)
     ustar(i)  = 0.5e00 * (scrch3(i) + scrch4(i))
     scrch1(i) = sign(1., ustar(i))
  enddo

! decide which state is located at the zone iterface based on the values 
! of the wave speeds.  This is just saying that if ustar > 0, then the state
! is U_L.  if ustar < 0, then the state on the axis is U_R.

  do i = 4, nzn4

     scrch2(i) = 0.5e00 * ( 1.e00 + scrch1(i))
     scrch3(i) = 0.5e00 * ( 1.e00 - scrch1(i))
     
     ps(i)    = plft(i)   * scrch2(i) + prght(i)  * scrch3(i)
     us(i)    = ulft(i)   * scrch2(i) + urght(i)  * scrch3(i)
     uts(i)   = utlft(i)  * scrch2(i) + utrght(i) * scrch3(i)
     vs(i)    = vlft(i)   * scrch2(i) + vrght(i)  * scrch3(i)
     
     rhos(i)  = 1.e00 / vs(i)
     rhos(i)  = max (smlrho, rhos(i))

     vs(i)    = 1.e00 / rhos(i)
     ws(i)    = wlft(i) * scrch2(i) + wrght(i) * scrch3(i)
     ces(i)   = sqrt (gamma * ps(i) * vs(i))

! compute rhostar, using the shock jump condition (Fry Eq. 80)
     
     vstar(i)  = vs(i) - (pstar(i) - ps(i)) / ws(i) / ws(i)
     rhostr(i) = 1.e00 / vstar(i)
     cestar(i) = sqrt (gamma * pstar(i) * vstar(i))

! compute some factors, Fry Eq. 81 and 82       

     wes(i)    = ces(i)    - scrch1(i) * us(i)
     westar(i) = cestar(i) - scrch1(i) * ustar(i)

     scrch4(i) = ws(i) * vs(i) - scrch1(i) * us(i)
     
     
     if (pstar(i) - ps(i) .ge. 0.e00) then
        wes(i)    = scrch4(i)
        westar(i) = scrch4(i)
     endif
     
  enddo
  
  
! compute correct state for rarefaction fan by linear interpolation

  do i = 4, nzn4
     scrch1(i) = max (wes(i) - westar(i), wes(i) + westar(i), smallu)
     scrch1(i) =     (wes(i) + westar(i)) / scrch1(i)
       
     scrch1(i) = 0.5e00 * (1.e00 + scrch1(i))
     scrch2(i) =           1.e00 - scrch1(i)
       
     rhoav(i)  = scrch1(i) * rhostr(i) + scrch2(i) * rhos (i)
     uav  (i)  = scrch1(i) * ustar(i)  + scrch2(i) * us(i)
     utav (i)  = uts(i)
     pav   (i) = scrch1(i) * pstar(i)  + scrch2(i) * ps(i)

  enddo
    
  do i = 4, nzn4
     if (westar(i) .ge. 0.e00) then
        rhoav(i)  = rhostr(i)
        uav(i)    = ustar(i)
        pav(i)    = pstar(i)
     endif
     
     if (wes(i) .lt. 0.e00) then
        rhoav(i)  = rhos(i)
        uav(i)    = us(i)
        pav(i)    = ps(i)
     endif
  enddo
  
  return
end subroutine riemann
