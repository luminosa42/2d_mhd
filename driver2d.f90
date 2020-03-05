! Driver program for 2D hydrodynamics problems (e.g., Sedov problem in 2D).
! Used for ASTR496CAC homework #2.
!

program driver2d

!------------------------------------------------------------------------------

implicit none

integer          :: Nx, Ny, nmax
real             :: tmax, cfl, gamma, Lx, Ly, E, rho_ambient, p_ambient, tout
real		 :: B_ambient, theta
character(len=7) :: method

real, allocatable, dimension(:,:) :: rho, P, Etot, ux, uy, Bx, By

integer          :: i, j, n, fnum, istat, Nbz
real             :: t, dt, tnext, dx, dy

external         :: plm1d, god1d, mhd1d

!------------------------------------------------------------------------------

! initialize

call read_input_file (Nx, Ny, nmax, tmax, tout, cfl, gamma, Lx, Ly, method, &
                      E, rho_ambient, p_ambient, B_ambient, theta)

if (method == "plm") then
  Nbz = 4  ! really more than we need, expected by riemann solver
else if (method == "godunov") then
  Nbz = 4  ! really more than we need, expected by riemann solver
else if (method == "mhd") then
  Nbz = 2
endif

allocate (rho(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), stat=istat)
if (istat /= 0) then
  write (*,*) 'Could not allocate memory for rho!'
  stop
endif
allocate (P(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), stat=istat)
if (istat /= 0) then
  write (*,*) 'Could not allocate memory for P!'
  stop
endif
allocate (Etot(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), stat=istat)
if (istat /= 0) then
  write (*,*) 'Could not allocate memory for Etot!'
  stop
endif
allocate (ux(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), stat=istat)
if (istat /= 0) then
  write (*,*) 'Could not allocate memory for ux!'
  stop
endif
allocate (uy(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), stat=istat)
if (istat /= 0) then
  write (*,*) 'Could not allocate memory for uy!'
  stop
endif
allocate (Bx(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), stat=istat)
if (istat /= 0) then
  write (*,*) 'Could not allocate memory for Bx!'
  stop
endif
allocate (By(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), stat=istat)
if (istat /= 0) then
  write (*,*) 'Could not allocate memory for By!'
  stop
endif

t = 0.
tnext = tout
fnum = 1
dx = Lx / Nx
dy = Ly / Ny

call init (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, gamma, &
           E, rho_ambient, p_ambient, B_ambient, theta)

call output (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, t, 0)

!------------------------------------------------------------------------------

! step

do n = 1, nmax
  call timestep (rho, P, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, dt, gamma, cfl, &
                 t, tnext, tmax)
  write(*,'(A15,I5,A7,ES13.6,A6,ES13.6)') &
       'Beginning step ', n, ':  t = ', t, ' dt = ', dt
  if (method == "plm") then
    call strang_split (rho, P, Etot, ux, uy, Nx, Ny, Nbz, dx, dy, &
                       t, dt, gamma, plm1d)
  else if (method == "godunov") then
    call strang_split (rho, P, Etot, ux, uy, Nx, Ny, Nbz, dx, dy, &
                       t, dt, gamma, god1d)
  else if (method == "mhd") then
     call strang_split_mhd (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, &
                       t, dt, gamma, mhd1d)
  endif
  if (t >= tnext) then
    call output (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, t, fnum)
    fnum = fnum + 1
    tnext = tnext + tout
  endif
  if (t >= tmax) exit
enddo

!------------------------------------------------------------------------------

! final output

call output (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, t, fnum)

!------------------------------------------------------------------------------

end

!==============================================================================

! Routine to read the input data file containing control parameters for the
! simulation.

! Input file format:

!   number-of-x-zones  number-of-y-zones
!   size-of-x-domain   size-of-y-domain
!   method-to-use ("godunov" or "plm")
!   cfl-parameter
!   gamma
!   maximum-time  maximum-number-of-steps
!   output-interval
!   explosion-energy ambient-density ambient-pressure
!   magnetic-strengh magnetic-angle

subroutine read_input_file (Nx, Ny, nmax, tmax, tout, cfl, gamma, Lx, Ly, &
                            method, E, rho_ambient, p_ambient, B_ambient, theta)

!------------------------------------------------------------------------------

integer, intent(out)          :: Nx, Ny, nmax
real, intent(out)             :: tmax, cfl, gamma, Lx, Ly, E, rho_ambient, &
                                 p_ambient, tout, B_ambient, theta
character(len=7), intent(out) :: method

!------------------------------------------------------------------------------

write (*,*) 'Reading input file inputs.dat...'

open (1, file='inputs.dat', status='old')
read (1,*) Nx, Ny
read (1,*) Lx, Ly
read (1,*) method
read (1,*) cfl
read (1,*) gamma
read (1,*) tmax, nmax
read (1,*) tout
read (1,*) E, rho_ambient, p_ambient
read (1,*) B_ambient, theta
close (1)

if ((method /= "godunov") .and. (method /= "plm") .and. (method /= "mhd")) then
  write (*,*) 'Invalid method specified:  ', method
  stop
endif

!------------------------------------------------------------------------------

return
end

!==============================================================================

! Set boundary conditions.

subroutine boundaries (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, t, dir)

!------------------------------------------------------------------------------

implicit none

integer, intent(in)   :: Nx, Ny, Nbz
real, intent(in)      :: dx, dy, t
character, intent(in) :: dir

real, dimension(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), intent(inout) :: rho, P, ux, uy, &
                                                             Etot, Bx, By

integer :: i, j

!------------------------------------------------------------------------------

if (dir == "x") then

! -x boundary

  do j = 1, Ny
    do i = 1-Nbz, 0
      rho(i,j) = rho(1,j)
      P(i,j)   = P(1,j)
      Etot(i,j)= Etot(1,j)
      ux(i,j)  = ux(1,j)
      uy(i,j)  = uy(1,j)
      Bx(i,j)  = Bx(1,j)
      By(i,j)  = By(1,j)
    enddo
  enddo

! +x boundary

  do j = 1, Ny
    do i = Nx+1, Nx+Nbz
      rho(i,j) = rho(Nx,j)
      P(i,j)   = P(Nx,j)
      Etot(i,j)= Etot(Nx,j)
      ux(i,j)  = ux(Nx,j)
      uy(i,j)  = uy(Nx,j)
      Bx(i,j)  = Bx(Nx,j)
      By(i,j)  = By(Nx,j)
    enddo
  enddo

else

! -y boundary

  do j = 1-Nbz, 0
    do i = 1, Nx
      rho(i,j) = rho(i,1)
      P(i,j)   = P(i,1)
      Etot(i,j)= Etot(i,1)
      ux(i,j)  = ux(i,1)
      uy(i,j)  = uy(i,1)
      Bx(i,j)  = Bx(i,1)
      By(i,j)  = By(i,1)
    enddo
  enddo

! +y boundary

  do j = Ny+1, Ny+Nbz
    do i = 1, Nx
      rho(i,j) = rho(i,Ny)
      P(i,j)   = P(i,Ny)
      Etot(i,j)= Etot(i,Ny)
      ux(i,j)  = ux(i,Ny)
      uy(i,j)  = uy(i,Ny)
      Bx(i,j)  = Bx(i,Ny)
      By(i,j)  = By(i,Ny)
    enddo
  enddo

endif

!------------------------------------------------------------------------------

return
end

!==============================================================================

! Set timestep using CFL criterion.

subroutine timestep (rho, P, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, dt, gamma, cfl, &
                     t, tnext, tmax)

!------------------------------------------------------------------------------

implicit none

integer, intent(in)   :: Nx, Ny, Nbz
real, intent(in)      :: dx, dy, gamma, cfl, t, tnext, tmax
real, intent(out)     :: dt

real, dimension(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), intent(in) :: rho, P, ux, uy, Bx, By

integer :: i, j
real    :: c, dtmin, temp, dxinv, dyinv
real	:: valf, vm, test

!------------------------------------------------------------------------------

dtmin = 1.E99
dxinv = 1./dx
dyinv = 1./dy
temp  = sqrt(dxinv**2 + dyinv**2)

do j = 1, Ny
  do i = 1, Nx
    c = sqrt( max(gamma*P(i,j)/rho(i,j), 0.) )
    valf = sqrt( max(Bx(i,j)**2.+ By(i,j)**2./rho(i,j), 0.) )
    vm = sqrt(c**2.+ valf**2.)
    test = 1./( abs(ux(i,j))*dxinv + abs(uy(i,j))*dyinv + vm*temp )
    dtmin = min(dtmin,test )
  enddo
enddo
dt = cfl*dtmin

! Make sure we pull up "exactly" to checkpoints and the final time, but go
! just a little bit past so that we trigger the appropriate if-tests and
! move past these points.

if ( t+2*dt > tnext ) dt = 0.5*abs(tnext*(1.+1.E-6) - t)
if ( t+2*dt > tmax ) dt = 0.5*abs(tmax*(1.+1.E-6) - t)

!------------------------------------------------------------------------------

return
end

!==============================================================================

! Write an output file.

! Output file format:

! #   t = time
! #   Nx= number-of-x-zones
! #   Ny= number-of-y-zones
! #   dx= x-zone-width
! #   dy= y-zone-width
! #
! #    i     j     rho     P     rho*E   ux     uy
!    1    1  rho(1,1)  P(1,1)  Etot(1,1)  ux(1,1)  uy(1,1)
!    2    1  rho(2,1)  P(2,1)  Etot(2,1)  ux(2,1)  uy(2,1)
!    ...
!    Nx   1  rho(Nx,1) P(Nx,1) Etot(Nx,1) ux(Nx,1) uy(Nx,1)
! (blank line)
!    1    2  rho(1,2)  P(1,2)  Etot(1,2)  ux(1,2)  uy(1,2)
!    2    2  rho(2,2)  P(2,2)  Etot(2,2)  ux(2,2)  uy(2,2)
!    ...
!    Nx   2  rho(Nx,2) P(Nx,2) Etot(Nx,2) ux(Nx,2) uy(Nx,2)
! (blank line)
!    1    3  rho(1,3)  P(1,3)  Etot(1,3)  ux(1,3)  uy(1,3)
!    ...
!    ...
!    Nx   Ny rho(Nx,Ny) P(Nx,Ny) Etot(Nx,Ny) ux(Nx,Ny) uy(Nx,Ny)


subroutine output (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, t, fnum)

!------------------------------------------------------------------------------

implicit none

integer, intent(in)   :: Nx, Ny, Nbz, fnum
real, intent(in)      :: dx, dy, t

real, dimension(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), intent(in) :: rho, P, ux, uy, Bx, By, Etot

integer          :: i, j
character(len=4) :: fnumstr

!------------------------------------------------------------------------------

write (fnumstr(1:4), '(I4.4)') fnum

write (*,*) 'Writing output file ', fnumstr, '...'

open (1, file='output_mhd_128_' // fnumstr)
write (1,'(A8,ES13.6)') '#   t= ', t
write (1,'(A8,I5)') '#   Nx= ', Nx
write (1,'(A8,I5)') '#   Ny= ', Ny
write (1,'(A8,ES13.6)') '#   dx= ', dx
write (1,'(A8,ES13.6)') '#   dy= ', dy
write (1,'(A1)') '#'
write (1,'(A)') '#   i     j     rho     P     rho*E   ux     uy	Bx	By'
do j = 1, Ny
  do i = 1, Nx
    write (1,'(I4,1X,I4,7(2X,E13.6))') i, j, rho(i,j), P(i,j), Etot(i,j), &
                                       ux(i,j), uy(i,j), Bx(i,j), By(i,j)
  enddo
  write (1,*)
enddo
close (1)

!------------------------------------------------------------------------------

return
end

!==============================================================================

! Use Strang splitting to create a 2D method from a 1D method

subroutine strang_split (rho, P, Etot, ux, uy, Nx, Ny, Nbz, dx, dy, &
                         t, dt, gamma, method1d)

!------------------------------------------------------------------------------

implicit none

integer, intent(in)   :: Nx, Ny, Nbz
real, intent(in)      :: dx, dy, dt, gamma
real, intent(inout)   :: t

real, dimension(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), &
      intent(inout)                         :: rho, P, ux, uy, Etot

real, dimension(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz) :: Bx, By

external method1d

!------------------------------------------------------------------------------
Bx(:,:) = 0.
By(:,:) = 0.

call boundaries (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, t, 'x')
call method1d (rho, P, Etot, ux, uy, Nx, Ny, Nbz, dx, dy, dt, gamma, 'x')
call boundaries (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, t, 'y')
call method1d (rho, P, Etot, ux, uy, Nx, Ny, Nbz, dx, dy, dt, gamma, 'y')
t = t + dt

call boundaries (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, t, 'y')
call method1d (rho, P, Etot, ux, uy, Nx, Ny, Nbz, dx, dy, dt, gamma, 'y')
call boundaries (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, t, 'x')
call method1d (rho, P, Etot, ux, uy, Nx, Ny, Nbz, dx, dy, dt, gamma, 'x')
t = t + dt

!------------------------------------------------------------------------------

return
end

subroutine strang_split_mhd(rho, P, Etot, ux, uy, Bx, By,Nx, Ny, Nbz, dx, dy, &
                         t, dt, gamma, method1d)

!------------------------------------------------------------------------------

implicit none

integer, intent(in)   :: Nx, Ny, Nbz
real, intent(in)      :: dx, dy, dt, gamma
real, intent(inout)   :: t

real, dimension(1-Nbz:Nx+Nbz,1-Nbz:Ny+Nbz), &
      intent(inout)                         :: rho, P, ux, uy, Bx, By, Etot

external method1d

!------------------------------------------------------------------------------


call boundaries (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, t, 'x')
call method1d (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, dt, gamma, 'x')
call boundaries (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, t, 'y')
call method1d (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, dt, gamma, 'y')
t = t + dt

call boundaries (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, t, 'y')
call method1d (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, dt, gamma, 'y')
call boundaries (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, t, 'x')
call method1d (rho, P, Etot, ux, uy, Bx, By, Nx, Ny, Nbz, dx, dy, dt, gamma, 'x')
t = t + dt

!------------------------------------------------------------------------------

return
end

!==============================================================================
