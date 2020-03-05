! Monte Carlo radiative transfer code
! For ASTR 510
! Yinghe Lu 4/7/2016

program radiation

implicit none

real, parameter :: pi = 3.141592653589793
real, parameter :: r = 1.

integer :: Nph = 1E3, Nbin = 100
integer :: i, j, seed, flag
real	:: x, y, k, ang, alpha, l, theta
real, dimension(101) :: bin, iten

iten(:) = 0.
do i = 1, Nbin+1
    bin(i) = 2*pi*(i-1)/Nbin
end do
    
do j = 1, Nph

    x = 0.
    y = -r
    theta = 0.
    alpha = 0.
    ang = 0.
    l = 0.
    flag = -1

    do while (sqrt(x**2.+y**2.) <= r) ! Intersects

        ! Determine scattering
	k = ran(0)
	if (k <= 0.5) then
	    ang = rand(0)*2.*pi
	    alpha = alpha + ang
	endif

	!call random_seed(seed)
        !call random_number(l)
	l = rand(0)
        x = x+l*cos(alpha)
        y = y+l*sin(alpha)
	!print *, l, alpha, x, y
	flag = flag + 1
    enddo

    ! Revert back
    x = x-l*cos(alpha)
    y = y-l*sin(alpha)
 
    ! If it escapes
    if (x > 0) then
        theta = pi/2. + atan(y/x)
    else if (x < 0) then
	theta = 3*pi/2. + atan(y/x)
    else if (flag > 0) then
	theta = pi
    else
	theta = 0.
    endif
    
    !print *, "photon", j, flag, "scatters to position", x, y, "has theta value of", theta
    do i = 1, Nbin
	if(theta .ne. 0 .and. theta >= bin(i) .and. theta < bin(i+1)) iten(i) = iten(i)+1
    enddo
    
enddo

open (1, file='radout.dat')
do i = 1, Nbin
    write (1,'(ES13.6,ES13.6)') bin(i), iten(i)/Nph
end do
close(1)

end program radiation
