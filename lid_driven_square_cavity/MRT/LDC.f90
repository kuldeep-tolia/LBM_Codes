!LBM code for Lid driven cavity using Multi Relaxation Time method using D2Q9 discrete velocity model
!Reference -- Lattice Boltzmann Method by A.A.Mohamad
!Re = 1000

subroutine collision(u, v, f, feq, rho, omega, w, rx, ry, n, m, tm, tminv, stminv)

implicit none

integer:: n, m, i, j, k, l
double precision:: f(0:8, 0:n, 0:m), feq(0:8, 0:n, 0:m), u(0:n, 0:m), v(0:n, 0:m), rho(0:n, 0:m)
double precision:: w(0:8), rx(0:8), ry(0:8)
double precision:: tm(0:8, 0:8), tminv(0:8, 0:8), stminv(0:8, 0:8)
double precision:: fmom(0:8, 0:n, 0:m), fmeq(0:8, 0:n, 0:m)
double precision:: omega, s

!Calculate equilibrium moments --> meq
do i = 0, n
	do j = 0, m
		fmeq(0, i, j) = rho(i, j)
		fmeq(1, i, j) = rho(i, j) * (-2.0d0 + 3.0d0 * rho(i, j) * (u(i, j) * u(i, j) + v(i, j) * v(i, j)))
		fmeq(2, i, j) = rho(i, j) * (1.0d0 - 3.0d0 * rho(i, j) * (u(i, j) * u(i, j) + v(i, j) * v(i, j)))
		fmeq(3, i, j) = rho(i, j) * u(i, j)
		fmeq(4, i, j) = -rho(i, j) * u(i, j)
		fmeq(5, i, j) = rho(i, j) * v(i, j)
		fmeq(6, i, j) = -rho(i, j) * v(i, j)
		fmeq(7, i, j) = (rho(i, j)**2) * (u(i, j) * u(i, j) - v(i, j) * v(i, j))
		fmeq(8, i, j) = (rho(i, j)**2) * u(i, j) * v(i, j)
	end do
end do

!Calculate moments --> m
do i = 0, n
	do j = 0, m
		do k = 0, 8
			s = 0.0d0
			do l = 0, 8
				s = s + tm(k, l) * f(l, i, j)
			end do
			fmom(k, i, j) = s
		end do
	end do
end do

!Collision in moment space
do i = 0, n
	do j = 0, m
		do k = 0, 8
			s = 0.0d0
			do l = 0, 8
				s = s + stminv(k, l) * (fmom(l, i, j) - fmeq(l, i, j))
			end do
			f(k, i, j) = f(k, i, j) - s
		end do
	end do
end do

end subroutine collision

subroutine streaming(f, n, m)

implicit none

integer:: n, m, i, j
double precision:: f(0:8, 0:n, 0:m)

do j = 0, m
	do i = n, 1, -1
		f(1, i, j) = f(1, i-1, j)
	end do
	
	do i = 0, n-1
		f(3, i, j) = f(3, i+1, j)
	end do
end do
	
do j = m, 1, -1
	do i = 0, n
		f(2, i, j) = f(2, i, j-1)
	end do
	
	do i = n, 1, -1
		f(5, i, j) = f(5, i-1, j-1) 
	end do
	
	do i = 0, n-1
		f(6, i, j) = f(6, i+1, j-1)
	end do
end do
	
do j = 0, m-1
	do i = 0, n
		f(4, i, j) = f(4, i, j+1)
	end do
		
	do i = 0, n-1
		f(7, i, j) = f(7, i+1, j+1)
	end do
	
	do i = n, 1, -1
		f(8, i, j) = f(8, i-1, j+1)
	end do
end do

end subroutine streaming

subroutine sfbound(f, n, m, uo)

implicit none

integer:: n, m, i, j
double precision:: f(0:8, 0:n, 0:m)
double precision:: uo, rho_n

!Boundary Conditions
do j = 0, m
	!West boundary - wall
	f(1, 0, j) = f(3, 0, j)
	f(5, 0, j) = f(7, 0, j)
	f(8, 0, j) = f(6, 0, j)
	
	!East boundary - wall
	f(3, n, j) = f(1, n, j)
	f(6, n, j) = f(8, n, j)
	f(7, n, j) = f(5, n, j)
end do
	
!South boundary - wall
do i = 0, n
	f(2, i, 0) = f(4, i, 0)
	f(5, i, 0) = f(7, i, 0)
	f(6, i, 0) = f(8, i, 0)
end do
	
!North boundary - moving wall
do i = 1, n-1
	rho_n = f(0, i, m) + f(1, i, m) + f(3, i, m) + (2.0d0 * (f(2, i, m) + f(6, i, m) + f(5, i, m)))
	f(4, i, m) = f(2, i, m)
	f(7, i, m) = f(5, i, m) - (rho_n * uo / 6.0d0)
	f(8, i, m) = f(6, i, m) + (rho_n * uo / 6.0d0)
end do

end subroutine sfbound

subroutine rhouv(u, v, f, rho, rx, ry, n, m)

implicit none

integer:: n, m, i, j, k
double precision:: f(0:8, 0:n, 0:m), rho(0:n, 0:m), u(0:n, 0:m), v(0:n, 0:m) 
double precision:: rx(0:8), ry(0:8)
double precision:: s, sx, sy

do j = 0, m
	do i = 0, n
		s = 0.0d0
		do k = 0, 8
			s = s + f(k, i, j)
		end do
		rho(i, j) = s
	end do
end do
	
do i = 0, n
	rho(i, m) = f(0, i, m) + f(1, i, m) + f(3, i, m) + (2.0d0 * (f(2, i, m) + f(6, i, m) + f(5, i, m)))
end do
	
do j = 0, m
	do i = 0, n
		sx = 0.0d0
		sy = 0.0d0
		do k = 0, 8
			sx = sx + f(k, i, j) * rx(k)
			sy = sy + f(k, i, j) * ry(k)
		end do
		u(i, j) = sx / rho(i, j)
		v(i, j) = sy / rho(i, j)
	end do
end do

end subroutine rhouv

program ldc

implicit none

integer, parameter:: n = 100, m = 100
double precision:: f(0:8, 0:n, 0:m), feq(0:8, 0:n, 0:m), u(0:n, 0:m), v(0:n, 0:m), rho(0:n, 0:m), strf(0:n, 0:m)
double precision:: w(0:8), rx(0:8), ry(0:8)
double precision:: nu = 0.001d0, dx = 1.0d0, dy = 1.0d0, dt = 1.0d0, uerr, verr, uo = 0.05d0, omega, Re
double precision:: tau, s, temp1, temp2
double precision:: tm(0:8, 0:8), tminv(0:8, 0:8), sm(0:8), stminv(0:8, 0:8), idt(0:8, 0:8)
integer:: i, j, k, t, steps = 100000

!Setting parameters
w(:) = (/ 4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 /)
rx(:) = (/ 0.0d0, 1.0d0, 0.0d0, -1.0d0, 0.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0 /)
ry(:) = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0 /)

tm(:, :) = reshape((/ 1.0d0, -4.0d0, 4.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
		& 1.0d0, -1.0d0, -2.0d0, 1.0d0, -2.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, &
		& 1.0d0, -1.0d0, -2.0d0, 0.0d0, 0.0d0, 1.0d0, -2.0d0, -1.0d0, 0.0d0, &
		& 1.0d0, -1.0d0, -2.0d0, -1.0d0, 2.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, &
		& 1.0d0, -1.0d0, -2.0d0, 0.0d0, 0.0d0, -1.0d0, 2.0d0, -1.0d0, 0.0d0, &
		& 1.0d0, 2.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0, 1.0d0, &
		& 1.0d0, 2.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0, 1.0d0, 0.0d0, -1.0d0, &
		& 1.0d0, 2.0d0, 1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, 0.0d0, 1.0d0, &
		& 1.0d0, 2.0d0, 1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0, 0.0d0, -1.0d0 /), (/ 9, 9 /))
		
tminv(:, :) = reshape((/ 4.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, &
		   & -4.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, &
		   & 4.0d0, -2.0d0, -2.0d0, -2.0d0, -2.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, &
		   & 0.0d0, 6.0d0, 0.0d0, -6.0d0, 0.0d0, 6.0d0, -6.0d0, -6.0d0, 6.0d0, &
		   & 0.0d0, -6.0d0, 0.0d0, 6.0d0, 0.0d0, 3.0d0, -3.0d0, -3.0d0, 3.0d0, &
		   & 0.0d0, 0.0d0, 6.0d0, 0.0d0, -6.0d0, 6.0d0, 6.0d0, -6.0d0, -6.0d0, &
		   & 0.0d0, 0.0d0, -6.0d0, 0.0d0, 6.0d0, 3.0d0, 3.0d0, -3.0d0, -3.0d0, &
		   & 0.0d0, 9.0d0, -9.0d0, 9.0d0, -9.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
		   & 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 9.0d0, -9.0d0, 9.0d0, -9.0d0 /), (/ 9, 9 /))

tminv(:, :) = tminv(:, :) / 36.0d0

do i = 0, 8
	do j = 0, 8
		s = 0.0d0
		do k = 0, 8
		
			s = s + tminv(i, k) * tm(k, j)
		end do
		idt(i, j) = s
	end do
end do

omega = 1.0/((3.0d0 * nu) + 0.50d0)
Re = uo * n / nu
tau = 1.0d0 / omega

print*, 'Reynolds number = ', Re, 'Omega = ', omega

sm(:) = (/ 1.0d0, 1.40d0, 1.40d0, 1.0d0, 1.20d0, 1.0d0, 1.20d0, tau, tau /)

do i = 0, 8
	do j = 0, 8
	
		stminv(i, j) = tminv(i, j) * sm(j)
	end do
end do

do j = 0, m
	do i = 0, n
		u(i, j) = 0.0d0
		v(i, j) = 0.0d0
		rho(i, j) = 1.0d0
	end do
end do

do i = 1, n-1
	u(i, m) = uo
	v(i, m) = 0.0d0
end do

do t = 1, steps

	call collision(u, v, f, feq, rho, omega, w, rx, ry, n, m, tm, tminv, stminv)
	call streaming(f, n, m)
	call sfbound(f, n, m, uo)
	call rhouv(u, v, f, rho, rx, ry, n, m)
	
	if(mod(t, 100) .eq. 0) then
		write(*, 1) t, u(0, m/2), v(0, m/2), rho(0, m/2), u(n/2, m/2), v(n/2, m/2), rho(n/2, m/2), u(n, m/2), v(n, m/2), rho(n, m/2)
		1 format(i5, 9f15.8) 
	end if
	
end do

open(17, file = 'uv_field.txt')
open(27, file = 'stream_fucntion.txt')
open(37, file = 'x_mid_plane.txt')
open(47, file = 'y_mid_plane.txt')

strf(0, 0) = 0.0d0
do i = 1, n
	temp1 = 0.50d0 * (rho(i-1, 0) + rho(i, 0))
	strf(i, 0) = strf(i-1, 0) - (temp1 * 0.50d0 * (v(i-1, 0) + v(i, 0)))
end do

do i = 0, n
	do j = 1, m
		temp2 = 0.50d0 * (rho(i, j) + rho(i, j-1))
		strf(i, j) = strf(i, j-1) + (temp2 * 0.50d0 * (u(i, j) + u(i, j-1)))
	end do
end do

do j = 0, m
	do i = 0, n
		write(17, *) i, j, u(i, j), v(i, j), rho(i, j)
		write(27, *) i, j, strf(i, j)
	end do
	write(17, *) ''
	write(27, *) ''
end do

do j = 0, m
	write(37, *) j / float(m), u(n/2, j) / uo, v(n/2, j) / uo
end do

do i = 0, n
	write(47, *) i / float(n), u(i, m/2) / uo, v(i, m/2) / uo
end do

end program ldc
