!LBM code for Couette flow
!Reference -- Lattice Boltzmann Method by A.A.Mohamad

program cf

implicit none

integer, parameter:: n = 200, m = 60
double precision:: f(0:8, 0:n, 0:m), feq(0:8, 0:n, 0:m), u(0:n, 0:m), v(0:n, 0:m), rho(0:n, 0:m)
double precision:: w(0:8), ex(0:8), ey(0:8)
double precision:: omega, rho_n, ut
integer:: i, j, k, t, steps = 40000
double precision:: nu, dx = 1.0d0, dy = 1.0d0, dt = 1.0d0, rho0 = 1.0d0, tau = 1.00d0, t1, t2, t3, t4, s, sx, sy
double precision:: uexact(0:m), utop(0:n)

!Setting parameters
w(:) = (/ 4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 /)
ex(:) = (/ 0.0d0, 1.0d0, 0.0d0, -1.0d0, 0.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0 /)
ey(:) = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0 /)

u(:, :) = 0.0d0
v(:, :) = 0.0d0
rho(:, :) = rho0

nu = (tau - 0.50d0) / 3.0d0
omega = 1.0d0/((3.0d0 * nu) + 0.50d0)

utop(:) = 0.10d0
u(:, m) = utop
ut = 0.10d0

!Initialising distribution function
do j = 0, m
	do i = 0, n
		t1 = u(i, j) * u(i, j) + v(i, j) * v(i, j)
		do k = 0, 8
			t2 = u(i, j) * ex(k) + v(i, j) * ey(k)
			feq(k, i, j) = rho(i, j) * w(k) * (1.0d0 + (3.0d0 * t2) + (4.50d0 * t2 * t2) - (1.50d0 * t1))
			f(k, i, j) = feq(k, i, j)
		end do
	end do
end do

do t = 1, steps
	
	!Collission
	do j = 0, m
		do i = 0, n
			t1 = u(i, j) * u(i, j) + v(i, j) * v(i, j)
			do k = 0, 8
				t2 = u(i, j) * ex(k) + v(i, j) * ey(k)
				feq(k, i, j) = rho(i, j) * w(k) * (1.0d0 + (3.0d0 * t2) + (4.50d0 * t2 * t2) - (1.50d0 * t1))
				f(k, i, j) = (omega * feq(k, i, j)) + ((1.0d0 - omega) * f(k, i, j))
			end do
		end do
	end do
	
	!Streaming
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
	
	!Boundary Conditions
	!South boundary - wall
	do i = 0, n
		f(2, i, 0) = f(4, i, 0)
		f(5, i, 0) = f(7, i, 0)
		f(6, i, 0) = f(8, i, 0)
	end do
	
	do j = 0, m
		!West boundary - periodicity
		f(1, 0, j) = f(1, n, j)
		f(5, 0, j) = f(5, n, j)
		f(8, 0, j) = f(8, n, j)
		
		!East boundary - periodicity
		f(3, n, j) = f(3, 0, j)
		f(6, n, j) = f(6, 0, j)
		f(7, n, j) = f(7, 0, j)
	end do
	
	!North boundary - moving wall
	do i = 0, n
		rho_n = f(0, i, m) + f(1, i, m) + f(3, i, m) + (2.0d0 * (f(2, i, m) + f(6, i, m) + f(5, i, m)))
		f(4, i, m) = f(2, i, m)
		f(7, i, m) = f(5, i, m) - (0.50d0 * rho_n * utop(i))
		f(8, i, m) = f(6, i, m) + (0.50d0 * rho_n * utop(i))
	end do
	
	do j = 0, m
		do i = 0, n
			s = 0.0d0
			do k = 0, 8
				s = s + f(k, i, j)
			end do
			rho(i, j) = s
		end do
	end do
	
	do j = 0, m
		do i = 0, n
			sx = 0.0d0
			sy = 0.0d0
			do k = 0, 8
				sx = sx + f(k, i, j) * ex(k)
				sy = sy + f(k, i, j) * ey(k)
			end do
			u(i, j) = (sx / rho(i, j)) 
			v(i, j) = sy / rho(i, j)
		end do
	end do
	
	u(:, m) = utop
	v(:, m) = 0.0d0
	u(:, 0) = 0.0d0
	v(:, 0) = 0.0d0
	
	if(mod(t, 100) .eq. 0) then
		write(*, 1) t, u(160, 0), rho(160, 0), u(160, 20), rho(160, 20), u(160, 40), rho(160, 40), u(160, m), rho(160, m)
		1 format(i5, 9f15.8) 
	end if	
	
end do

do j = 0, m
	uexact(j) = j * ut / m
end do

open(17, file = 'uv_field.txt')
open(27, file = 'u_prof.txt')

do j = 0, m
	do i = 0, n
		write(17, *) i, j, u(i, j), v(i, j), rho(i, j)
	end do
	write(17, *) ''
end do

do j = 0, m
	write(27, *) j / float(m), u(150, j), v(150, j), uexact(j)
end do

end program cf
	
	
	
	
	
				

