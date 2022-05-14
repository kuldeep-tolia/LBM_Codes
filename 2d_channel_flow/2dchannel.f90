!LBM code for 2d Flow in Channel

program cf

implicit none

integer, parameter:: n = 2000, m = 80
double precision:: f(0:8,0:n,0:m), feq(0:8,0:n,0:m), u(0:n,0:m), v(0:n,0:m), rho(0:n,0:m), w(0:8), cx(0:8), cy(0:8)
double precision:: uin(0:m), rhow, rhoe, s, sx, sy, omega, t1, t2, u0 = 0.10d0, nu = 0.020, dx = 1.0d0, dy = 1.0d0, dt = 1.0d0, Re
integer:: i, j, k, t, steps = 7200

!Setting parameters

w(:) = (/ 4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 /)
cx(:) = (/ 0.0d0, 1.0d0, 0.0d0, -1.0d0, 0.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0 /)
cy(:) = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0 /)
u(:, :) = 0.0d0
v(:, :) = 0.0d0
rho(:, :) = 1.0d0
omega = 1.0d0 / ((3.0d0 * nu) + 0.50d0)
f(:, :, :) = 0.0d0
feq(:, :, :) = 0.0d0
uin(0) = 0.0d0
uin(1:m-1) = u0
uin(m) = 0.0d0

do i = 0, n
	do j = 0, m
		t1 = u(i,j) * u(i,j) + v(i,j) * v(i,j)
		do k = 0, 8
			t2 = u(i,j) * cx(k) + v(i,j) * cy(k)
			feq(k, i, j) = rho(i,j) * w(k) * (1.0d0 + 3.0d0 * t2 + 4.50d0 * t2 * t2 - 1.50d0 * t1)
			f(k, i, j) = feq(k, i, j)
		end do
	end do
end do
u(0, 1:m-1) = u0

Re = u0 * m / nu
print*, 'Reynolds number = ', Re, 'Omega = ', omega

do t = 1, steps
	!Collision
	do j = 0, m
		do i = 0, n
			t1 = u(i,j) * u(i,j) + v(i,j) * v(i,j)
			do k = 0, 8
				t2 = u(i,j) * cx(k) + v(i,j) * cy(k)
				feq(k, i, j) = rho(i, j) * w(k) * (1.0d0 + 3.0d0 * t2 + 4.50d0 * t2 * t2 - 1.50d0 * t1)
				f(k, i, j) = (1.0d0 - omega) * f(k, i, j) + omega * feq(k, i, j)
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
	!West boundary - velocity inlet
	do j = 0, m
		rhow = f(0, 0, j) + f(2, 0, j) + f(4, 0, j) + 2.0d0 * (f(3, 0, j) + f(6, 0, j) + f(7, 0, j))
		rhow = rhow / (1.0d0 - uin(j))
		f(1, 0, j) = f(3, 0, j) + rhow * uin(j) * 2.0d0 / 3.0d0
		f(5, 0, j) = f(7, 0, j) + rhow * uin(j) / 6.0d0
		f(8, 0, j) = f(6, 0, j) + rhow * uin(j) / 6.0d0
	end do

	!North boundary - wall
	do i = 0, n
		f(4, i, m) = f(2, i, m)
		f(8, i, m) = f(6, i, m)
		f(7, i, m) = f(5, i, m)
	end do
	
	!South boundary - wall
	do i = 0, n
		f(2, i, 0) = f(4, i, 0)
		f(5, i, 0) = f(7, i, 0)
		f(6, i, 0) = f(8, i, 0)
	end do

	!East boundary - outlet
	do j = 0, m
		rhoe = 1.0d0
		u(n, j) = -1.0d0 + ((f(0, n, j) + f(2, n, j) + f(4, n, j) + 2.0d0 * (f(1, n, j) + f(5, n, j) + f(8, n, j))) / rhoe)
		f(3, n, j) = f(1, n, j) - (2.0d0 * rhoe * u(n, j) / 3.0d0)
		f(7, n, j) = f(5, n, j) - (rhoe * u(n, j) / 6.0d0)
		f(6, n, j) = f(8, n, j) - (rhoe * u(n, j) / 6.0d0)
	end do

	do j = 0, m	
		do i = 0, n
			s = 0.0d0
			sx = 0.0d0
			sy = 0.0d0
			do k = 0, 8
				s = s + f(k, i, j)
				sx = sx + f(k, i, j) * cx(k)
				sy = sy + f(k, i, j) * cy(k)
			end do
			rho(i,j) = s	
			u(i, j) = sx / s
			v(i, j) = sy / s
		end do
	end do
	
	do j = 0, m
		v(n, j) = 0.0d0
		v(0, j) = 0.0d0
	end do

	if(mod(t,100) .eq. 0) then
		write(*, 1) t, u(0, m/2), v(0, m/2), rho(0, m/2), u(n/2, m/2), v(n/2, m/2), rho(n/2, m/2), u(n, m/2), v(n, m/2), rho(n, m/2)
		1 format(i5, 9f15.8)
	end if

end do

open(17, file = 'uv_field.txt')
open(37, file = 'u_vel_prof.txt')

do j = 0, m
	do i = 0, n
		write(17, *) i, j, u(i, j), v(i, j), rho(i, j)
	end do
	write(17, *) ''
end do

do j = 0, m
	write(37, *) j / float(m), u(200, j) / u0, u(400, j) / u0, u(1000, j) / u0, u(1600, j) / u0
end do

end program cf
