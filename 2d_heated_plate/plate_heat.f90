!LBM code for diffusion problem in 2D plate using D2Q9

program heated_plate

implicit none

	integer, parameter :: n = 800, m = 400
	double precision :: f(0:8, 0:n, 0:m), feq(0:8, 0:n, 0:m), theta(0:n, 0:m)
	double precision :: w(0:8), cx(0:8), cy(0:8)
	double precision :: dx = 1.0d0, dy = 1.0d0, dt = 1.0d0, omega1, omega2, omega_i, twall1 = 1.0d0, twall2 = 1.0d-1
	double precision :: s, alpha2 = 2.50d-1, alpha1, alpha_i, alpha_ratio = 1.0d1
	integer :: i, j, k, t, steps = 200000
	
	!Setting parameters
	w(:) = (/ 4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 /)
	cx(:) = (/ 0.0d0, 1.0d0, 0.0d0, -1.0d0, 0.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0 /)
	cy(:) = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0 /)
	
	!Initialising temperature field
	theta(:, :) = 0.0d0
	
	!Imposing wall boundary condition
	theta(0, :) = twall1
	theta(n, :) = twall2
	
	alpha1 = alpha_ratio * alpha2
	alpha_i = (2.0d0 * alpha1 * alpha2) / (alpha1 + alpha2)
	omega1 = 1.0d0 / ((3.0d0 * alpha1) + 5.0d-1)
	omega2 = 1.0d0 / ((3.0d0 * alpha2) + 5.0d-1)
	omega_i = 1.0d0 / ((3.0d0 * alpha_i) + 5.0d-1)		!Method-1 to calculate omega_i
!	omega_i = (2.0d0 * omega1 * omega2) / (omega1 + omega2)	!Method-2 to calculate omega_i
	print*, "Value of Omega1 = ", omega1
	print*, "Value of Omega2 = ", omega2
	print*, "Value of Omega at interface = ", omega_i
	
	!Initialising density and equilibrium distribution functions
	do j = 0, m
		do i = 0, n
			do k = 0, 8
				f(k, i, j) = w(k) * theta(i, j)
				feq(k, i, j) = f(k, i, j)
			end do
		end do
	end do
	
	do t = 1, steps
	
		!Calculating macro-variable theta
		do j = 0, m
			do i = 0, n
				s = 0.0d0
				do k = 0, 8
					s = s + f(k, i, j)
				end do
				theta(i, j) = s
			end do
		end do
		
		!Collision
		do j = 0, m
			do i = 0, (n/2 - 1)
				do k = 0 , 8
					feq(k, i, j) = w(k) * theta(i, j) 
					f(k, i, j) = (omega1 * feq(k, i, j)) + ((1.0d0 - omega1) * f(k, i, j))
				end do
			end do
			
			do k = 0, 8
				feq(k, n/2, j) = w(k) * theta(n/2, j)
				f(k, n/2, j) = (omega_i * feq(k, n/2, j)) + ((1.0d0 - omega_i) * f(k, n/2, j))
			end do
			
			do i = (n/2 + 1), n
				do k = 0, 8
					feq(k, i, j) = w(k) * theta(i, j) 
					f(k, i, j) = (omega2 * feq(k, i, j)) + ((1.0d0 - omega2) * f(k, i, j))
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
		
		!Boundary Conditions -- for problem 1
!		do j = 0, m
			!Left boundary - Constant temperature (theta = 1)
!			f(1, 0, j) = w(1) * twall + w(3) * twall - f(3, 0, j)
!			f(5, 0, j) = w(5) * twall + w(7) * twall - f(7, 0, j)
!			f(8, 0, j) = w(8) * twall + w(6) * twall - f(6, 0, j)
						
			!Right boundary - Constant temperatre (theta = 0)
!			f(3, n, j) = -1.0d0 * f(1, n, j)
!			f(7, n, j) = -1.0d0 * f(5, n, j)
!			f(6, n, j) = -1.0d0 * f(8, n, j)
!		end do
		
!		do i = 0, n
			!Top boundary - Constant temperature (theta = 0)
!			f(4, i, m) = -1.0d0 * f(2, i, m)
!			f(7, i, m) = -1.0d0 * f(5, i, m)
!			f(8, i, m) = -1.0d0 * f(6, i, m)
			
			!Bottom boundary - Adiabatic condition
!			do k = 0, 8
!				f(k, i, 0) = f(k, i, 1)
!			end do
!		end do

		!Boundary Conditions -- for problem 2
!		do j = 0, m
			!Left boundary - Constant temperature (theta = 0)
!			f(1, 0, j) = -1.0d0 * f(3, 0, j)
!			f(5, 0, j) = -1.0d0 * f(7, 0, j)
!			f(8, 0, j) = -1.0d0 * f(6, 0, j)
			
			!Right boundary -- Constant temperature (theta = 0)
!			f(3, n, j) = -1.0d0 * f(1, n, j)
!			f(7, n, j) = -1.0d0 * f(5, n, j)
!			f(6, n, j) = -1.0d0 * f(8, n, j)
!		end do
		
!		do i = 0, n
			!Top boundary -- Constant temperature (theta = twall)
!			f(4, i, m) = w(4) * twall + w(2) * twall - f(2, i, m)
!			f(7, i, m) = w(7) * twall + w(5) * twall - f(5, i, m)
!			f(8, i, m) = w(8) * twall + w(6) * twall - f(6, i, m)
			
			!Bottom boundary -- Constant temperature (theta = 0)
!			f(2, i, 0) = -1.0d0 * f(4, i, 0)
!			f(5, i, 0) = -1.0d0 * f(7, i, 0)
!			f(6, i, 0) = -1.0d0 * f(8, i, 0)
!		end do
		
		!Boundary Conditions -- for problem 3 (k_A / k_B = 10)
		do j = 0, m
			!Left boundary - Constant temperature (theta = 1)
			f(1, 0, j) = w(1) * twall1 + w(3) * twall1 - f(3, 0, j)
			f(5, 0, j) = w(5) * twall1 + w(7) * twall1 - f(7, 0, j)
			f(8, 0, j) = w(8) * twall1 + w(6) * twall1 - f(6, 0, j)
			
			!Right boundary -- Constant temperature (theta = 0.1)
			f(3, n, j) = w(3) * twall2 + w(1) * twall2 - f(1, n, j)
			f(7, n, j) = w(7) * twall2 + w(5) * twall2 - f(5, n, j)
			f(6, n, j) = w(6) * twall2 + w(8) * twall2 - f(8, n, j)
		end do
		
		do i = 0, n
			!Top boundary -- Adiabatic condition
			do k = 0, 8
				f(k, i, m) = f(k, i, m-1)
			end do
			
			!Bottom boundary -- Constant temperature (theta = 0)
			f(2, i, 0) = -1.0d0 * f(4, i, 0)
			f(5, i, 0) = -1.0d0 * f(7, i, 0)
			f(6, i, 0) = -1.0d0 * f(8, i, 0)
		end do
			
		
		!Printing some theta values
		if (mod(t, 10000) .eq. 0) then
			write(*, 1) t, theta(0, m/2), theta(n/4, m/2), theta(n/2, m/2), theta(3*n/4, m/2), theta(n, m/2)
			1 format(i6, 5f15.8)
		end if
	
	end do
	
	do j = 0, m
		do i = 0, n
			s = 0.0d0
			do k = 0, 8
				s = s + f(k, i, j)
			end do
			theta(i, j) = s
		end do
	end do
	
	open(17, file = "adiabatic_boundary_temp_profile.txt")
	open(27, file = "interface_boundary_temp_profile.txt")
	
	do i = 0, n
		write(17, *) i / float(n), theta(i, m)
	end do
	
	do j = 0, m
		write(27, *) theta(n/2, j), j / float(n)
	end do
	
end program heated_plate
