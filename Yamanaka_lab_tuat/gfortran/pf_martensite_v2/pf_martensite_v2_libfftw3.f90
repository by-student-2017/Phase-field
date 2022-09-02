!
! Purpose   :: Phase-field model for martensitic transformation in 2D
! Reference :: Y. U. Wang and A. G. Khachaturyan, Acta Materialia, 45 (1999), p. 759.
!		  * Please refere the following textbook for equations numbered in this source code
!		  T. Takaki and A. Yamanaka, Phase-field method, Yokendo-Ltd., (2012), Sections 7 & 8 (in Japanese)
! Programmer:: Akinori Yamanaka (Tokyo Univ. Agri. Tech.)
! Date	:: Jan., 2019
!
! << Please note that programmer does not take any responsibility or liability for any damage or loss caused by using this program >>
!

!================================
! Main program
!================================
program martensite
	implicit none

!================================
! libfftw3, gfortran *.f90 -I/usr/include -lfftw3
!================================
	include 'fftw3.f'
	integer(8)::  plan
	integer(8):: iplan
	complex(8), allocatable:: in(:,:), out(:,:)

!================================
! Setting of model parameters
!================================
	integer:: ig
	integer:: grid				!! number of computational grids (must be 2^ig)
	integer:: ndm
	integer:: half_grid
	integer:: step_end			!! number of time steps
	integer:: step_out			!! interval of output
	double precision:: length	!! length of computationa domain [m]
	double precision:: dx		!! spacing of computational grids
	double precision:: dt		!! time increment [dimensionless time]
	double precision:: temp		!! temperature [K]
	double precision, parameter:: pi = 3.141592
	double precision, parameter:: rr = 8.3145	!! gas constant
	double precision:: mobility	!! mobility of phase-field
	double precision:: vm0	 	!! molar volume [m3/mol]
	double precision:: aa0	 	!! magnitude of driving force [dimensionless]
	double precision:: aa		!! constant A for chemical free energy
	double precision:: bb		!! B
	double precision:: cc		!! C
	double precision:: grad_energ_coeff	!! gradient energy coefficient [dimensionless]
	double precision:: stress_unit	!! stress is regularized by this number
	double precision:: c11		!! elastic constant: C11=140GPa
	double precision:: c22
	double precision:: c33
	double precision:: c44		!! C44=(C11-C12)/2 because of isotropic material
	double precision:: c55
	double precision:: c66
	double precision:: c12		!! C12=84GPa
	double precision:: c21
	double precision:: c13
	double precision:: c31
	double precision:: c23
	double precision:: c32
	double precision:: ram0		!! Lame's constant ¥lambda
	double precision:: mu0		!! Lame's constant ¥mu
	double precision:: nu0		!! Poisson's ratio
	double precision:: munu0
	double precision:: sig22_a	!! applied stress along y direction
!	double precision, parameter::  sig22_a=2.e+06*vm0/rr/temp
	double precision:: ep11_a	!! applied strain due to the applied stress
	double precision:: ep22_a

!==========================================================
! Declaration of arrays and variables for phase field model
!==========================================================
	double precision, allocatable:: p1(:,:), p2(:,:)	!! phase-field variables 1 and 2
	double precision:: chem_pot1,elast_pot1,grad_pot1	!! chemical, elastic, and gradient potential for p1
	double precision:: chem_pot2,elast_pot2,grad_pot2
	double precision:: dpdt1,dpdt2						!! right hand side of Allen-Chan equation for p1 and p2
	double precision, allocatable:: eigen11_r(:,:), eigen22_r(:,:)	!! eigen strain in real space
	double precision, allocatable:: eigen12_r(:,:), eigen21_r(:,:)
	double precision, allocatable:: eigen11_f_real(:,:), eigen11_f_imag(:,:)	!! real and imagenary part of eigen strain in fourier space
	double precision, allocatable:: eigen22_f_real(:,:), eigen22_f_imag(:,:)
	double precision, allocatable:: eigen12_f_real(:,:), eigen12_f_imag(:,:)
	double precision, allocatable:: eigen21_f_real(:,:), eigen21_f_imag(:,:)
	double precision, allocatable:: inhomo_strain11(:,:), inhomo_strain22(:,:)	!! inhomogeneous strains
	double precision, allocatable:: inhomo_strain12(:,:), inhomo_strain21(:,:)
	double precision, allocatable:: stress11(:,:), stress22(:,:)				!! Cauchy stress
	double precision, allocatable:: stress33(:,:), stress12(:,:)
	double precision, allocatable:: elastic_strain11(:,:), elastic_strain22(:,:)!! elastic strain
	double precision, allocatable:: elastic_strain12(:,:), elastic_strain21(:,:)
	double precision, dimension(4,4):: eigen0_1,eigen0_2						!! eigen strain for each variants
	double precision:: homo_strain11,homo_strain22,homo_strain12,homo_strain21	!! homogeneous strain
	double precision:: el_strain11,el_strain22
	double precision:: eigen11_p1,eigen22_p1,eigen11_p2,eigen22_p2
	double precision:: sigma11_r,sigma22_r,sigma11_i,sigma22_i,sigma12_r,sigma12_i
	double precision:: sum11,sum22,sum12,sum21
	double precision:: s1,s2,rand

!==========================================================
! Declaration of arrays and variables for Fourier transform
!==========================================================
	double precision, allocatable:: xi(:,:), xr(:,:)	!! imagenary and real part of variable to be transformed
	double precision:: kxx,kyy,kxy,nnn,k_mag			!! wave vectors in Fourier space

	integer:: i,j,k,l,ii,jj,kk,iii,jjj,m,n,ip,im,jp,jm
	integer:: nstep,iout
	real(8):: random

	real(8)::data(40)
	open(1,file="parameters.txt")
	read(1,'(11x,e10.3)') (data(i),i=1,31)
	close(1)
	ig       = int(data(1))
	grid     = int(data(2))
	step_end = int(data(3))
	step_out = int(data(4))
	length   = data(5)
	dt       = data(6)
	temp     = data(7)
	mobility = data(8)
	vm0      = data(9)
	aa       = data(10)
	bb       = data(11)
	cc       = data(12)
	stress_unit=(1.0e+11)*vm0/rr/temp
	c11      = data(13)*stress_unit
	c22      = data(14)*stress_unit
	c33      = data(15)*stress_unit
	c44      = data(16)*stress_unit
	c55      = data(17)*stress_unit
	c66      = data(18)*stress_unit
	c12      = data(19)*stress_unit
	c13      = data(20)*stress_unit
	c21      = data(21)*stress_unit
	c23      = data(22)*stress_unit
	c31      = data(23)*stress_unit
	c32      = data(24)*stress_unit
	ram0     = data(25)*stress_unit ! ram0=c12   !! Lame's constant ¥lambda
	mu0      = data(26)*stress_unit ! mu0=c44    !! Lame's constant ¥mu
	nu0      = data(27) ! nu0=ram0/2./(ram0+mu0) !! Poisson's ratio
	munu0    = data(28) ! munu0=2.*mu0*(1.-nu0)
	sig22_a  = data(29) ! applied stress along y direction
	ep11_a   = data(30) ! ep11_a=-c12/(c11-c12)/(c11+2.*c12)*sig22_a
	ep22_a   = data(31) ! ep22_a=(c11+c12)/(c11-c12)/(c11+2.*c12)*sig22_a

	ndm = grid-1
	half_grid = grid/2
	dx = length/grid
	aa0 = (3.5e+08*(405.-temp)/405.)*vm0/rr/temp
	!bb = 3.*aa+12.
	!cc = 2.*aa+12.
	grad_energ_coeff=(1.0e-15)/rr/temp/dx/dx
	!stress_unit=(1.0e+11)*vm0/rr/temp
	!nu0=ram0/2./(ram0+mu0)
	!munu0=2.*mu0*(1.-nu0)
	ep11_a=-c12/(c11-c12)/(c11+2.*c12)*sig22_a
	ep22_a=(c11+c12)/(c11-c12)/(c11+2.*c12)*sig22_a

	allocate (p1(0:grid,0:grid))
	allocate (p2(0:grid,0:grid))
	allocate (eigen11_r(0:grid,0:grid))
	allocate (eigen22_r(0:grid,0:grid))
	allocate (eigen12_r(0:grid,0:grid))
	allocate (eigen21_r(0:grid,0:grid))
	allocate (eigen11_f_real(0:grid,0:grid))
	allocate (eigen11_f_imag(0:grid,0:grid))
	allocate (eigen22_f_real(0:grid,0:grid))
	allocate (eigen22_f_imag(0:grid,0:grid))
	allocate (eigen12_f_real(0:grid,0:grid))
	allocate (eigen12_f_imag(0:grid,0:grid))
	allocate (eigen21_f_real(0:grid,0:grid))
	allocate (eigen21_f_imag(0:grid,0:grid))
	allocate (inhomo_strain11(0:grid,0:grid))
	allocate (inhomo_strain22(0:grid,0:grid))
	allocate (inhomo_strain12(0:grid,0:grid))
	allocate (inhomo_strain21(0:grid,0:grid))
	allocate (stress11(0:grid,0:grid))
	allocate (stress22(0:grid,0:grid))
	allocate (stress33(0:grid,0:grid))
	allocate (stress12(0:grid,0:grid))
	allocate (elastic_strain11(0:grid,0:grid))
	allocate (elastic_strain22(0:grid,0:grid))
	allocate (elastic_strain12(0:grid,0:grid))
	allocate (elastic_strain21(0:grid,0:grid))
	allocate (xi(0:grid,0:grid))
	allocate (xr(0:grid,0:grid))
	
! settings for FFT and IFFT
	allocate ( in(0:ndm,0:ndm)) !! ndm = grid-1
	allocate (out(0:ndm,0:ndm)) !! ndm = grid-1
	call dfftw_plan_dft_2d( plan,grid,grid,in,out,FFTW_FORWARD, FFTW_ESTIMATE) !! forward FFT (FFT)
	call dfftw_plan_dft_2d(iplan,grid,grid,in,out,FFTW_BACKWARD,FFTW_ESTIMATE) !! inverse FFT (IFFT)

! input eigen strain for each variant
	eigen0_1 = 0
	eigen0_2 = 0
	eigen0_1(1,1)=-0.1994 ! eigen strain for variant 1 (phase field 1)
	eigen0_1(2,2)= 0.1322
	eigen0_2(1,1)= 0.1322 ! eigen strain for variant 2 (phase field 2)
	eigen0_2(2,2)=-0.1994

! setting initial distribution of phase-field variable
	do i=0,ndm
		do j=0,ndm
			call random_number(random)
			p1(i,j)=random	 ! rand() generates random numbers ranging from 0 to 1
			p2(i,j)=1-p1(i,j)
		enddo
	enddo

! output of initial condition
	iout = 0
	call output(iout,grid,ndm,vm0,rr,temp,p1,p2,stress11,stress22,stress12,stress33)

! start iteration
	do nstep=1,step_end

! calculation of eigen strain 11 and Fourier transform
		do i=0,ndm
			do j=0,ndm
				xr(i,j)=eigen0_1(1,1)*p1(i,j)+eigen0_2(1,1)*p2(i,j)  !! see Eq.(8.8)
				xi(i,j)=0.
				eigen11_r(i,j)=xr(i,j)
				in(i,j)= cmplx(xr(i,j),xi(i,j))
			enddo
		enddo

		call dfftw_execute_dft(plan, in, out) !! forward FFT

		do i=0,ndm
			do j=0,ndm
				!eigen11_f_real(i,j)=xr(i,j)  ! real part of eigen strain in Fourier space
				!eigen11_f_imag(i,j)=xi(i,j)
				eigen11_f_real(i,j)= dble( out(i,j) )
				eigen11_f_imag(i,j)=aimag( out(i,j) )
			enddo
		enddo
		eigen11_f_real(0,0)=0.
		eigen11_f_imag(0,0)=0.

! calculation of eigen strain 22 and Fourier transform
		do i=0,ndm
			do j=0,ndm
				xr(i,j)=eigen0_1(2,2)*p1(i,j)+eigen0_2(2,2)*p2(i,j)
				xi(i,j)=0.
				eigen22_r(i,j)=xr(i,j)
				in(i,j)= cmplx(xr(i,j),xi(i,j))
			enddo
		enddo

		call dfftw_execute_dft(plan, in, out) !! forward FFT

		do i=0,ndm
			do j=0,ndm
				!eigen22_f_real(i,j)=xr(i,j)
				!eigen22_f_imag(i,j)=xi(i,j)
				eigen22_f_real(i,j)= dble( out(i,j) )
				eigen22_f_imag(i,j)=aimag( out(i,j) )
			enddo
		enddo
		eigen22_f_real(0,0)=0.
		eigen22_f_imag(0,0)=0.

! calculation of eigen strain 12 and Fourier transform
		do i=0,ndm
			do j=0,ndm
				xr(i,j)=eigen0_1(1,2)*p1(i,j)+eigen0_2(1,2)*p2(i,j)
				xi(i,j)=0.
				eigen12_r(i,j)=xr(i,j)
				in(i,j)= cmplx(xr(i,j),xi(i,j))
			enddo
		enddo

		call dfftw_execute_dft(plan, in, out) !! forward FFT

		do i=0,ndm
			do j=0,ndm
				!eigen12_f_real(i,j)=xr(i,j)
				!eigen12_f_imag(i,j)=xi(i,j)
				eigen12_f_real(i,j)= dble( out(i,j) )
				eigen12_f_imag(i,j)=aimag( out(i,j) )
			enddo
		enddo
		eigen12_f_real(0,0)=0.
		eigen12_f_imag(0,0)=0.
		eigen21_f_real = eigen12_f_real  !! strain tensor is symmetry
		eigen21_f_imag = eigen12_f_imag

! calculation of homogeneous strain (this program assumes free surface)
		sum11=0.
		sum22=0.
		sum12=0.
		sum21=0.
		do i=0,ndm
			do j=0,ndm
				sum11=sum11+eigen11_r(i,j)
				sum22=sum22+eigen22_r(i,j)
				sum12=sum12+eigen12_r(i,j)
				sum21=sum21+eigen21_r(i,j)
			enddo
		enddo
		homo_strain11=sum11/grid/grid  !! homogeneous strain is volume average of eigen strain in this program
		homo_strain22=sum22/grid/grid  !! see Eq. (7.18)
		homo_strain12=sum12/grid/grid
		homo_strain21=homo_strain12

! calculation of inhomogeneous strain 11
		do i=0,ndm
			if(i.le.half_grid-1) ii=i
			if(i.ge.half_grid  ) ii=i-grid
			do j=0,ndm
				if(j.le.half_grid-1) jj=j
				if(j.ge.half_grid  ) jj=j-grid
				k_mag=ii*ii+jj*jj
				nnn=dsqrt(k_mag)	   !! magnutude of wave vector |k|
				if(nnn.eq.0.) nnn=1.
				kxx=(ii/nnn)*(ii/nnn)
				kyy=(jj/nnn)*(jj/nnn)
				kxy=(ii/nnn)*(jj/nnn)
				!! calculation of eigen stress (see Eq.(7.30))
				sigma11_r=ram0*(eigen11_f_real(i,j)+eigen22_f_real(i,j))+2.*mu0*eigen11_f_real(i,j)
				sigma22_r=ram0*(eigen11_f_real(i,j)+eigen22_f_real(i,j))+2.*mu0*eigen22_f_real(i,j)
				sigma12_r= mu0*(eigen12_f_real(i,j)+eigen21_f_real(i,j))
				sigma11_i=ram0*(eigen11_f_imag(i,j)+eigen22_f_imag(i,j))+2.*mu0*eigen11_f_imag(i,j)
				sigma22_i=ram0*(eigen11_f_imag(i,j)+eigen22_f_imag(i,j))+2.*mu0*eigen22_f_imag(i,j)
				sigma12_i= mu0*(eigen12_f_imag(i,j)+eigen21_f_imag(i,j))
				!! real part of inhomogeneous strain in Fourier space
				xr(i,j)=kxx*(1./mu0-kxx/munu0)*sigma11_r-kxy*kxy/munu0*sigma22_r &
					   +kxy*(1./mu0-kxx/munu0)*sigma12_r-kxy*kxx/munu0*sigma12_r
				!! imagenary part of inhomogeneous strain in Fourier space
				xi(i,j)=kxx*(1./mu0-kxx/munu0)*sigma11_i-kxy*kxy/munu0*sigma22_i &
					   +kxy*(1./mu0-kxx/munu0)*sigma12_i-kxy*kxx/munu0*sigma12_i
				in(i,j)= cmplx(xr(i,j),xi(i,j))
			enddo
		enddo

		call dfftw_execute_dft(iplan, in, out) !! inverse FFT (IFFT)

		!inhomo_strain11 = xr    !! inhomogeneous strain in real space
		do i=0,ndm
			do j=0,ndm
				inhomo_strain11(i,j)=dble( out(i,j) )/dble(grid*grid)
			enddo
		enddo

! calculation of inhomogeneous strain 22
		do i=0,ndm
			if(i.le.half_grid-1) ii=i
			if(i.ge.half_grid  ) ii=i-grid
			do j=0,ndm
				if(j.le.half_grid-1) jj=j
				if(j.ge.half_grid  ) jj=j-grid
				k_mag=ii*ii+jj*jj
				nnn=dsqrt(k_mag)
				if(nnn.eq.0.) nnn=1.
				kxx=(ii/nnn)*(ii/nnn)
				kyy=(jj/nnn)*(jj/nnn)
				kxy=(ii/nnn)*(jj/nnn)
				sigma11_r=ram0*(eigen11_f_real(i,j)+eigen22_f_real(i,j))+2.*mu0*eigen11_f_real(i,j)
				sigma22_r=ram0*(eigen11_f_real(i,j)+eigen22_f_real(i,j))+2.*mu0*eigen22_f_real(i,j)
				sigma12_r= mu0*(eigen12_f_real(i,j)+eigen21_f_real(i,j))
				sigma11_i=ram0*(eigen11_f_imag(i,j)+eigen22_f_imag(i,j))+2.*mu0*eigen11_f_imag(i,j)
				sigma22_i=ram0*(eigen11_f_imag(i,j)+eigen22_f_imag(i,j))+2.*mu0*eigen22_f_imag(i,j)
				sigma12_i= mu0*(eigen12_f_imag(i,j)+eigen21_f_imag(i,j))
				xr(i,j)=-kxy*kxy/munu0*sigma11_r+kyy*(1./mu0-kyy/munu0)*sigma22_r &
						-kxy*kyy/munu0*sigma12_r+kxy*(1./mu0-kyy/munu0)*sigma12_r
				xi(i,j)=-kxy*kxy/munu0*sigma11_i+kyy*(1./mu0-kyy/munu0)*sigma22_i &
						-kxy*kyy/munu0*sigma12_i+kxy*(1./mu0-kyy/munu0)*sigma12_i
				in(i,j)= cmplx(xr(i,j),xi(i,j))
			enddo
		enddo

		call dfftw_execute_dft(iplan, in, out) !! inverse FFT (IFFT)

		!inhomo_strain22 =xr
		do i=0,ndm
			do j=0,ndm
				inhomo_strain22(i,j)=dble( out(i,j) )/dble(grid*grid)
			enddo
		enddo

! calculation of inhomogeneous strain 12
		do i=0,ndm
			if(i.le.half_grid-1) ii=i
			if(i.ge.half_grid  ) ii=i-grid
			do j=0,ndm
				if(j.le.half_grid-1) jj=j
				if(j.ge.half_grid  ) jj=j-grid
				k_mag=ii*ii+jj*jj
				nnn=dsqrt(k_mag)
				if(nnn.eq.0.) nnn=1.
				kxx=(ii/nnn)*(ii/nnn)
				kyy=(jj/nnn)*(jj/nnn)
				kxy=(ii/nnn)*(jj/nnn)
				sigma11_r=ram0*(eigen11_f_real(i,j)+eigen22_f_real(i,j))+2.*mu0*eigen11_f_real(i,j)
				sigma22_r=ram0*(eigen11_f_real(i,j)+eigen22_f_real(i,j))+2.*mu0*eigen22_f_real(i,j)
				sigma12_r= mu0*(eigen12_f_real(i,j)+eigen21_f_real(i,j))
				sigma11_i=ram0*(eigen11_f_imag(i,j)+eigen22_f_imag(i,j))+2.*mu0*eigen11_f_imag(i,j)
				sigma22_i=ram0*(eigen11_f_imag(i,j)+eigen22_f_imag(i,j))+2.*mu0*eigen22_f_imag(i,j)
				sigma12_i= mu0*(eigen12_f_imag(i,j)+eigen21_f_imag(i,j))
				xr(i,j)=kxy*(1./2./mu0-kxx/munu0)*sigma11_r+(1./2./mu0*kxy-kyy*kxy/munu0)*sigma22_r  &
					   +kyy*(1./2./mu0-kxx/munu0)*sigma12_r+(1./2./mu0*kxx-kxy*kxy/munu0)*sigma12_r
				xi(i,j)=kxy*(1./2./mu0-kxx/munu0)*sigma11_i+(1./2./mu0*kxy-kyy*kxy/munu0)*sigma22_i  &
					   +kyy*(1./2./mu0-kxx/munu0)*sigma12_i+(1./2./mu0*kxx-kxy*kxy/munu0)*sigma12_i
				in(i,j)= cmplx(xr(i,j),xi(i,j))
			enddo
		enddo

		call dfftw_execute_dft(iplan, in, out) !! inverse FFT (IFFT)

		!inhomo_strain12 = xr
		do i=0,ndm
			do j=0,ndm
				inhomo_strain12(i,j)=dble( out(i,j) )/dble(grid*grid)
			enddo
		enddo
		inhomo_strain21 = inhomo_strain12

! calculation of Cauchy stress
		do i=0,ndm
			do j=0,ndm
				stress11(i,j)=(ram0+2.*mu0)*(homo_strain11+inhomo_strain11(i,j)-eigen11_r(i,j)) &   !! Cauchy stress 11
									  +ram0*(homo_strain22+inhomo_strain22(i,j)-eigen22_r(i,j))
				stress22(i,j)=ram0*(homo_strain11+inhomo_strain11(i,j)-eigen11_r(i,j)) &
				    +(ram0+2.*mu0)*(homo_strain22+inhomo_strain22(i,j)-eigen22_r(i,j))
				stress33(i,j)=ram0*((homo_strain11+inhomo_strain11(i,j)-eigen11_r(i,j)) &
								   +(homo_strain22+inhomo_strain22(i,j)-eigen22_r(i,j)))
				stress12(i,j)=mu0*((homo_strain12+inhomo_strain12(i,j)-eigen12_r(i,j)) &
								  +(homo_strain21+inhomo_strain21(i,j)-eigen21_r(i,j)))

				elastic_strain11(i,j)=homo_strain11+inhomo_strain11(i,j)+ep11_a-eigen11_r(i,j)     !! elastic strain 11
				elastic_strain22(i,j)=homo_strain22+inhomo_strain22(i,j)+ep22_a-eigen22_r(i,j)
				elastic_strain12(i,j)=homo_strain12+inhomo_strain12(i,j)-eigen12_r(i,j)
				elastic_strain21(i,j)=elastic_strain12(i,j)
			enddo
		enddo

! calculation of potentials (derivatives each energy density w.r.t. phase field variable)
		do i=0,ndm
			do j=0,ndm

				ip=i+1
				im=i-1
				jp=j+1
				jm=j-1
				if(i.eq.ndm) ip=0	!! periodic boundary condition is assumed
				if(i.eq.0  ) im=ndm
				if(j.eq.ndm) jp=0
				if(j.eq.0  ) jm=ndm

				! potential for gradient energy
				grad_pot1=-grad_energ_coeff*(p1(ip,j)+p1(im,j)+p1(i,jp)+p1(i,jm)-4.*p1(i,j))
				grad_pot2=-grad_energ_coeff*(p2(ip,j)+p2(im,j)+p2(i,jp)+p2(i,jm)-4.*p2(i,j))

				! potential for chemical free energy
				s1=p1(i,j)
				s2=p2(i,j)
				chem_pot1=aa0*s1*(aa-bb*s1+cc*(s1*s1+s2*s2))
				chem_pot2=aa0*s2*(aa-bb*s2+cc*(s1*s1+s2*s2))

				! potential for elastic strain energy (see Eq.(7.39))
				eigen11_p1=eigen0_1(1,1)
				eigen22_p1=eigen0_1(2,2)
				eigen11_p2=eigen0_2(1,1)
				eigen22_p2=eigen0_2(2,2)
				el_strain11= -elastic_strain11(i,j)
				el_strain22= -elastic_strain22(i,j)
				elast_pot1=ram0*(el_strain11+el_strain22)*(eigen11_p1+eigen22_p1) &
						+2.*mu0*(el_strain11*eigen11_p1+el_strain22*eigen22_p1)
				elast_pot2=ram0*(el_strain11+el_strain22)*(eigen11_p2+eigen22_p2) &
						+2.*mu0*(el_strain11*eigen11_p2+el_strain22*eigen22_p2)

! calculation of right hand side of Allen-Cahn equation (see Eq.(8.11))
		 		dpdt1=-mobility*(chem_pot1+grad_pot1+elast_pot1)
		 		dpdt2=-mobility*(chem_pot2+grad_pot2+elast_pot2)

! time integration of phase field variables
	 			p1(i,j)=p1(i,j)+dpdt1*dt
	 			p2(i,j)=p2(i,j)+dpdt2*dt

! adjust phase-field variables within 0~1
		 		if(p1(i,j).ge.1.) p1(i,j)=1.
		 		if(p1(i,j).le.0.) p1(i,j)=0.
	 			if(p2(i,j).ge.1.) p2(i,j)=1.
	 			if(p2(i,j).le.0.) p2(i,j)=0.

			enddo
		enddo

! output result
		if(mod(nstep,step_out).eq.0)then
			iout = iout + 1
			write(*,*) 'output = ',iout
			call output(iout,grid,ndm,vm0,rr,temp,p1,p2,stress11,stress22,stress12,stress33)
		endif

	enddo

	deallocate (p1)
	deallocate (p2)
	deallocate (eigen11_r)
	deallocate (eigen22_r)
	deallocate (eigen12_r)
	deallocate (eigen21_r)
	deallocate (eigen11_f_real)
	deallocate (eigen11_f_imag)
	deallocate (eigen22_f_real)
	deallocate (eigen22_f_imag)
	deallocate (eigen12_f_real)
	deallocate (eigen12_f_imag)
	deallocate (eigen21_f_real)
	deallocate (eigen21_f_imag)
	deallocate (inhomo_strain11)
	deallocate (inhomo_strain22)
	deallocate (inhomo_strain12)
	deallocate (inhomo_strain21)
	deallocate (stress11)
	deallocate (stress22)
	deallocate (stress33)
	deallocate (stress12)
	deallocate (elastic_strain11)
	deallocate (elastic_strain22)
	deallocate (elastic_strain12)
	deallocate (elastic_strain21)
	deallocate (xi)
	deallocate (xr)

	call dfftw_destroy_plan( plan)
	call dfftw_destroy_plan(iplan)

end program martensite

!===================================
! file output of results
! VTK file can be visualized by ParaView software.
! ParaView can be downloaded at https://www.paraview.org/
!===================================
	subroutine output(iout,grid,ndm,vm0,rr,temp,p_1,p_2,stress11,stress22,stress12,stress33)
	implicit none

	character(len=20) filename
	integer:: grid,ndm
	double precision:: vm0,rr,temp
	double precision, dimension(0:grid,0:grid):: p_1,p_2,stress11,stress22,stress12,stress33
	integer:: iout,l,m

	write(filename,'(a,i3.3,a)') 'mt_result',iout,'.vtk'
	open(1,file=filename)
	write(1,'(a)') '# vtk DataFile Version 3.0'
	write(1,'(a)') 'output.vtk'
	write(1,'(a)') 'ASCII'
	write(1,'(a)') 'DATASET STRUCTURED_POINTS'
	write(1,'(a,3i5)') 'DIMENSIONS', grid,grid,1
	write(1,'(a,3f4.1)') 'ORIGIN ', 0.0, 0.0, 0.0
	write(1,'(a,3i2)') 'ASPECT_RATIO', 1, 1, 1
	write(1,'(a,1i11)') 'POINT_DATA', grid*grid*1
	write(1,'(a)') 'SCALARS phase_field_1 float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do m=0,ndm
		do l=0,ndm
			write(1,'(1f10.6)') p_1(l,m)
		enddo
	enddo
	write(1,'(a)') 'SCALARS phase_field_2 float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do m=0,ndm
		do l=0,ndm
			write(1,'(1f10.6)') p_2(l,m)
		enddo
	enddo
	write(1,'(a)') 'SCALARS stress_11_MPa float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do m=0,ndm
		do l=0,ndm
			write(1,'(1f18.8)') stress11(l,m)/(1.e6*vm0/rr/temp) ! unit of stress [MPa]
		enddo
	enddo
	write(1,'(a)') 'SCALARS stress_22_MPa float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do m=0,ndm
		do l=0,ndm
			write(1,'(1f18.8)') stress22(l,m)/(1.e6*vm0/rr/temp) ! unit of stress [MPa]
		enddo
	enddo
	write(1,'(a)') 'SCALARS stress_33_MPa float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do m=0,ndm
		do l=0,ndm
			write(1,'(1f18.8)') stress33(l,m)/(1.e6*vm0/rr/temp) ! unit of stress [MPa]
		enddo
	enddo
	write(1,'(a)') 'SCALARS stress_12_MPa float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do m=0,ndm
		do l=0,ndm
			write(1,'(1f18.8)') stress12(l,m)/(1.e6*vm0/rr/temp) ! unit of stress [MPa]
		enddo
	enddo
	close(1)

	return
	end
