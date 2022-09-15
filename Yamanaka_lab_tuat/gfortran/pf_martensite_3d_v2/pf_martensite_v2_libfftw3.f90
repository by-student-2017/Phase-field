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
	complex(8), allocatable:: in(:,:,:), out(:,:,:)

!================================
! Setting of model parameters
!================================
	integer:: grid				!! number of computational grids (must be 2^ig)
	integer:: ndm
	integer:: half_grid
	integer:: step_end			!! number of time steps
	integer:: step_out			!! interval of output
	double precision:: length	!! length of computationa domain (m)
	double precision:: dx		!! spacing of computational grids
	double precision:: dt		!! time increment (dimensionless time)
	double precision:: temp		!! temperature (K)
	double precision, parameter:: pi = 3.141592
	double precision, parameter:: rr = 8.3145	!! gas constant
	double precision:: mobility	!! mobility of phase-field
	double precision:: vm0	 	!! molar volume (m3/mol)
	double precision:: aa0,aa0c	!! magnitude of driving force (dimensionless)
	double precision:: aa		!! constant A for chemical free energy
	double precision:: bb		!! B
	double precision:: cc		!! C
	double precision:: grad_energ_coeff, grad0c	!! gradient energy coefficient (dimensionless)
	double precision:: stress_unit	!! stress is regularized by this number
	double precision:: c11		!! elastic constant: C11=140GPa
	double precision:: c22
	double precision:: c33
	double precision:: c44		!! C44=(C11-C12)/2 because of isotropic material
	double precision:: c55
	double precision:: c66
	double precision:: c12		!! C12=84GPa
	double precision:: c13
	double precision:: c23
!	double precision:: ram0		!! Lame's constant ¥lambda
!	double precision:: mu0		!! Lame's constant ¥mu
!	double precision:: nu0		!! Poisson's ratio
!	double precision:: munu0
	double precision:: sig22_a	!! applied stress along y direction
!	double precision, parameter::  sig22_a=2.e+06*vm0/rr/temp
!	double precision:: ep11_a, ep22_a, ep33_a	!! applied strain due to the applied stress
!	double precision:: ep12_a, ep13_a, ep23_a

!==========================================================
! Declaration of arrays and variables for phase field model
!==========================================================
	double precision, allocatable:: p1(:,:,:), p2(:,:,:)	!! phase-field variables 1 and 2
	double precision:: chem_pot1,elast_pot1,grad_pot1	!! chemical, elastic, and gradient potential for p1
	double precision:: chem_pot2,elast_pot2,grad_pot2
	double precision:: dpdt1,dpdt2						!! right hand side of Allen-Chan equation for p1 and p2
	!
	double precision, allocatable:: eigen11_r(:,:,:), eigen22_r(:,:,:), eigen33_r(:,:,:)	!! eigen strain in real space
	double precision, allocatable:: eigen12_r(:,:,:), eigen21_r(:,:,:)
	double precision, allocatable:: eigen13_r(:,:,:), eigen31_r(:,:,:)
	double precision, allocatable:: eigen23_r(:,:,:), eigen32_r(:,:,:)
	!
	double precision, allocatable:: eigen11_f_real(:,:,:), eigen11_f_imag(:,:,:)	!! real and imagenary part of eigen strain in fourier space
	double precision, allocatable:: eigen22_f_real(:,:,:), eigen22_f_imag(:,:,:)
	double precision, allocatable:: eigen33_f_real(:,:,:), eigen33_f_imag(:,:,:)
	double precision, allocatable:: eigen12_f_real(:,:,:), eigen12_f_imag(:,:,:)
	double precision, allocatable:: eigen21_f_real(:,:,:), eigen21_f_imag(:,:,:)
	double precision, allocatable:: eigen13_f_real(:,:,:), eigen13_f_imag(:,:,:)
	double precision, allocatable:: eigen31_f_real(:,:,:), eigen31_f_imag(:,:,:)
	double precision, allocatable:: eigen23_f_real(:,:,:), eigen23_f_imag(:,:,:)
	double precision, allocatable:: eigen32_f_real(:,:,:), eigen32_f_imag(:,:,:)
	!
	double precision, allocatable:: inhomo_strain11(:,:,:), inhomo_strain22(:,:,:), inhomo_strain33(:,:,:)	!! inhomogeneous strains
	double precision, allocatable:: inhomo_strain12(:,:,:), inhomo_strain21(:,:,:)
	double precision, allocatable:: inhomo_strain13(:,:,:), inhomo_strain31(:,:,:)
	double precision, allocatable:: inhomo_strain23(:,:,:), inhomo_strain32(:,:,:)
	!
	double precision, allocatable:: stress11(:,:,:), stress22(:,:,:), stress33(:,:,:)	!! Cauchy stress
	double precision, allocatable:: stress12(:,:,:), stress13(:,:,:), stress23(:,:,:)
	!
	double precision, allocatable:: elastic_strain11(:,:,:), elastic_strain22(:,:,:), elastic_strain33(:,:,:)!! elastic strain
	double precision, allocatable:: elastic_strain12(:,:,:), elastic_strain21(:,:,:)
	double precision, allocatable:: elastic_strain13(:,:,:), elastic_strain31(:,:,:)
	double precision, allocatable:: elastic_strain23(:,:,:), elastic_strain32(:,:,:)
	!
	double precision, dimension(3,3):: eigen0_1		!! eigen strain for variant 1 (phase field 1)
	double precision, dimension(3,3):: eigen0_2		!! eigen strain for variant 1 (phase field 2)
	!
	double precision:: homo_strain11,homo_strain22,homo_strain33	!! homogeneous strain
	double precision:: homo_strain12,homo_strain21
	double precision:: homo_strain13,homo_strain31
	double precision:: homo_strain23,homo_strain32
	!
	double precision:: el_strain11,el_strain22,el_strain33
	double precision:: el_strain12,el_strain13,el_strain23
	!
	double precision:: eigen11_p1,eigen22_p1,eigen33_p1
	double precision:: eigen11_p2,eigen22_p2,eigen33_p2
	!
	double precision, dimension(3,3):: sigma_r	! sigma11_r,sigma22_r,sigma12_r...
	double precision, dimension(3,3):: sigma_i	! sigma11_i,sigma22_i,sigma12_i...
	double precision, dimension(3,3,3,3):: cec	! C ijkl
	!
	double precision, allocatable:: sigma11_r(:,:,:), sigma22_r(:,:,:), sigma33_r(:,:,:)	!! eigen strain in real space
	double precision, allocatable:: sigma12_r(:,:,:), sigma13_r(:,:,:), sigma23_r(:,:,:)
	!
	double precision, allocatable:: sigma11_i(:,:,:), sigma22_i(:,:,:), sigma33_i(:,:,:)	!! eigen strain in imaginary space
	double precision, allocatable:: sigma12_i(:,:,:), sigma13_i(:,:,:), sigma23_i(:,:,:)
	!
	double precision:: sum11,sum22,sum33
	double precision:: sum12,sum21
	double precision:: sum13,sum31
	double precision:: sum23,sum32
	!
	double precision:: s1,s2
	double precision:: rand
	!
	double precision, dimension(3,3):: ep_a	!! applied strain due to the applied stress

!==========================================================
! Declaration of arrays and variables for Fourier transform
!==========================================================
	double precision, allocatable:: xi(:,:,:), xr(:,:,:)	!! imagenary and real part of variable to be transformed
	double precision:: kxx,kyy,kxy,nnn,k_mag			!! wave vectors in Fourier space

	integer:: i,j,k,l,ii,jj,kk,iii,jjj,m,n
	integer:: ip,im,jp,jm,kp,km
	integer:: nstep,iout
	real(8):: random

	real(8)::data(40)
	open(1,file="parameters.txt")
	read(1,'(11x,e10.3)') (data(i),i=1,26)
	close(1)
	grid     = int(data(1))
	step_end = int(data(2))
	step_out = int(data(3))
	length   = data(4)
	dt       = data(5)
	temp     = data(6)
	mobility = data(7)
	aa0c     = data(8)
	vm0      = data(9)
	grad0c   = data(10)
	aa       = data(11)
	bb       = data(12)
	cc       = data(13)
	stress_unit=(1.0e+11)*vm0/rr/temp
	c11      = data(14)*stress_unit
	c22      = data(15)*stress_unit
	c33      = data(16)*stress_unit
	c44      = data(17)*stress_unit
	c55      = data(18)*stress_unit
	c66      = data(19)*stress_unit
	c12      = data(20)*stress_unit
	c13      = data(21)*stress_unit
	c23      = data(22)*stress_unit
	sig22_a  = data(23) ! applied stress along y direction
	ep_a = 0.0
	ep_a(1,1)= data(24) ! ep11_a=-c12/(c11-c12)/(c11+2.*c12)*sig22_a
	ep_a(2,2)= data(25) ! ep22_a=(c11+c12)/(c11-c12)/(c11+2.*c12)*sig22_a
	ep_a(3,3)= data(26) ! ep33_a=-c12/(c11-c12)/(c11+2.*c12)*sig22_a
	!ep11_a   = data(24) ! ep11_a=-c12/(c11-c12)/(c11+2.*c12)*sig22_a
	!ep22_a   = data(25) ! ep22_a=(c11+c12)/(c11-c12)/(c11+2.*c12)*sig22_a
	!ep33_a   = data(26) ! ep33_a=-c12/(c11-c12)/(c11+2.*c12)*sig22_a

	!ig=int(log(dble(grid+0.5))/log(2.0))
	!write(*,*) 'ig = ',ig
	
	ndm = grid-1
	half_grid = grid/2
	dx = length/grid
	aa0 = (aa0c*(405.-temp)/405.)*vm0/rr/temp
	!bb = 3.*aa+12.
	!cc = 2.*aa+12.
	grad_energ_coeff=grad0c/rr/temp/dx/dx
	!stress_unit=(1.0e+11)*vm0/rr/temp
	!nu0=ram0/2./(ram0+mu0)
	!munu0=2.*mu0*(1.-nu0)
	!ep11_a=-c12/(c11-c12)/(c11+2.*c12)*sig22_a
	!ep22_a=(c11+c12)/(c11-c12)/(c11+2.*c12)*sig22_a
	!ep33_a=-c12/(c11-c12)/(c11+2.*c12)*sig22_a

	allocate (p1(0:grid,0:grid,0:grid))
	allocate (p2(0:grid,0:grid,0:grid))
	!
	allocate (eigen11_r(0:grid,0:grid,0:grid))
	allocate (eigen22_r(0:grid,0:grid,0:grid))
	allocate (eigen33_r(0:grid,0:grid,0:grid))
	allocate (eigen12_r(0:grid,0:grid,0:grid))
	allocate (eigen21_r(0:grid,0:grid,0:grid))
	allocate (eigen13_r(0:grid,0:grid,0:grid))
	allocate (eigen31_r(0:grid,0:grid,0:grid))
	allocate (eigen23_r(0:grid,0:grid,0:grid))
	allocate (eigen32_r(0:grid,0:grid,0:grid))
	!
	! data after forward FFT
	allocate (eigen11_f_real(0:grid,0:grid,0:grid))
	allocate (eigen11_f_imag(0:grid,0:grid,0:grid))
	allocate (eigen22_f_real(0:grid,0:grid,0:grid))
	allocate (eigen22_f_imag(0:grid,0:grid,0:grid))
	allocate (eigen33_f_real(0:grid,0:grid,0:grid))
	allocate (eigen33_f_imag(0:grid,0:grid,0:grid))
	allocate (eigen12_f_real(0:grid,0:grid,0:grid))
	allocate (eigen12_f_imag(0:grid,0:grid,0:grid))
	allocate (eigen21_f_real(0:grid,0:grid,0:grid))
	allocate (eigen21_f_imag(0:grid,0:grid,0:grid))
	allocate (eigen13_f_real(0:grid,0:grid,0:grid))
	allocate (eigen13_f_imag(0:grid,0:grid,0:grid))
	allocate (eigen31_f_real(0:grid,0:grid,0:grid))
	allocate (eigen31_f_imag(0:grid,0:grid,0:grid))
	allocate (eigen23_f_real(0:grid,0:grid,0:grid))
	allocate (eigen23_f_imag(0:grid,0:grid,0:grid))
	allocate (eigen32_f_real(0:grid,0:grid,0:grid))
	allocate (eigen32_f_imag(0:grid,0:grid,0:grid))
	!
	allocate (inhomo_strain11(0:grid,0:grid,0:grid))
	allocate (inhomo_strain22(0:grid,0:grid,0:grid))
	allocate (inhomo_strain33(0:grid,0:grid,0:grid))
	allocate (inhomo_strain12(0:grid,0:grid,0:grid))
	allocate (inhomo_strain21(0:grid,0:grid,0:grid))
	allocate (inhomo_strain13(0:grid,0:grid,0:grid))
	allocate (inhomo_strain31(0:grid,0:grid,0:grid))
	allocate (inhomo_strain23(0:grid,0:grid,0:grid))
	allocate (inhomo_strain32(0:grid,0:grid,0:grid))
	!
	allocate (stress11(0:grid,0:grid,0:grid))
	allocate (stress22(0:grid,0:grid,0:grid))
	allocate (stress33(0:grid,0:grid,0:grid))
	allocate (stress12(0:grid,0:grid,0:grid))
	allocate (stress13(0:grid,0:grid,0:grid))
	allocate (stress23(0:grid,0:grid,0:grid))
	!
	allocate (elastic_strain11(0:grid,0:grid,0:grid))
	allocate (elastic_strain22(0:grid,0:grid,0:grid))
	allocate (elastic_strain33(0:grid,0:grid,0:grid))
	allocate (elastic_strain12(0:grid,0:grid,0:grid))
	allocate (elastic_strain21(0:grid,0:grid,0:grid))
	allocate (elastic_strain13(0:grid,0:grid,0:grid))
	allocate (elastic_strain31(0:grid,0:grid,0:grid))
	allocate (elastic_strain23(0:grid,0:grid,0:grid))
	allocate (elastic_strain32(0:grid,0:grid,0:grid))
	!
	allocate (xi(0:grid,0:grid,0:grid))
	allocate (xr(0:grid,0:grid,0:grid))
	!
	allocate (sigma11_r(0:grid,0:grid,0:grid))
	allocate (sigma22_r(0:grid,0:grid,0:grid))
	allocate (sigma33_r(0:grid,0:grid,0:grid))
	allocate (sigma12_r(0:grid,0:grid,0:grid))
	allocate (sigma13_r(0:grid,0:grid,0:grid))
	allocate (sigma23_r(0:grid,0:grid,0:grid))
	!
	allocate (sigma11_i(0:grid,0:grid,0:grid))
	allocate (sigma22_i(0:grid,0:grid,0:grid))
	allocate (sigma33_i(0:grid,0:grid,0:grid))
	allocate (sigma12_i(0:grid,0:grid,0:grid))
	allocate (sigma13_i(0:grid,0:grid,0:grid))
	allocate (sigma23_i(0:grid,0:grid,0:grid))
	
! initialization
	sigma_r = 0.0
	sigma_i = 0.0
	
! settings for FFT and IFFT
	allocate ( in(0:ndm,0:ndm,0:ndm)) !! ndm = grid-1
	allocate (out(0:ndm,0:ndm,0:ndm)) !! ndm = grid-1
	call dfftw_plan_dft_3d( plan,grid,grid,grid,in,out,FFTW_FORWARD, FFTW_ESTIMATE) !! forward FFT (FFT)
	call dfftw_plan_dft_3d(iplan,grid,grid,grid,in,out,FFTW_BACKWARD,FFTW_ESTIMATE) !! inverse FFT (IFFT)

! cec(i,j,k,l) = C ijkl
	cec = 0 ! initialization
	!
	! Memo 1 (cec(i,j,k,l)=cec(i,j,l,k)=cec(j,i,k,l)=cec(k,l,i,j))
	! cec(1,1,1,1)=c11;  cec(1,1,2,2)=c12;  cec(1,1,3,3)=c13;  c1123=c1131=c1112=0, c2311=c3111=c1211=0
	! cec(2,2,1,1)=c21;  cec(2,2,2,2)=c22;  cec(2,2,3,3)=c23;  c2223=c2231=c2212=0, c2322=c3122=c1222=0
	! cec(3,3,1,1)=c31;  cec(3,3,2,2)=c32;  cec(3,3,3,3)=c33;  c3323=c2231=c2212=0, c2333=c3133=c1232=0
	! cec(2,3,2,3)=c44;  c2331=c2312=0
	! cec(3,1,3,1)=c55;  c3123=c3112=0
	! cec(1,2,1,2)=c66;  c1223=c1231=0
	!
	cec(1,1,1,1)=c11
	cec(2,2,2,2)=c22
	cec(3,3,3,3)=c33
	cec(2,3,2,3)=c44;  cec(2,3,3,2)=c44; cec(3,2,2,3)=c44; cec(3,2,3,2)=c44
	cec(1,3,1,3)=c55;  cec(1,3,3,1)=c55; cec(3,1,1,3)=c55; cec(3,1,3,1)=c55
	cec(1,2,1,2)=c66;  cec(1,2,2,1)=c66; cec(2,1,1,2)=c66; cec(2,1,2,1)=c66
	cec(1,1,2,2)=c12; cec(2,2,1,1)=c12
	cec(1,1,3,3)=c13; cec(3,3,1,1)=c13
	cec(2,2,3,3)=c23; cec(3,3,2,2)=c23

! input eigen strain for each variant (Lattice mismatch, epsilon00(i,j))
	! eigen strain for variant 1 (phase field 1)
	eigen0_1 = 0
	eigen0_1(1,1)=-0.1994
	eigen0_1(2,2)= 0.1322
	eigen0_1(3,3)= 0.1322
	!
	! eigen strain for variant 2 (phase field 2)
	eigen0_2 = 0
	eigen0_2(1,1)= 0.1322
	eigen0_2(2,2)=-0.1994
	eigen0_2(3,3)= 0.1322

! setting initial distribution of phase-field variable
	do i=0,ndm
		do j=0,ndm
			do k=0,ndm
				call random_number(random)
				p1(i,j,k)=random	 ! rand() generates random numbers ranging from 0 to 1
				p2(i,j,k)=1-p1(i,j,k)
			enddo
		enddo
	enddo

! output of initial condition
	iout = 0
	call output(iout,grid,ndm,vm0,rr,temp,p1,p2,stress11,stress22,stress33,stress12,stress13,stress23)

! start iteration
	do nstep=1,step_end

! calculation of eigen strain 11 and Fourier transform
! eigen strain 11: eigen0_1(1,1), eigen0_2(1,1)
! Eq.(8.8): epsilon_0 ij = epsilon_01 ij * phi_1 + epsilon_02 ij * phi_2 + ...
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					xr(i,j,k)=eigen0_1(1,1)*p1(i,j,k)+eigen0_2(1,1)*p2(i,j,k)  !! see Eq.(8.8)
					xi(i,j,k)=0.
					eigen11_r(i,j,k)=xr(i,j,k)
					in(i,j,k)= cmplx(xr(i,j,k),xi(i,j,k))
				enddo
			enddo
		enddo

		call dfftw_execute_dft(plan, in, out) !! forward FFT

		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					!eigen11_f_real(i,j,k)=xr(i,j,k)  ! real part of eigen strain in Fourier space
					!eigen11_f_imag(i,j,k)=xi(i,j,k)
					eigen11_f_real(i,j,k)= dble( out(i,j,k) )
					eigen11_f_imag(i,j,k)=aimag( out(i,j,k) )
				enddo
			enddo
		enddo
		eigen11_f_real(0,0,0)=0.0
		eigen11_f_imag(0,0,0)=0.0

! calculation of eigen strain 22 and Fourier transform
! eigen strain 22: eigen0_1(2,2), eigen0_2(2,2)
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					xr(i,j,k)=eigen0_1(2,2)*p1(i,j,k)+eigen0_2(2,2)*p2(i,j,k)  !! see Eq.(8.8)
					xi(i,j,k)=0.
					eigen22_r(i,j,k)=xr(i,j,k)
					in(i,j,k)= cmplx(xr(i,j,k),xi(i,j,k))
				enddo
			enddo
		enddo

		call dfftw_execute_dft(plan, in, out) !! forward FFT

		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					!eigen22_f_real(i,j,k)=xr(i,j,k)
					!eigen22_f_imag(i,j,k)=xi(i,j,k)
					eigen22_f_real(i,j,k)= dble( out(i,j,k) )
					eigen22_f_imag(i,j,k)=aimag( out(i,j,k) )
				enddo
			enddo
		enddo
		eigen22_f_real(0,0,0)=0.0
		eigen22_f_imag(0,0,0)=0.0

! calculation of eigen strain 33 and Fourier transform
! eigen strain 33: eigen0_1(3,3), eigen0_2(3,3)
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					xr(i,j,k)=eigen0_1(3,3)*p1(i,j,k)+eigen0_2(3,3)*p2(i,j,k)  !! see Eq.(8.8)
					xi(i,j,k)=0.
					eigen33_r(i,j,k)=xr(i,j,k)
					in(i,j,k)= cmplx(xr(i,j,k),xi(i,j,k))
				enddo
			enddo
		enddo

		call dfftw_execute_dft(plan, in, out) !! forward FFT

		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					!eigen33_f_real(i,j,k)=xr(i,j,k)
					!eigen33_f_imag(i,j,k)=xi(i,j,k)
					eigen33_f_real(i,j,k)= dble( out(i,j,k) )
					eigen33_f_imag(i,j,k)=aimag( out(i,j,k) )
				enddo
			enddo
		enddo
		eigen33_f_real(0,0,0)=0.0
		eigen33_f_imag(0,0,0)=0.0

! calculation of eigen strain 12 and Fourier transform
! eigen strain 12: eigen0_1(1,2), eigen0_2(1,2)
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					xr(i,j,k)=eigen0_1(1,2)*p1(i,j,k)+eigen0_2(1,2)*p2(i,j,k)  !! see Eq.(8.8)
					xi(i,j,k)=0.
					eigen12_r(i,j,k)=xr(i,j,k)
					in(i,j,k)= cmplx(xr(i,j,k),xi(i,j,k))
				enddo
			enddo
		enddo

		call dfftw_execute_dft(plan, in, out) !! forward FFT

		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					!eigen12_f_real(i,j,k)=xr(i,j,k)
					!eigen12_f_imag(i,j,k)=xi(i,j,k)
					eigen12_f_real(i,j,k)= dble( out(i,j,k) )
					eigen12_f_imag(i,j,k)=aimag( out(i,j,k) )
				enddo
			enddo
		enddo
		eigen12_f_real(0,0,0)=0.0
		eigen12_f_imag(0,0,0)=0.0
		eigen21_f_real = eigen12_f_real  !! strain tensor is symmetry
		eigen21_f_imag = eigen12_f_imag

! calculation of eigen strain 13 and Fourier transform
! eigen strain 13: eigen0_1(1,3), eigen0_2(1,3)
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					xr(i,j,k)=eigen0_1(1,3)*p1(i,j,k)+eigen0_2(1,3)*p2(i,j,k)  !! see Eq.(8.8)
					xi(i,j,k)=0.
					eigen13_r(i,j,k)=xr(i,j,k)
					in(i,j,k)= cmplx(xr(i,j,k),xi(i,j,k))
				enddo
			enddo
		enddo

		call dfftw_execute_dft(plan, in, out) !! forward FFT

		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					!eigen13_f_real(i,j,k)=xr(i,j,k)
					!eigen13_f_imag(i,j,k)=xi(i,j,k)
					eigen13_f_real(i,j,k)= dble( out(i,j,k) )
					eigen13_f_imag(i,j,k)=aimag( out(i,j,k) )
				enddo
			enddo
		enddo
		eigen13_f_real(0,0,0)=0.0
		eigen13_f_imag(0,0,0)=0.0
		eigen31_f_real = eigen13_f_real  !! strain tensor is symmetry
		eigen31_f_imag = eigen13_f_imag

! calculation of eigen strain 23 and Fourier transform
! eigen strain 23: eigen0_1(2,3), eigen0_2(2,3)
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					xr(i,j,k)=eigen0_1(2,3)*p1(i,j,k)+eigen0_2(2,3)*p2(i,j,k)  !! see Eq.(8.8)
					xi(i,j,k)=0.
					eigen23_r(i,j,k)=xr(i,j,k)
					in(i,j,k)= cmplx(xr(i,j,k),xi(i,j,k))
				enddo
			enddo
		enddo

		call dfftw_execute_dft(plan, in, out) !! forward FFT

		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					!eigen23_f_real(i,j,k)=xr(i,j,k)
					!eigen23_f_imag(i,j,k)=xi(i,j,k)
					eigen23_f_real(i,j,k)= dble( out(i,j,k) )
					eigen23_f_imag(i,j,k)=aimag( out(i,j,k) )
				enddo
			enddo
		enddo
		eigen23_f_real(0,0,0)=0.0
		eigen23_f_imag(0,0,0)=0.0
		eigen32_f_real = eigen13_f_real  !! strain tensor is symmetry
		eigen32_f_imag = eigen13_f_imag

! calculation of homogeneous strain (this program assumes free surface)
		sum11=0.0; sum22=0.0; sum33=0.0
		sum12=0.0; sum21=0.0
		sum13=0.0; sum31=0.0
		sum23=0.0; sum32=0.0
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					sum11=sum11+eigen11_r(i,j,k)
					sum22=sum22+eigen22_r(i,j,k)
					sum33=sum33+eigen33_r(i,j,k)
					!
					sum12=sum12+eigen12_r(i,j,k); sum21=sum21+eigen21_r(i,j,k)
					sum13=sum13+eigen13_r(i,j,k); sum31=sum31+eigen31_r(i,j,k)
					sum23=sum23+eigen23_r(i,j,k); sum32=sum32+eigen32_r(i,j,k)
				enddo
			enddo
		enddo
		homo_strain11=sum11/dble(grid*grid*grid)	!! homogeneous strain is volume average of eigen strain in this program
		homo_strain22=sum22/dble(grid*grid*grid)	!! see Eq. (7.18) = (1/V)* integral epsilon_0 dV
		homo_strain33=sum33/dble(grid*grid*grid)
		homo_strain12=sum12/dble(grid*grid*grid)
		homo_strain21=homo_strain12
		homo_strain13=sum13/dble(grid*grid*grid)
		homo_strain31=homo_strain13
		homo_strain23=sum23/dble(grid*grid*grid)
		homo_strain32=homo_strain23
	
! calculation of sigma
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					!
					sigma11_r(i,j,k)=cec(1,1,1,1)*eigen11_f_real(i,j,k)+cec(1,1,2,2)*eigen22_f_real(i,j,k)+cec(1,1,3,3)*eigen33_f_real(i,j,k)
					sigma22_r(i,j,k)=cec(2,2,1,1)*eigen11_f_real(i,j,k)+cec(2,2,2,2)*eigen22_f_real(i,j,k)+cec(2,2,3,3)*eigen33_f_real(i,j,k)
					sigma33_r(i,j,k)=cec(3,3,1,1)*eigen11_f_real(i,j,k)+cec(3,3,2,2)*eigen22_f_real(i,j,k)+cec(3,3,3,3)*eigen33_f_real(i,j,k)
					sigma12_r(i,j,k)=cec(1,2,1,2)*2.0*eigen12_f_real(i,j,k)
					sigma13_r(i,j,k)=cec(1,3,1,3)*2.0*eigen13_f_real(i,j,k)
					sigma23_r(i,j,k)=cec(2,3,2,3)*2.0*eigen23_f_real(i,j,k)
					!
					sigma11_i(i,j,k)=cec(1,1,1,1)*eigen11_f_imag(i,j,k)+cec(1,1,2,2)*eigen22_f_imag(i,j,k)+cec(1,1,3,3)*eigen33_f_imag(i,j,k)
					sigma22_i(i,j,k)=cec(2,2,1,1)*eigen11_f_imag(i,j,k)+cec(2,2,2,2)*eigen22_f_imag(i,j,k)+cec(2,2,3,3)*eigen33_f_imag(i,j,k)
					sigma33_i(i,j,k)=cec(3,3,1,1)*eigen11_f_imag(i,j,k)+cec(3,3,2,2)*eigen22_f_imag(i,j,k)+cec(3,3,3,3)*eigen33_f_imag(i,j,k)
					sigma12_i(i,j,k)=cec(1,2,1,2)*2.0*eigen12_f_imag(i,j,k)
					sigma13_i(i,j,k)=cec(1,3,1,3)*2.0*eigen13_f_imag(i,j,k)
					sigma23_i(i,j,k)=cec(2,3,2,3)*2.0*eigen23_f_imag(i,j,k)
					!
				enddo
			enddo
		enddo

! calculation of inhomogeneous strain 11
		iii=1; jjj=1 ! strain iii jjj = strain 11
		do i=0,ndm
			if(i.le.half_grid-1) ii=i			!! periodic boundary conditions after forward FFT (FFT)
			if(i.ge.half_grid  ) ii=i-grid		!! periodic boundary conditions after forward FFT (FFT)
			do j=0,ndm
				if(j.le.half_grid-1) jj=j		!! periodic boundary conditions after forward FFT (FFT)
				if(j.ge.half_grid  ) jj=j-grid	!! periodic boundary conditions after forward FFT (FFT)
				do k=0,ndm
					if(k.le.half_grid-1) kk=k		!! periodic boundary conditions after forward FFT (FFT)
					if(k.ge.half_grid  ) kk=k-grid	!! periodic boundary conditions after forward FFT (FFT)
					!
					!! calculation of eigen stress (see Eq.(7.30)) (other text 4-14)
					!! ram0=lambda=C12, mu0=mu=C44
					!! C11 C12 C13   0   0   0  = lambda+2mu lambda     lambda     0  0  0  -> sigma11
					!! C21 C22 C23   0   0   0  = lambda     lambda+2mu lambda     0  0  0  -> sigma22
					!! C31 C32 C33   0   0   0  = lambda     lambda     lambda+2mu 0  0  0  -> sigma33
					!!   0   0   0 C44   0   0  = 0          0          0          mu 0  0  -> sigma23
					!!   0   0   0   0 C55   0  = 0          0          0          0  mu 0  -> sigma31
					!!   0   0   0   0   0 C66  = 0          0          0          0  0  mu -> sigma12
					!! [sigma(i,j)]=[C(i,j,k,l)][epsilon(k,l)]^t
					!! [epsilon(k,l)] = [epsilon11 epsilon22 epsilon33 2*epsilon23 2*epsilon31 2*epsilon12]
							! old version 1
							!sigma_r(1,1)=ram0*(eigen11_f_real(i,j,k)+eigen22_f_real(i,j,k)+eigen33_f_real(i,j,k)) + 2.0*mu0*eigen11_f_real(i,j,k)
							!sigma_r(2,2)=ram0*(eigen11_f_real(i,j,k)+eigen22_f_real(i,j,k)+eigen33_f_real(i,j,k)) + 2.0*mu0*eigen22_f_real(i,j,k)
							!sigma_r(3,3)=ram0*(eigen11_f_real(i,j,k)+eigen22_f_real(i,j,k)+eigen33_f_real(i,j,k)) + 2.0*mu0*eigen33_f_real(i,j,k)
							!sigma_r(1,2)= mu0*(eigen12_f_real(i,j,k)+eigen21_f_real(i,j,k))
							!sigma_r(1,3)= mu0*(eigen13_f_real(i,j,k)+eigen31_f_real(i,j,k))
							!sigma_r(2,3)= mu0*(eigen23_f_real(i,j,k)+eigen32_f_real(i,j,k))
						! old version 2
						!sigma_r(1,1)=cec(1,1,1,1)*eigen11_f_real(i,j,k)+cec(1,1,2,2)*eigen22_f_real(i,j,k)+cec(1,1,3,3)*eigen33_f_real(i,j,k)
						!sigma_r(2,2)=cec(2,2,1,1)*eigen11_f_real(i,j,k)+cec(2,2,2,2)*eigen22_f_real(i,j,k)+cec(2,2,3,3)*eigen33_f_real(i,j,k)
						!sigma_r(3,3)=cec(3,3,1,1)*eigen11_f_real(i,j,k)+cec(3,3,2,2)*eigen22_f_real(i,j,k)+cec(3,3,3,3)*eigen33_f_real(i,j,k)
							!sigma_r(1,2)=cec(1,2,1,2)*eigen12_f_real(i,j,k)+cec(1,2,2,1)*eigen21_f_real(i,j,k)
							!sigma_r(1,3)=cec(1,3,1,3)*eigen13_f_real(i,j,k)+cec(1,3,3,1)*eigen31_f_real(i,j,k)
							!sigma_r(2,3)=cec(2,3,2,3)*eigen23_f_real(i,j,k)+cec(2,3,3,2)*eigen32_f_real(i,j,k)
						!sigma_r(1,2)=cec(1,2,1,2)*2.0*eigen12_f_real(i,j,k)
						!sigma_r(1,3)=cec(1,3,1,3)*2.0*eigen13_f_real(i,j,k)
						!sigma_r(2,3)=cec(2,3,2,3)*2.0*eigen23_f_real(i,j,k)
					sigma_r(1,1)=sigma11_r(i,j,k); sigma_r(2,2)=sigma22_r(i,j,k); sigma_r(3,3)=sigma33_r(i,j,k)
					sigma_r(1,2)=sigma12_r(i,j,k); sigma_r(1,3)=sigma13_r(i,j,k); sigma_r(2,3)=sigma23_r(i,j,k)
					!
							! old version 1
							!sigma_i(1,1)=ram0*(eigen11_f_imag(i,j,k)+eigen22_f_imag(i,j,k)+eigen33_f_imag(i,j,k)) + 2.0*mu0*eigen11_f_imag(i,j,k)
							!sigma_i(2,2)=ram0*(eigen11_f_imag(i,j,k)+eigen22_f_imag(i,j,k)+eigen33_f_imag(i,j,k)) + 2.0*mu0*eigen22_f_imag(i,j,k)
							!sigma_i(3,3)=ram0*(eigen11_f_imag(i,j,k)+eigen22_f_imag(i,j,k)+eigen33_f_imag(i,j,k)) + 2.0*mu0*eigen33_f_imag(i,j,k)
							!sigma_i(1,2)= mu0*(eigen12_f_imag(i,j,k)+eigen21_f_imag(i,j,k))
							!sigma_i(1,3)= mu0*(eigen13_f_imag(i,j,k)+eigen31_f_imag(i,j,k))
							!sigma_i(2,3)= mu0*(eigen23_f_imag(i,j,k)+eigen32_f_imag(i,j,k))
						!sigma_i(1,1)=cec(1,1,1,1)*eigen11_f_imag(i,j,k)+cec(1,1,2,2)*eigen22_f_imag(i,j,k)+cec(1,1,3,3)*eigen33_f_imag(i,j,k)
						!sigma_i(2,2)=cec(2,2,1,1)*eigen11_f_imag(i,j,k)+cec(2,2,2,2)*eigen22_f_imag(i,j,k)+cec(2,2,3,3)*eigen33_f_imag(i,j,k)
						!sigma_i(3,3)=cec(3,3,1,1)*eigen11_f_imag(i,j,k)+cec(3,3,2,2)*eigen22_f_imag(i,j,k)+cec(3,3,3,3)*eigen33_f_imag(i,j,k)
							!sigma_i(1,2)=cec(1,2,1,2)*eigen12_f_imag(i,j,k)+cec(1,2,2,1)*eigen21_f_imag(i,j,k)
							!sigma_i(1,3)=cec(1,3,1,3)*eigen13_f_imag(i,j,k)+cec(1,3,3,1)*eigen31_f_imag(i,j,k)
							!sigma_i(2,3)=cec(2,3,2,3)*eigen23_f_imag(i,j,k)+cec(2,3,3,2)*eigen32_f_imag(i,j,k)
						!sigma_i(1,2)=cec(1,2,1,2)*2.0*eigen12_f_imag(i,j,k)
						!sigma_i(1,3)=cec(1,3,1,3)*2.0*eigen13_f_imag(i,j,k)
						!sigma_i(2,3)=cec(2,3,2,3)*2.0*eigen23_f_imag(i,j,k)
					sigma_i(1,1)=sigma11_i(i,j,k); sigma_i(2,2)=sigma22_i(i,j,k); sigma_i(3,3)=sigma33_i(i,j,k)
					sigma_i(1,2)=sigma12_i(i,j,k); sigma_i(1,3)=sigma13_i(i,j,k); sigma_i(2,3)=sigma23_i(i,j,k)
					!
					! lambda = v*E/((1+v)*(1-2*v)), mu=E/(2*(1+v))
					!! munu0=2.0*mu0*(1.0-nu0)
					!! nu0=ram0/2.0/(ram0+mu0) Poisson's ratio
					!
					k_mag=ii*ii+jj*jj+kk*kk
					nnn=dsqrt(k_mag)	   !! magnutude of wave vector |k|
					if(nnn.eq.0.0) nnn=1.0
					!
					!! real part of inhomogeneous strain in Fourier space
					xr(i,j,k)=zcij(ii, jj, kk, nnn, iii, jjj, sigma_r, cec)
					!! imaginary part of inhomogeneous strain in Fourier space
					xi(i,j,k)=zcij(ii, jj, kk, nnn, iii, jjj, sigma_i, cec)
					!
					in(i,j,k)= cmplx(xr(i,j,k),xi(i,j,k))
				enddo
			enddo
		enddo

		call dfftw_execute_dft(iplan, in, out) !! inverse FFT (IFFT)

		!inhomo_strain11 = xr    !! inhomogeneous strain in real space
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					inhomo_strain11(i,j,k)=dble( out(i,j,k) )/dble(grid*grid*grid)
				enddo
			enddo
		enddo

! calculation of inhomogeneous strain 22
		iii=2; jjj=2 ! strain iii jjj = strain 22
		do i=0,ndm
			if(i.le.half_grid-1) ii=i			!! periodic boundary conditions after forward FFT (FFT)
			if(i.ge.half_grid  ) ii=i-grid		!! periodic boundary conditions after forward FFT (FFT)
			do j=0,ndm
				if(j.le.half_grid-1) jj=j		!! periodic boundary conditions after forward FFT (FFT)
				if(j.ge.half_grid  ) jj=j-grid	!! periodic boundary conditions after forward FFT (FFT)
				do k=0,ndm
					if(k.le.half_grid-1) kk=k		!! periodic boundary conditions after forward FFT (FFT)
					if(k.ge.half_grid  ) kk=k-grid	!! periodic boundary conditions after forward FFT (FFT)
					!
					sigma_r(1,1)=sigma11_r(i,j,k); sigma_r(2,2)=sigma22_r(i,j,k); sigma_r(3,3)=sigma33_r(i,j,k)
					sigma_r(1,2)=sigma12_r(i,j,k); sigma_r(1,3)=sigma13_r(i,j,k); sigma_r(2,3)=sigma23_r(i,j,k)
					!
					sigma_i(1,1)=sigma11_i(i,j,k); sigma_i(2,2)=sigma22_i(i,j,k); sigma_i(3,3)=sigma33_i(i,j,k)
					sigma_i(1,2)=sigma12_i(i,j,k); sigma_i(1,3)=sigma13_i(i,j,k); sigma_i(2,3)=sigma23_i(i,j,k)
					!
					k_mag=ii*ii+jj*jj+kk*kk
					nnn=dsqrt(k_mag)	   !! magnutude of wave vector |k|
					if(nnn.eq.0.0) nnn=1.0
					!
					!! real part of inhomogeneous strain in Fourier space
					xr(i,j,k)=zcij(ii, jj, kk, nnn, iii, jjj, sigma_r, cec)
					!! imaginary part of inhomogeneous strain in Fourier space
					xi(i,j,k)=zcij(ii, jj, kk, nnn, iii, jjj, sigma_i, cec)
					!
					in(i,j,k)= cmplx(xr(i,j,k),xi(i,j,k))
				enddo
			enddo
		enddo

		call dfftw_execute_dft(iplan, in, out) !! inverse FFT (IFFT)

		!inhomo_strain22 =xr
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					inhomo_strain22(i,j,k)=dble( out(i,j,k) )/dble(grid*grid*grid)
				enddo
			enddo
		enddo

! calculation of inhomogeneous strain 33
		iii=3; jjj=3 ! strain iii jjj = strain 33
		do i=0,ndm
			if(i.le.half_grid-1) ii=i			!! periodic boundary conditions after forward FFT (FFT)
			if(i.ge.half_grid  ) ii=i-grid		!! periodic boundary conditions after forward FFT (FFT)
			do j=0,ndm
				if(j.le.half_grid-1) jj=j		!! periodic boundary conditions after forward FFT (FFT)
				if(j.ge.half_grid  ) jj=j-grid	!! periodic boundary conditions after forward FFT (FFT)
				do k=0,ndm
					if(k.le.half_grid-1) kk=k		!! periodic boundary conditions after forward FFT (FFT)
					if(k.ge.half_grid  ) kk=k-grid	!! periodic boundary conditions after forward FFT (FFT)
					!
					sigma_r(1,1)=sigma11_r(i,j,k); sigma_r(2,2)=sigma22_r(i,j,k); sigma_r(3,3)=sigma33_r(i,j,k)
					sigma_r(1,2)=sigma12_r(i,j,k); sigma_r(1,3)=sigma13_r(i,j,k); sigma_r(2,3)=sigma23_r(i,j,k)
					!
					sigma_i(1,1)=sigma11_i(i,j,k); sigma_i(2,2)=sigma22_i(i,j,k); sigma_i(3,3)=sigma33_i(i,j,k)
					sigma_i(1,2)=sigma12_i(i,j,k); sigma_i(1,3)=sigma13_i(i,j,k); sigma_i(2,3)=sigma23_i(i,j,k)
					!
					k_mag=ii*ii+jj*jj+kk*kk
					nnn=dsqrt(k_mag)	   !! magnutude of wave vector |k|
					if(nnn.eq.0.0) nnn=1.0
					!
					!! real part of inhomogeneous strain in Fourier space
					xr(i,j,k)=zcij(ii, jj, kk, nnn, iii, jjj, sigma_r, cec)
					!! imaginary part of inhomogeneous strain in Fourier space
					xi(i,j,k)=zcij(ii, jj, kk, nnn, iii, jjj, sigma_i, cec)
					!
					in(i,j,k)= cmplx(xr(i,j,k),xi(i,j,k))
				enddo
			enddo
		enddo

		call dfftw_execute_dft(iplan, in, out) !! inverse FFT (IFFT)

		!inhomo_strain33 =xr
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					inhomo_strain33(i,j,k)=dble( out(i,j,k) )/dble(grid*grid*grid)
				enddo
			enddo
		enddo

! calculation of inhomogeneous strain 12
		iii=1; jjj=2 ! strain iii jjj = strain 12
		do i=0,ndm
			if(i.le.half_grid-1) ii=i
			if(i.ge.half_grid  ) ii=i-grid
			do j=0,ndm
				if(j.le.half_grid-1) jj=j
				if(j.ge.half_grid  ) jj=j-grid
				do k=0,ndm
					if(k.le.half_grid-1) kk=k		!! periodic boundary conditions after forward FFT (FFT)
					if(k.ge.half_grid  ) kk=k-grid	!! periodic boundary conditions after forward FFT (FFT)
					!
					sigma_r(1,1)=sigma11_r(i,j,k); sigma_r(2,2)=sigma22_r(i,j,k); sigma_r(3,3)=sigma33_r(i,j,k)
					sigma_r(1,2)=sigma12_r(i,j,k); sigma_r(1,3)=sigma13_r(i,j,k); sigma_r(2,3)=sigma23_r(i,j,k)
					!
					sigma_i(1,1)=sigma11_i(i,j,k); sigma_i(2,2)=sigma22_i(i,j,k); sigma_i(3,3)=sigma33_i(i,j,k)
					sigma_i(1,2)=sigma12_i(i,j,k); sigma_i(1,3)=sigma13_i(i,j,k); sigma_i(2,3)=sigma23_i(i,j,k)
					!
					k_mag=ii*ii+jj*jj+kk*kk
					nnn=dsqrt(k_mag)	   !! magnutude of wave vector |k|
					if(nnn.eq.0.0) nnn=1.0
					!
					!! real part of inhomogeneous strain in Fourier space
					xr(i,j,k)=zcij(ii, jj, kk, nnn, iii, jjj, sigma_r, cec)
					!! imaginary part of inhomogeneous strain in Fourier space
					xi(i,j,k)=zcij(ii, jj, kk, nnn, iii, jjj, sigma_i, cec)
					!
					in(i,j,k)= cmplx(xr(i,j,k),xi(i,j,k))
				enddo
			enddo
		enddo

		call dfftw_execute_dft(iplan, in, out) !! inverse FFT (IFFT)

		!inhomo_strain12 = xr
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					inhomo_strain12(i,j,k)=dble( out(i,j,k) )/dble(grid*grid*grid)
				enddo
			enddo
		enddo
		inhomo_strain21 = inhomo_strain12

! calculation of inhomogeneous strain 13
		iii=1; jjj=3 ! strain iii jjj = strain 13
		do i=0,ndm
			if(i.le.half_grid-1) ii=i
			if(i.ge.half_grid  ) ii=i-grid
			do j=0,ndm
				if(j.le.half_grid-1) jj=j
				if(j.ge.half_grid  ) jj=j-grid
				do k=0,ndm
					if(k.le.half_grid-1) kk=k		!! periodic boundary conditions after forward FFT (FFT)
					if(k.ge.half_grid  ) kk=k-grid	!! periodic boundary conditions after forward FFT (FFT)
					!
					sigma_r(1,1)=sigma11_r(i,j,k); sigma_r(2,2)=sigma22_r(i,j,k); sigma_r(3,3)=sigma33_r(i,j,k)
					sigma_r(1,2)=sigma12_r(i,j,k); sigma_r(1,3)=sigma13_r(i,j,k); sigma_r(2,3)=sigma23_r(i,j,k)
					!
					sigma_i(1,1)=sigma11_i(i,j,k); sigma_i(2,2)=sigma22_i(i,j,k); sigma_i(3,3)=sigma33_i(i,j,k)
					sigma_i(1,2)=sigma12_i(i,j,k); sigma_i(1,3)=sigma13_i(i,j,k); sigma_i(2,3)=sigma23_i(i,j,k)
					!
					k_mag=ii*ii+jj*jj+kk*kk
					nnn=dsqrt(k_mag)	   !! magnutude of wave vector |k|
					if(nnn.eq.0.0) nnn=1.0
					!
					!! real part of inhomogeneous strain in Fourier space
					xr(i,j,k)=zcij(ii, jj, kk, nnn, iii, jjj, sigma_r, cec)
					!! imaginary part of inhomogeneous strain in Fourier space
					xi(i,j,k)=zcij(ii, jj, kk, nnn, iii, jjj, sigma_i, cec)
					!
					in(i,j,k)= cmplx(xr(i,j,k),xi(i,j,k))
				enddo
			enddo
		enddo

		call dfftw_execute_dft(iplan, in, out) !! inverse FFT (IFFT)

		!inhomo_strain13 = xr
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					inhomo_strain13(i,j,k)=dble( out(i,j,k) )/dble(grid*grid*grid)
				enddo
			enddo
		enddo
		inhomo_strain31 = inhomo_strain13


! calculation of inhomogeneous strain 23
		iii=2; jjj=3 ! strain iii jjj = strain 23
		do i=0,ndm
			if(i.le.half_grid-1) ii=i
			if(i.ge.half_grid  ) ii=i-grid
			do j=0,ndm
				if(j.le.half_grid-1) jj=j
				if(j.ge.half_grid  ) jj=j-grid
				do k=0,ndm
					if(k.le.half_grid-1) kk=k		!! periodic boundary conditions after forward FFT (FFT)
					if(k.ge.half_grid  ) kk=k-grid	!! periodic boundary conditions after forward FFT (FFT)
					!
					sigma_r(1,1)=sigma11_r(i,j,k); sigma_r(2,2)=sigma22_r(i,j,k); sigma_r(3,3)=sigma33_r(i,j,k)
					sigma_r(1,2)=sigma12_r(i,j,k); sigma_r(1,3)=sigma13_r(i,j,k); sigma_r(2,3)=sigma23_r(i,j,k)
					!
					sigma_i(1,1)=sigma11_i(i,j,k); sigma_i(2,2)=sigma22_i(i,j,k); sigma_i(3,3)=sigma33_i(i,j,k)
					sigma_i(1,2)=sigma12_i(i,j,k); sigma_i(1,3)=sigma13_i(i,j,k); sigma_i(2,3)=sigma23_i(i,j,k)
					!
					k_mag=ii*ii+jj*jj+kk*kk
					nnn=dsqrt(k_mag)	   !! magnutude of wave vector |k|
					if(nnn.eq.0.0) nnn=1.0
					!
					!! real part of inhomogeneous strain in Fourier space
					xr(i,j,k)=zcij(ii, jj, kk, nnn, iii, jjj, sigma_r, cec)
					!! imaginary part of inhomogeneous strain in Fourier space
					xi(i,j,k)=zcij(ii, jj, kk, nnn, iii, jjj, sigma_i, cec)
					!
					in(i,j,k)= cmplx(xr(i,j,k),xi(i,j,k))
				enddo
			enddo
		enddo

		call dfftw_execute_dft(iplan, in, out) !! inverse FFT (IFFT)

		!inhomo_strain23 = xr
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					inhomo_strain23(i,j,k)=dble( out(i,j,k) )/dble(grid*grid*grid)
				enddo
			enddo
		enddo
		inhomo_strain32 = inhomo_strain23

! calculation of Cauchy stress
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					!
						! old version
						!stress11(i,j,k)=(ram0+2.0*mu0)*(homo_strain11+inhomo_strain11(i,j,k)-eigen11_r(i,j,k)) &   !! Cauchy stress 11
						!			+ram0*(homo_strain22+inhomo_strain22(i,j,k)-eigen22_r(i,j,k)) &
						!			+ram0*(homo_strain33+inhomo_strain33(i,j,k)-eigen33_r(i,j,k))
						!stress22(i,j,k)=ram0*(homo_strain11+inhomo_strain11(i,j,k)-eigen11_r(i,j,k)) &
						!			+(ram0+2.0*mu0)*(homo_strain22+inhomo_strain22(i,j,k)-eigen22_r(i,j,k)) &
						!			+ram0*(homo_strain33+inhomo_strain33(i,j,k)-eigen33_r(i,j,k))
						!stress33(i,j,k)=ram0*(homo_strain11+inhomo_strain11(i,j,k)-eigen11_r(i,j,k)) &
						!			+(homo_strain22+inhomo_strain22(i,j,k)-eigen22_r(i,j,k)) &
						!			+(ram0+2.0*mu0)*(homo_strain33+inhomo_strain33(i,j,k)-eigen33_r(i,j,k))
						!stress12(i,j,k)=mu0*(homo_strain12+inhomo_strain12(i,j,k)-eigen12_r(i,j,k)) &
						!			+(homo_strain21+inhomo_strain21(i,j,k)-eigen21_r(i,j,k))
						!stress13(i,j,k)=mu0*(homo_strain13+inhomo_strain13(i,j,k)-eigen13_r(i,j,k)) &
						!			+(homo_strain31+inhomo_strain31(i,j,k)-eigen31_r(i,j,k))
						!stress23(i,j,k)=mu0*(homo_strain23+inhomo_strain23(i,j,k)-eigen23_r(i,j,k)) &
						!			+(homo_strain32+inhomo_strain32(i,j,k)-eigen32_r(i,j,k))
					!
					!! homogeneous strain is volume average of eigen strain in this program
					!! see Eq. (7.18) = (1/V)* integral epsilon_0 dV
					stress11(i,j,k)=cec(1,1,1,1)*( homo_strain11 + inhomo_strain11(i,j,k) - eigen11_r(i,j,k) ) &
								   +cec(1,1,2,2)*( homo_strain22 + inhomo_strain22(i,j,k) - eigen22_r(i,j,k) ) &
								   +cec(1,1,3,3)*( homo_strain33 + inhomo_strain33(i,j,k) - eigen33_r(i,j,k) )
					stress22(i,j,k)=cec(2,2,1,1)*( homo_strain11 + inhomo_strain11(i,j,k) - eigen11_r(i,j,k) ) &
								   +cec(2,2,2,2)*( homo_strain22 + inhomo_strain22(i,j,k) - eigen22_r(i,j,k) ) &
								   +cec(2,2,3,3)*( homo_strain33 + inhomo_strain33(i,j,k) - eigen33_r(i,j,k) )
					stress33(i,j,k)=cec(3,3,1,1)*( homo_strain11 + inhomo_strain11(i,j,k) - eigen11_r(i,j,k) ) &
								   +cec(3,3,2,2)*( homo_strain22 + inhomo_strain22(i,j,k) - eigen22_r(i,j,k) ) &
								   +cec(3,3,3,3)*( homo_strain33 + inhomo_strain33(i,j,k) - eigen33_r(i,j,k) )
					stress12(i,j,k)=cec(1,2,1,2)*( homo_strain12 + inhomo_strain12(i,j,k) - eigen12_r(i,j,k) ) &
								   +cec(1,2,2,1)*( homo_strain21 + inhomo_strain21(i,j,k) - eigen21_r(i,j,k) )
					stress13(i,j,k)=cec(1,3,1,3)*( homo_strain13 + inhomo_strain13(i,j,k) - eigen13_r(i,j,k) ) &
								   +cec(1,3,3,1)*( homo_strain31 + inhomo_strain31(i,j,k) - eigen31_r(i,j,k) )
					stress23(i,j,k)=cec(2,3,2,3)*( homo_strain23 + inhomo_strain23(i,j,k) - eigen23_r(i,j,k) ) &
								   +cec(2,3,3,2)*( homo_strain32 + inhomo_strain32(i,j,k) - eigen32_r(i,j,k) )
					!
					!! epsilonT ij = epsilon0 ij - average epsilon0 ij - d(epsilonC ij) - epsilonA ij
					elastic_strain11(i,j,k)=eigen11_r(i,j,k) - ( homo_strain11 + inhomo_strain11(i,j,k) + ep_a(1,1) )	!! elastic strain 11
					elastic_strain22(i,j,k)=eigen22_r(i,j,k) - ( homo_strain22 + inhomo_strain22(i,j,k) + ep_a(2,2) )	!! elastic strain 22
					elastic_strain33(i,j,k)=eigen33_r(i,j,k) - ( homo_strain33 + inhomo_strain33(i,j,k) + ep_a(3,3) )	!! elastic strain 33
					elastic_strain12(i,j,k)=eigen12_r(i,j,k) - ( homo_strain12 + inhomo_strain12(i,j,k) + ep_a(1,2) )
					elastic_strain21(i,j,k)=elastic_strain12(i,j,k)
					elastic_strain13(i,j,k)=eigen13_r(i,j,k) - ( homo_strain13 + inhomo_strain13(i,j,k) + ep_a(1,3) )
					elastic_strain31(i,j,k)=elastic_strain13(i,j,k)
					elastic_strain23(i,j,k)=eigen23_r(i,j,k) - ( homo_strain23 + inhomo_strain23(i,j,k) + ep_a(2,3) )
					elastic_strain32(i,j,k)=elastic_strain23(i,j,k)
					!
				enddo
			enddo
		enddo

! calculation of potentials (derivatives each energy density w.r.t. phase field variable)
		do i=0,ndm
			do j=0,ndm
				do k=0,ndm
					!
					ip=i+1; im=i-1
					jp=j+1; jm=j-1
					kp=k+1; km=k-1
					if(i.eq.0) im=ndm; if(i.eq.ndm) ip=0	!! periodic boundary condition is assumed
					if(j.eq.0) jm=ndm; if(j.eq.ndm) jp=0
					if(k.eq.0) km=ndm; if(k.eq.ndm) kp=0
					!
					! potential for gradient energy
					grad_pot1=-grad_energ_coeff*(p1(ip,j,k)+p1(im,j,k)+p1(i,jp,k)+p1(i,jm,k)+p1(i,j,kp)+p1(i,j,km)-6.*p1(i,j,k))
					grad_pot2=-grad_energ_coeff*(p2(ip,j,k)+p2(im,j,k)+p2(i,jp,k)+p2(i,jm,k)+p2(i,j,kp)+p2(i,j,km)-6.*p2(i,j,k))

					! potential for chemical free energy
					s1=p1(i,j,k)
					s2=p2(i,j,k)
					
					! Landau polynomial expansion
					! fchem = df*{(a/2)*sum(sp^2)-(b/3)*sum(sp^3)-(c/4)*(sum(sp^2))^2}
					! chem_pot = d(fchem)/d(sp)
					chem_pot1=aa0*s1*(aa-bb*s1+cc*(s1*s1+s2*s2))
					chem_pot2=aa0*s2*(aa-bb*s2+cc*(s1*s1+s2*s2))

					! potential for elastic strain energy (see Eq.(7.39))
						!eigen11_p1=eigen0_1(1,1)	! epsilon00 11 (p=p1)
						!eigen22_p1=eigen0_1(2,2)	! epsilon00 22 (p=p1)
						!eigen33_p1=eigen0_1(3,3)	! epsilon00 33 (p=p1)
						!eigen11_p2=eigen0_2(1,1)	! epsilon00 11 (p=p2)
						!eigen22_p2=eigen0_2(2,2)	! epsilon00 22 (p=p2)
						!eigen33_p2=eigen0_2(3,3)	! epsilon00 33 (p=p2)
					!
					el_strain11= elastic_strain11(i,j,k)	! -epsilonT ij = -(epsilon0 ij - average epsilon0 ij - d(epsilonC ij) - epsilonA ij)
					el_strain22= elastic_strain22(i,j,k)	! -epsilonT ij = -(epsilon0 ij - average epsilon0 ij - d(epsilonC ij) - epsilonA ij)
					el_strain33= elastic_strain33(i,j,k)	! -epsilonT ij = -(epsilon0 ij - average epsilon0 ij - d(epsilonC ij) - epsilonA ij)
					el_strain12= elastic_strain12(i,j,k)
					el_strain13= elastic_strain13(i,j,k)
					el_strain23= elastic_strain23(i,j,k)
					!
						!elast_pot1=c12*(el_strain11+el_strain22+el_strain33)*(eigen11_p1+eigen22_p1+eigen33_p1) &
						!		+2.0*c44*(el_strain11*eigen11_p1+el_strain22*eigen22_p1+el_strain33*eigen33_p1)
						!elast_pot2=c12*(el_strain11+el_strain22+el_strain33)*(eigen11_p2+eigen22_p2+eigen33_p2) &
						!		+2.0*c44*(el_strain11*eigen11_p2+el_strain22*eigen22_p2+el_strain33*eigen33_p2)
					!! d(Estr)/d(sp1) = C(i,j,k,l) * [el_strain(k,l)] * eigen(i,j)_p1
					elast_pot1=cec(1,1,1,1)*eigen0_1(1,1)*el_strain11 &
							  +cec(2,2,2,2)*eigen0_1(2,2)*el_strain22 &
							  +cec(3,3,3,3)*eigen0_1(3,3)*el_strain33 &
							  +cec(1,1,2,2)*eigen0_1(1,1)*el_strain22 + cec(2,2,1,1)*eigen0_1(2,2)*el_strain11 &
							  +cec(1,1,3,3)*eigen0_1(1,1)*el_strain33 + cec(3,3,1,1)*eigen0_1(3,3)*el_strain11 &
							  +cec(2,2,3,3)*eigen0_1(2,2)*el_strain33 + cec(3,3,2,2)*eigen0_1(3,3)*el_strain22 &
							  +cec(2,3,2,3)*eigen0_1(2,3)*el_strain23*4.0 &
							  +cec(1,3,1,3)*eigen0_1(1,3)*el_strain13*4.0 &
							  +cec(1,2,1,2)*eigen0_1(1,2)*el_strain12*4.0
					elast_pot2=cec(1,1,1,1)*eigen0_2(1,1)*el_strain11 &
							  +cec(2,2,2,2)*eigen0_2(2,2)*el_strain22 &
							  +cec(3,3,3,3)*eigen0_2(3,3)*el_strain33 &
							  +cec(1,1,2,2)*eigen0_2(1,1)*el_strain22 + cec(2,2,1,1)*eigen0_2(2,2)*el_strain11 &
							  +cec(1,1,3,3)*eigen0_2(1,1)*el_strain33 + cec(3,3,1,1)*eigen0_2(3,3)*el_strain11 &
							  +cec(2,2,3,3)*eigen0_2(2,2)*el_strain33 + cec(3,3,2,2)*eigen0_2(3,3)*el_strain22 &
							  +cec(2,3,2,3)*eigen0_2(2,3)*el_strain23*4.0 &
							  +cec(1,3,1,3)*eigen0_2(1,3)*el_strain13*4.0 &
							  +cec(1,2,1,2)*eigen0_2(1,2)*el_strain12*4.0

! calculation of right hand side of Allen-Cahn equation (see Eq.(8.11))
			 		dpdt1=-mobility*(chem_pot1+grad_pot1+elast_pot1)
			 		dpdt2=-mobility*(chem_pot2+grad_pot2+elast_pot2)

! time integration of phase field variables
	 				p1(i,j,k)=p1(i,j,k)+dpdt1*dt
	 				p2(i,j,k)=p2(i,j,k)+dpdt2*dt

! adjust phase-field variables within 0~1
		 			if(p1(i,j,k).ge.1.0) p1(i,j,k)=1.0; if(p1(i,j,k).le.0.0) p1(i,j,k)=0.0
		 			if(p2(i,j,k).ge.1.0) p2(i,j,k)=1.0; if(p2(i,j,k).le.0.0) p2(i,j,k)=0.0
		 			!
				enddo
			enddo
		enddo

! output result
		if(mod(nstep,step_out).eq.0)then
			iout = iout + 1
			write(*,*) 'output = ',iout
			call output(iout,grid,ndm,vm0,rr,temp,p1,p2,stress11,stress22,stress33,stress12,stress13,stress23)
		endif

	enddo

	deallocate (p1)
	deallocate (p2)
	!
	deallocate (eigen11_r)
	deallocate (eigen22_r)
	deallocate (eigen33_r)
	deallocate (eigen12_r)
	deallocate (eigen21_r)
	deallocate (eigen13_r)
	deallocate (eigen31_r)
	deallocate (eigen23_r)
	deallocate (eigen32_r)
	!
	deallocate (eigen11_f_real)
	deallocate (eigen11_f_imag)
	deallocate (eigen22_f_real)
	deallocate (eigen22_f_imag)
	deallocate (eigen33_f_real)
	deallocate (eigen33_f_imag)
	deallocate (eigen12_f_real)
	deallocate (eigen12_f_imag)
	deallocate (eigen21_f_real)
	deallocate (eigen21_f_imag)
	deallocate (eigen13_f_real)
	deallocate (eigen13_f_imag)
	deallocate (eigen31_f_real)
	deallocate (eigen31_f_imag)
	deallocate (eigen23_f_real)
	deallocate (eigen23_f_imag)
	deallocate (eigen32_f_real)
	deallocate (eigen32_f_imag)
	!
	deallocate (inhomo_strain11)
	deallocate (inhomo_strain22)
	deallocate (inhomo_strain33)
	deallocate (inhomo_strain12)
	deallocate (inhomo_strain21)
	deallocate (inhomo_strain13)
	deallocate (inhomo_strain31)
	deallocate (inhomo_strain23)
	deallocate (inhomo_strain32)
	!
	deallocate (stress11)
	deallocate (stress22)
	deallocate (stress33)
	deallocate (stress12)
	deallocate (stress13)
	deallocate (stress23)
	!
	deallocate (elastic_strain11)
	deallocate (elastic_strain22)
	deallocate (elastic_strain33)
	deallocate (elastic_strain12)
	deallocate (elastic_strain21)
	deallocate (elastic_strain13)
	deallocate (elastic_strain31)
	deallocate (elastic_strain23)
	deallocate (elastic_strain32)
	!
	deallocate (xi)
	deallocate (xr)
	!
	deallocate (sigma11_r)
	deallocate (sigma22_r)
	deallocate (sigma33_r)
	deallocate (sigma12_r)
	deallocate (sigma13_r)
	deallocate (sigma23_r)
	!
	deallocate (sigma11_i)
	deallocate (sigma22_i)
	deallocate (sigma33_i)
	deallocate (sigma12_i)
	deallocate (sigma13_i)
	deallocate (sigma23_i)

	call dfftw_destroy_plan( plan)
	call dfftw_destroy_plan(iplan)

	stop
	!
contains
!===================================
! Eq.(7.30)
! Omega ik = (C ijkl * nj * nl)^-1
!===================================
	real(8) function zcij(ii, jj, kk, nnn, iii, jjj, sigma, cec)
	implicit none
	
	integer:: ii, jj, kk
	double precision:: nnn
	integer:: iii, jjj
	double precision, dimension(3,3):: sigma
	double precision, dimension(3,3,3,3):: cec
	
	integer:: m, n
	double precision:: nx, ny, nz
	double precision, dimension(3):: nec
	!
	double precision:: a11, a22, a33, a12, a13, a23, a21, a31, a32
	double precision:: b11, b22, b33, b12, b13, b23, b21, b31, b32
	double precision:: det1, zij
	!
	double precision, dimension(3,3):: om
	
	nec(1)=ii/nnn; nx=nec(1)	! n is the unit vector in the k direction, n=k/|k|
	nec(2)=jj/nnn; ny=nec(2)	! n is the unit vector in the k direction, n=k/|k|
	nec(3)=kk/nnn; nz=nec(3)	! n is the unit vector in the k direction, n=k/|k|

	! C(i,k,j,l)*n(j)*n(l), C: elastic modulus, n: unit vector
	! C = cec
	a11=cec(1,1,1,1)*nx*nx+cec(1,2,1,2)*ny*ny+cec(1,3,1,3)*nz*nz
	a22=cec(1,2,1,2)*nx*nx+cec(2,2,2,2)*ny*ny+cec(2,3,2,3)*nz*nz
	a33=cec(3,1,3,1)*nx*nx+cec(2,3,2,3)*ny*ny+cec(3,3,3,3)*nz*nz
	a12=(cec(1,1,2,2)+cec(1,2,1,2))*nx*ny
	a23=(cec(2,2,3,3)+cec(2,3,2,3))*ny*nz
	a31=(cec(3,3,1,1)+cec(3,1,3,1))*nx*nz
	a21=a12
	a32=a23
	a13=a31

	! cofactor
	b11=a22*a33-a23*a32
	b22=a11*a33-a13*a31
	b33=a11*a22-a12*a21
	b12=-(a21*a33-a23*a31)
	b23=-(a11*a32-a12*a31)
	b31=-(a22*a31-a21*a32)
	b21=b12
	b32=b23
	b13=b31

	! det (C(i,k,j,l)*n(j)*n(l))
	det1=a11*a22*a33+a12*a23*a31+a13*a32*a21-a13*a31*a22-a11*a23*a32-a33*a12*a21;
	if(det1==0.0) det1=1.0

	! inverse matrix
	om(1,1)=b11/det1
	om(2,2)=b22/det1
	om(3,3)=b33/det1
	om(1,2)=b12/det1; om(2,1)=om(1,2)
	om(2,3)=b23/det1; om(3,2)=om(2,3)
	om(3,1)=b31/det1; om(1,3)=om(3,1)
	
	! sigma(i,j) = cec(i,j,k,l)*eta_c(k,l)*(Kronecker delta(k,l))
	! sigma: Eigen stress
	zij=0.0
	do m=1,3
		do n=1,3
			! /zij=zij+0.5*( om(m,iii)*nec(n)*nec(jjj) + om(m,jjj)*nec(n)*nec(iii) )*sigma(m,n); // eq.(5.26) or eq.(II 3.5)
			zij=zij+0.5*( nec(jjj)*om(m,iii) + nec(iii)*om(m,jjj) )*nec(n)*sigma(m,n) ! Eq.(7.30)
		enddo
	enddo
	
	zcij = zij ! return data
	
	return
	end function zcij
	!
end program martensite

!===================================
! file output of results
! VTK file can be visualized by ParaView software.
! ParaView can be downloaded at https://www.paraview.org/
!===================================
	subroutine output(iout,grid,ndm,vm0,rr,temp,p_1,p_2,stress11,stress22,stress33,stress12,stress13,stress23)
	implicit none

	character(len=20) filename
	integer:: grid,ndm
	double precision:: vm0,rr,temp
	double precision, dimension(0:grid,0:grid,0:grid):: p_1,p_2
	double precision, dimension(0:grid,0:grid,0:grid):: stress11,stress22,stress33
	double precision, dimension(0:grid,0:grid,0:grid):: stress12,stress13,stress23
	integer:: iout,i,j,k

	write(filename,'(a,i3.3,a)') 'mt_result',iout,'.vtk'
	open(1,file=filename)
	write(1,'(a)') '# vtk DataFile Version 3.0'
	write(1,'(a)') 'output.vtk'
	write(1,'(a)') 'ASCII'
	write(1,'(a)') 'DATASET STRUCTURED_POINTS'
	write(1,'(a,3i5)') 'DIMENSIONS', grid,grid,grid
	write(1,'(a,3f4.1)') 'ORIGIN ', 0.0, 0.0, 0.0
	write(1,'(a,3i2)') 'ASPECT_RATIO', 1, 1, 1
	write(1,'(a,1i11)') 'POINT_DATA', grid*grid*grid
	write(1,'(a)') 'SCALARS phase_field_1 float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do k=0,ndm
		do j=0,ndm
			do i=0,ndm
				write(1,'(1f10.6)') p_1(i,j,k)
			enddo
		enddo
	enddo
	!
	write(1,'(a)') 'SCALARS phase_field_2 float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do k=0,ndm
		do j=0,ndm
			do i=0,ndm
				write(1,'(1f10.6)') p_2(i,j,k)
			enddo
		enddo
	enddo
	!
	write(1,'(a)') 'SCALARS stress_11_MPa float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do k=0,ndm
		do j=0,ndm
			do i=0,ndm
				write(1,'(1f18.8)') stress11(i,j,k)/(1.e6*vm0/rr/temp) ! unit of stress [MPa]
			enddo
		enddo
	enddo
	!
	write(1,'(a)') 'SCALARS stress_22_MPa float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do k=0,ndm
		do j=0,ndm
			do i=0,ndm
				write(1,'(1f18.8)') stress22(i,j,k)/(1.e6*vm0/rr/temp) ! unit of stress [MPa]
			enddo
		enddo
	enddo
	!
	write(1,'(a)') 'SCALARS stress_33_MPa float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do k=0,ndm
		do j=0,ndm
			do i=0,ndm
				write(1,'(1f18.8)') stress33(i,j,k)/(1.e6*vm0/rr/temp) ! unit of stress [MPa]
			enddo
		enddo
	enddo
	!
	write(1,'(a)') 'SCALARS stress_12_MPa float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do k=0,ndm
		do j=0,ndm
			do i=0,ndm
				write(1,'(1f18.8)') stress12(i,j,k)/(1.e6*vm0/rr/temp) ! unit of stress [MPa]
			enddo
		enddo
	enddo
	!
	write(1,'(a)') 'SCALARS stress_13_MPa float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do k=0,ndm
		do j=0,ndm
			do i=0,ndm
				write(1,'(1f18.8)') stress13(i,j,k)/(1.e6*vm0/rr/temp) ! unit of stress [MPa]
			enddo
		enddo
	enddo
	!
	write(1,'(a)') 'SCALARS stress_23_MPa float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do k=0,ndm
		do j=0,ndm
			do i=0,ndm
				write(1,'(1f18.8)') stress23(i,j,k)/(1.e6*vm0/rr/temp) ! unit of stress [MPa]
			enddo
		enddo
	enddo
	close(1)

	return
	end subroutine output
