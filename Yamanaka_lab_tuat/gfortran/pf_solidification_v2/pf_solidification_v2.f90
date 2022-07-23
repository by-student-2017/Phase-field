!
! Purpose   :: Phase-field model for dendrite solification in binary-alloy based on Warren-Boettinger-McFadden(WBM) model
! Reference :: J.A. Warren and W. J. Boettinger, Acta Metall. Mater., 43 (1995), p. 689
!              * Please refere the following textbook for equations numbered in this source code
!              T. Takaki and A. Yamanaka, Phase-field method, Yokendo-Ltd., (2012), Section 5.3 (in Japanese)
! Programmer:: Akinori Yamanaka (Tokyo Univ. Agri. Tech.)
! Date      :: Jan., 2019
!
! << Please note that programmer does not take any responsibility or liability for any damage or loss caused by using this program >>
!

!================================
! Setting of model parameters
!================================
module variables
 implicit none
 !contains

 integer :: step_end   ! total number of time steps
 integer :: step_out   ! interval of output
 integer :: grid_x     ! number of computational grids along x-axis
 integer :: grid_y

 double precision :: dt  ! time increment [s]
 double precision :: dx  ! spacing of computational grids along x-axis [m]
 double precision :: dy  ! spacing of computational grids along y-axis [m]
 double precision, parameter :: pi = 3.141592
 double precision :: latent_heat_A      ! latent heat of element A [J/m3]
 double precision :: latent_heat_B
 double precision :: melting_temp_A     ! melting temperature of element A [K]
 double precision :: melting_temp_B
 double precision :: int_energy_A       ! interfacial energy of element A [J/m2]
 double precision :: int_energy_B
 double precision :: kinetic_coeff_A    ! kinetic coefficient [m/K/s]
 double precision :: kinetic_coeff_B
 double precision :: diff_coeff_liquid  ! diffusion coefficient of element B in liquid [m2/s]
 double precision :: diff_coeff_solid
 double precision :: vmol               ! molar volume [m3/mol]
 double precision, parameter :: Rgas = 8.31451   ! gas constant [J/K/mol]
 double precision :: aniso_mode         ! anisotropic mode number of interfacial anisotropy
 double precision :: aniso_stre         ! strength of interfacial anisotropy (must be smaller than 1/((aniso_mode)^2-1)
 double precision :: magnitude_noise    ! magnitude of noise
 double precision :: theta_0            ! preferential growth direction
 double precision :: initial_concentration   ! initial concentration of element B in liquid []
 double precision :: initial_temperature     ! initial temperature of the system [K]
 double precision :: thick_diff_interface    ! thickness of diffuse interface [m]
 double precision :: bbb
 double precision :: tinb
 double precision :: grad_energ_coeff        ! gradient energy coefficient
 double precision :: W_A                     ! double-well potential height
 double precision :: W_B
 double precision :: M_A                     ! mobility of phase-field
 double precision :: M_B

end module variables

subroutine input_data
 use variables
 implicit none

 integer i 
 real(8)::data(25)
 open(1,file="parameters.txt")
 read(1,'(23x,e10.3)') (data(i),i=1,23)
 close(1)

 step_end =    int(data(1))
 step_out =    int(data(2))
 grid_x   =    int(data(3))
 grid_y   =    int(data(4))
 dt       =        data(5)
 dx       =        data(6)
 latent_heat_A   = data(7)
 latent_heat_B   = data(8)
 melting_temp_A  = data(9)
 melting_temp_B  = data(10)
 int_energy_A    = data(11)
 int_energy_B    = data(12)
 kinetic_coeff_A = data(13)
 kinetic_coeff_B = data(14)
 diff_coeff_liquid = data(15)
 diff_coeff_solid  = data(16)
 vmol            = data(17)
 aniso_mode      = data(18)
 aniso_stre      = data(19) 
 magnitude_noise = data(20)
 theta_0         = data(21)
 initial_concentration = data(22)
 initial_temperature   = data(23)

 dy = dx
 theta_0  = theta_0*pi/180.0
 thick_diff_interface  = 6.*dx ! thickness of diffuse interface [m]
 bbb  = 2.*log((1.+(1.-2.*0.1))/(1.-(1.-2.*0.1)))/2.
 tinb = int_energy_A/int_energy_B*melting_temp_B/melting_temp_A*thick_diff_interface
 grad_energ_coeff = sqrt(3.*thick_diff_interface*int_energy_A/(bbb*melting_temp_A)) ! gradient energy coefficient
 W_A = 6.*int_energy_A*bbb/(thick_diff_interface*melting_temp_A)                    ! double-well potential height
 W_B = 6.*int_energy_B*bbb/(tinb*melting_temp_B)
 M_A = bbb*melting_temp_A*melting_temp_A*kinetic_coeff_A/(3.*thick_diff_interface*latent_heat_A) ! mobility of phase-field
end subroutine input_data

!================================
! Main program
!================================
program pf
 use variables
 implicit none

 integer :: nstep,iout   ! countors of simulation step and output step
 double precision, allocatable:: pp(:,:), ppp(:,:) ! phase-field variable at time t and t+dt
 double precision, allocatable:: cc(:,:), ccc(:,:) ! concentration variable at time t and t+dt

 call input_data
 allocate (pp(grid_x,grid_y))
 allocate (ppp(grid_x,grid_y))
 allocate (cc(grid_x,grid_y))
 allocate (ccc(grid_x,grid_y))

 call initial_data(pp,cc,ppp,ccc) ! set initial condition

 iout  = 0
 call output(iout,cc,pp)  ! output the initial condition

  do nstep = 1, step_end  ! start iteration
   call solve(nstep,pp,ppp,cc,ccc) ! solve time evolution equations
	 if(mod(nstep,step_out)==0) then
      iout=iout+1
      call output(iout,cc,pp)
 	    write(*,*) 'output No. = ', iout
   endif
  enddo

 deallocate (pp)
 deallocate (ppp)
 deallocate (cc)
 deallocate (ccc)

 end

!===================================
! initial setting of field variables
!===================================
subroutine initial_data(pp,cc,ppp,ccc)
 use variables
 implicit none

 double precision, dimension(grid_x,grid_y) :: pp, ppp
 double precision, dimension(grid_x,grid_y) :: cc, ccc
 integer :: i,j
 double precision :: x,y,r,radius_seed

 radius_seed = 3.*dx   ! radius of initial nuclei
 pp = 0.               ! set phase-field variable as 0
 cc = initial_concentration  ! set initial concentration of element B

	  do i=1,grid_x
          do j=1,grid_y
              x = dx*real(i-grid_x/2)  ! an initial nuclei is put at the center of region
              y = dy*real(j-grid_y/2)
              r = sqrt(x*x + y*y)
              pp(i,j) = 0.5*(1.-tanh(sqrt(2.*W_A)/(2.*grad_energ_coeff)*(r-radius_seed))) ! see Eq.(2.33)
              if(pp(i,j).le.1.e-05) pp(i,j)=0.
          enddo
	   enddo

  ppp = pp
  ccc = cc

 return
end subroutine initial_data

!===================================
! solve time evolution equations
!===================================
subroutine solve(nstep,pp,ppp,cc,ccc)
 use variables
 implicit none

 double precision, dimension(grid_x,grid_y) :: pp, ppp
 double precision, dimension(grid_x,grid_y) :: cc, ccc
 integer :: i,ip,im,j,jp,jm,k,iran,nstep
 double precision, dimension(4) :: dphidx,dphidy,dcdx,dcdy
 double precision, dimension(4) :: grad_energ_coeff_at,dgrad_energ_coeff_at
 double precision, dimension(4) :: dc,dp
 double precision :: dpt,dcc,phi,con
 double precision :: dpdx,dpdy,ddpp,theta,tcos,tsin
 double precision :: eepy1,eepy2,eepx1,eepx2,eppx1,eppx2,eppy1,eppy2
 double precision :: function_p,doublwell,H_A,H_B,dca1,dca2,dcb1,dcb2
 double precision :: ag,uran

 iran = nstep  ! set a seed for generating uniform random numbers

 do j=1,grid_y
     do i=1,grid_x

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      if(ip==grid_x+1) ip = grid_x  ! zero neumann boundary condition
      if(im==0       ) im = 1
      if(jp==grid_y+1) jp = grid_y
      if(jm==0       ) jm = 1

     !----- solve Allen-Cahn equation for phase-field variable -----
      phi = pp(i,j)
 	    con = cc(i,j)
      call calc_properties(phi,function_p,doublwell,H_A,H_B)

      dphidx(1)=(pp(ip ,j  )-pp(i  ,j  ))/dx  ! spacial gradient of phase-field variable at (i+1/2,j)
      dphidx(2)=(pp(i  ,j  )-pp(im ,j  ))/dx  !                                             (i-1/2,j)
      dphidx(3)=(pp(ip ,j  )+pp(ip ,jp )-pp(im ,j  )-pp(im ,jp ))/(4.*dx)
      dphidx(4)=(pp(ip ,j  )+pp(ip ,jm )-pp(im ,j  )-pp(im ,jm ))/(4.*dx)
      dphidy(1)=(pp(i  ,jp )+pp(ip ,jp )-pp(i  ,jm )-pp(ip ,jm ))/(4.*dy)
      dphidy(2)=(pp(i  ,jp )+pp(im ,jp )-pp(i  ,jm )-pp(im ,jm ))/(4.*dy)
      dphidy(3)=(pp(i  ,jp )-pp(i  ,j  ))/dy
      dphidy(4)=(pp(i  ,j  )-pp(i  ,jm ))/dy

      ! calculate normal direction of the interface and anisotropic interfacial energy
      do k=1,4
         dpdx = dphidx(k)  ! see Eq.(3.1)
         dpdy = dphidy(k)
	       ddpp = sqrt(dpdx*dpdx+dpdy*dpdy)
         if(dpdx.eq.0..and.dpdy.eq.0.) then
           theta = 0.
         else
           tcos=-dpdx/ddpp
           tsin=-dpdy/ddpp
           if(tsin.ge.0.) theta=      acos(tcos)     ! angle of the interface normal (see Eq.(3.2))
	         if(tsin.le.0.) theta=2.*pi-acos(tcos)
	         if(dpdx.eq.0..and.(-dpdy).gt.0.) theta=pi/2.
	         if(dpdx.eq.0..and.(-dpdy).lt.0.) theta=3.*pi/2.
	         if(dpdy.eq.0..and.(-dpdx).gt.0.) theta=0.
	         if(dpdy.eq.0..and.(-dpdx).lt.0.) theta=pi
         endif
          grad_energ_coeff_at(k) = grad_energ_coeff*(1.+aniso_stre*cos(aniso_mode*(theta-theta_0)))         ! gradient energy coefficnet (see Eq.(3.20))
         dgrad_energ_coeff_at(k) = -aniso_mode*grad_energ_coeff*aniso_stre*sin(aniso_mode*(theta-theta_0))  ! derivative of grad_energ_coeff w.r.t. phase-field variable
      enddo

      eepy1 = grad_energ_coeff_at(1)*dgrad_energ_coeff_at(1)*dphidy(1)
      eepy2 = grad_energ_coeff_at(2)*dgrad_energ_coeff_at(2)*dphidy(2)
      eepx1 = grad_energ_coeff_at(3)*dgrad_energ_coeff_at(3)*dphidx(3)
      eepx2 = grad_energ_coeff_at(4)*dgrad_energ_coeff_at(4)*dphidx(4)
      eppx1 = grad_energ_coeff_at(1)*grad_energ_coeff_at(1)*dphidx(1)
      eppx2 = grad_energ_coeff_at(2)*grad_energ_coeff_at(2)*dphidx(2)
      eppy1 = grad_energ_coeff_at(3)*grad_energ_coeff_at(3)*dphidy(3)
      eppy2 = grad_energ_coeff_at(4)*grad_energ_coeff_at(4)*dphidy(4)

      call noise(ag,iran)  ! get noise to generate secondary arm of dendrite

      ! right hand side of Allen-Cahn equation (see Eq.(5.32))
      dpt = ((1.-con)*M_A+con*M_B)*(-(eepy1-eepy2)/dx+(eepx1-eepx2)/dy &
                                    +(eppx1-eppx2)/dx+(eppy1-eppy2)/dy &
            -(1.-16.*ag*doublwell)*((1.-con)*H_A+con*H_B))

      ! time integration (calculate the phase-field variable at time t+dt)
      ppp(i,j) = pp(i,j) + dpt*dt
      if(ppp(i,j).le.1.e-20) ppp(i,j)=0.

      ! --------- solve Cahn-Hilliard equation for concentration field ----------
      dcdx(1)=(cc(ip,j  )-cc(i ,j ))/dx     ! spacial gradient of concentration variable at (i+1/2,j)
      dcdx(2)=(cc(i  ,j )-cc(im,j ))/dx
      dcdx(3)=(cc(ip,j  )+cc(ip,jp)-cc(im,j  )-cc(im,jp))/(4.*dx)
      dcdx(4)=(cc(ip,j  )+cc(ip,jm)-cc(im,j  )-cc(im,jm))/(4.*dx)
      dcdy(1)=(cc(i  ,jp)+cc(ip,jp)-cc(i  ,jm)-cc(ip,jm))/(4.*dy)
      dcdy(2)=(cc(i  ,jp)+cc(im,jp)-cc(i  ,jm)-cc(im,jm))/(4.*dy)
      dcdy(3)=(cc(i  ,jp)-cc(i  ,j ))/dy
      dcdy(4)=(cc(i  ,j )-cc(i  ,jm))/dy

      phi = (pp(ip,j  )+pp(i,j))/2.
      con = (cc(ip,j  )+cc(i,j))/2.
      call calc_properties(phi,function_p,doublwell,H_A,H_B)
      dc(1) = diff_coeff_liquid+function_p*(diff_coeff_solid-diff_coeff_liquid)
      dp(1) = dc(1)*con*(1.-con)*vmol/Rgas*(H_A-H_B)    ! First term of Eq.(5.36) at (i+1/2,j)

      phi=(pp(im,j  )+pp(i,j))/2.
      con=(cc(im,j  )+cc(i,j))/2.
      call calc_properties(phi,function_p,doublwell,H_A,H_B)
      dc(2) = diff_coeff_liquid+function_p*(diff_coeff_solid-diff_coeff_liquid)
      dp(2) = dc(2)*con*(1.-con)*vmol/Rgas*(H_A-H_B)    ! First term of Eq.(5.36) at (i-1/2,j)

      phi=(pp(i  ,jp)+pp(i,j))/2.
      con=(cc(i  ,jp)+cc(i,j))/2.
      call calc_properties(phi,function_p,doublwell,H_A,H_B)
      dc(3)=diff_coeff_liquid+function_p*(diff_coeff_solid-diff_coeff_liquid)
      dp(3) =dc(3)*con*(1.-con)*vmol/Rgas*(H_A-H_B)     ! First term of Eq.(5.36) at (i,j+1/2)

      phi=(pp(i  ,jm)+pp(i,j))/2.
      con=(cc(i  ,jm)+cc(i,j))/2.
      call calc_properties(phi,function_p,doublwell,H_A,H_B)
      dc(4)=diff_coeff_liquid+function_p*(diff_coeff_solid-diff_coeff_liquid)
      dp(4) =dc(4)*con*(1.-con)*vmol/Rgas*(H_A-H_B)     ! First term of Eq.(5.36) at (i,j-1/2)

      dca1=dp(1)*dphidx(1)-dc(1)*dcdx(1)    ! Second term of Eq. (5.36) at (i+1/2,j)
      dca2=dp(2)*dphidx(2)-dc(2)*dcdx(2)
      dcb1=dp(3)*dphidy(3)-dc(3)*dcdy(3)
      dcb2=dp(4)*dphidy(4)-dc(4)*dcdy(4)

      ! right hand side of Cahn-Hilliard equation (see. Eq.(5.36))
      dcc=-((dca1-dca2)/dx+(dcb1-dcb2)/dy)

      ! time integration (calculate concentration variable variable at time t+dt)
      ccc(i,j) = cc(i,j) + dcc*dt

   enddo
  enddo

  ! update of order parameters(phase-field and concentration variables)
  pp = ppp
  cc = ccc

 return
end

!===================================
! calculate thermal properties
!===================================
subroutine calc_properties(phi,function_p,doublwell,H_A,H_B)
  use variables
  implicit none

  double precision:: phi
  double precision:: function_p,doublwell,dpdphi,dqdphi_A,dqdphi_B,H_A,H_B

  function_p = phi*phi*phi*(10.-15.*phi+6.*phi*phi)  ! Energy density function (see. Eq.(5.23))
  doublwell = phi*phi*(1.-phi)*(1.-phi)              ! Double well function (see Eq.(5.24))
  dpdphi    = 30.*doublwell                          ! derivative of energy density function w.r.p. phase-field variable
  dqdphi_A  = W_A*2.*phi*(1.-phi)*(1.-2.*phi)        ! second term of right hand side of Eq.(5.30)
  dqdphi_B  = W_B*2.*phi*(1.-phi)*(1.-2.*phi)
  H_A = dqdphi_A+dpdphi*latent_heat_A*(initial_temperature-melting_temp_A)/(initial_temperature*melting_temp_A)  ! see Eq.(5.30)
  H_B = dqdphi_B+dpdphi*latent_heat_B*(initial_temperature-melting_temp_B)/(initial_temperature*melting_temp_B)

  return
end

!===================================
! generation of noise
!===================================
subroutine noise(ag,iran)
  use variables
  implicit none

  integer :: iran,l1,l2,l3
  real(8) :: t3
  double precision :: uran,ag

    iran=iran+1
    l1=843314861
    l2=453816693
    l3=2**30
    t3=2.0**31
    iran=l1*iran+l2
    if(iran.lt.0) iran=(iran+l3)+l3
    uran=real(iran)/t3
    uran=2.*(uran-0.5)
    ag=magnitude_noise*uran

  return
end

!===================================
! file output of results
! VTK file can be visualized by ParaView software.
! ParaView can be downloaded at https://www.paraview.org/
!===================================
subroutine output(iout,cc,pp)
 use variables
 implicit none

 integer :: iout
 double precision, dimension(grid_x,grid_y) :: pp
 double precision, dimension(grid_x,grid_y) :: cc

 character*30::filename
 integer :: i,j

 write(filename,'(a,i4.4,a)') 'pf_result',iout,'.vtk'
 open(101,file=filename)
 write(101,'(a)') '# vtk DataFile Version 3.0'
 write(101,'(a)') 'output.vtk'
 write(101,'(a)') 'ASCII'
 write(101,'(a)') 'DATASET STRUCTURED_POINTS'
 write(101,'(a,3i5)') 'DIMENSIONS',grid_x,grid_y,1
 write(101,'(a,3f4.1)')'ORIGIN' ,0.0,0.0,0.0
 write(101,'(a,3i2)')'ASPECT_RATIO',1,1,1
 write(101,'(a,1i11)')'POINT_DATA',grid_x*grid_y*1
 write(101,'(a)')'SCALARS concentration double'
 write(101,'(a)')'LOOKUP_TABLE default'
  do j=1,grid_y
     do i=1,grid_x
         write(101,*) cc(i,j)
     end do
  end do
 write(101,'(a)')'SCALARS phase_field double'
 write(101,'(a)')'LOOKUP_TABLE default'
 do j=1,grid_y
     do i=1,grid_x
		write(101,*) pp(i,j)
	 end do
  end do
  close(101)

 return
end
