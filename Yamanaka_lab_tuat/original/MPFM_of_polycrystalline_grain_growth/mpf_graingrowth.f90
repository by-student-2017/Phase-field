!
! Purpose   :: Multi-phase-field model for ideal grain growth
! Reference :: I. Steinbach et al., Physica D, 134 (1999), pp. 385-393.
!              * Please refere the following textbook for equations numbered in this source code
!              T. Takaki and A. Yamanaka, Phase-field method, Yokendo-Ltd., (2012), Section 9 (in Japanese)
! Programmer:: Akinori Yamanaka (Tokyo Univ. Agri. Tech.)
! Date      :: Jan., 2019
!
! << Please note that programmer does not take any responsibility or liability for any damage or loss caused by using this program >>
!

!================================
! Main program
!================================
program mpf
implicit none

!================================
! Setting of model parameters
!================================
integer, parameter:: num_nucleus= 30       !! number of initial nucleus
integer, parameter:: num_grains = 1+num_nucleus  !! total number of grains
integer, parameter:: grid_x = 256          !! number of computational grids along x axis
integer, parameter:: grid_y = 256
integer, parameter:: step_end = 20000      !! total number of time steps
integer, parameter:: step_out = 2000       !! interval of output
double precision, parameter:: dx = 1.0e-06 !! spacing of finite difference grid along x axis [m]
double precision, parameter:: dy = dx
double precision, parameter:: dt = 1.0e-5                   !! time increment for a time step [s]
double precision, parameter:: r0 = 3.*dx                    !! radius of the initial nucleus
double precision, parameter:: thick_diff_interface = 5.*dx  !! thickness of diffuse interface
double precision, parameter:: gb_energy   = 0.5             !! grain boundary energy [J/m2]
double precision, parameter:: gb_mobility = 1.0e-8          !! grain boundary mobilitiy
double precision, parameter:: pi = 3.141592654
double precision, parameter:: grad_energ_coeff = 2./pi*sqrt(2.*thick_diff_interface*gb_energy)  !! gradient energy coefficient (see Eq.(9.35))
double precision, parameter:: height_W = 4.*gb_energy/thick_diff_interface                      !! height of double obstacle potential (see Eq.(9.36))
double precision, parameter:: pf_mobility = gb_mobility*pi*pi/(8.*thick_diff_interface)         !! phase-field mobility (see Eq.(9.37))
double precision, parameter:: eee = 5.0e+05                 !! magnitude of driving force [Pa]
character(len=12) filename

!====================================
! Declaration of arrays and variables
!====================================
integer, dimension(num_grains,grid_x,grid_y):: grain_id             !! index of phase-field variable
integer, dimension(grid_x,grid_y):: no_of_grain                     !! grain number at a local coordinate
double precision, dimension(num_grains,grid_x,grid_y):: p,pp        !! phase-field variables at time t and t+dt
double precision, dimension(num_grains,num_grains):: aij,wij,mij    !! phase-field parameters
double precision, dimension(num_grains):: nuclei_x,nuclei_y         !! x and y position of the initial grains
double precision, dimension(grid_x,grid_y):: e                      !! driving force to make the initial structure

integer:: i,j,k,l,lp,lm,m,mp,mm,n,n1,n2,n3,nf,iout,step
double precision:: phi,r,a,p_tmp,dpi,ppp,randx,randy

!==================================================================
! initial setting of phase field parameters (see Eqs.(9.42)-(9.45))
!==================================================================
wij = height_W
aij = grad_energ_coeff
mij = pf_mobility
e   = 0.0

do i=1,num_grains
  do j=1,num_grains
     if(i.eq.j) then
       wij(i,j)=0.
       aij(i,j)=0.
       mij(i,j)=0.
         e(i,j)=0.
     endif
     if(i.eq.num_grains.or.j.eq.num_grains) e(i,j) = eee
     if(i.gt.j) e(i,j) = -e(i,j)
  enddo
enddo

!===================================
! initial setting of field variables
!===================================
p  = 0.
pp = 0.

p(num_grains,:,:) = 1.          !! phase-field of mother phase
grain_id(1,:,:)   = num_grains  !! grain number of mother phase
no_of_grain(:,:)  = 1

do i=1,num_nucleus              !! position of initial nucleus are determined by random numbers
  call random_number(randx)
  call random_number(randy)
  nuclei_x(i)=int(randx*grid_x)
  nuclei_y(i)=int(randy*grid_y)
enddo

do i=1,num_nucleus
 do m=1,grid_y
  do l=1,grid_x
      r = sqrt((real(l-nuclei_x(i))*dx)**2+(real(m-nuclei_y(i))*dy)**2)
      a = sqrt(2.*height_W)/grad_energ_coeff*(r-r0)
      p_tmp = 0.5*(1.-sin(a))  !! see Eq.(9.41)
      if(a.ge. pi/2.) p_tmp=0.
      if(a.le.-pi/2.) p_tmp=1.
      if(p_tmp.gt.0.) then
        nf = no_of_grain(l,m)+1
        no_of_grain(l,m) = nf
        grain_id(nf,l,m) = i
        p(i,l,m) = p_tmp
        p(num_grains,l,m) = p(num_grains,l,m)-p_tmp
      endif
  enddo
 enddo
enddo

iout = 0
call output(iout,grid_x,grid_y,num_grains,p,no_of_grain,grain_id)  ! output of initial condition

! start iteration
step=0
do step = 1,step_end

! calculate number of grains at all grids
do m=1,grid_y
	do l=1,grid_x

 	 lp=l+1
	 lm=l-1
	 mp=m+1
	 mm=m-1
	 if(l.eq.grid_x) lp = 1       ! periodic boundary condition
   if(l.eq.1     ) lm = grid_x
	 if(m.eq.grid_y) mp = 1
	 if(m.eq.1     ) mm = grid_y

	 n=0
   do i=1,num_grains
     if(p(i,l  ,m ).gt.0..or.  &
       (p(i,l  ,m ).eq.0..and. &
        p(i,lp ,m ).gt.0..or.  &
        p(i,lm, m ).gt.0..or.  &
        p(i,l  ,mp).gt.0..or.  &
        p(i,l  ,mm).gt.0.))then
        n=n+1
        grain_id(n,l,m)=i
      endif
   enddo
   no_of_grain(l,m)=n

  enddo
enddo

! solve Allen-Cahn equation
do m=1,grid_y
 do l=1,grid_x

	lp=l+1
	lm=l-1
	mp=m+1
	mm=m-1
	if(l.eq.grid_x) lp = 1       ! periodic boundary condition
  if(l.eq.1     ) lm = grid_x
	if(m.eq.grid_y) mp = 1
	if(m.eq.1     ) mm = grid_y

  do n1 = 1,no_of_grain(l,m)   ! loop for index i
     i = grain_id(n1,l,m)
     dpi = 0.
    do n2 = 1,no_of_grain(l,m)   ! loop for index j
       j = grain_id(n2,l,m)
       ppp = 0.

       do n3 = 1,no_of_grain(l,m)   ! loop for index k
          k = grain_id(n3,l,m)
		  ! calculation of gradient term
		  ppp = ppp + (aij(i,k)*aij(i,k)-aij(j,k)*aij(j,k))/2.      &
                     *((p(k,lp,m )-2.*p(k,l,m)+p(k,lm,m ))/(dx*dx)  &
                      +(p(k,l ,mp)-2.*p(k,l,m)+p(k,l ,mm))/(dy*dy)) &
                     + (wij(i,k)-wij(j,k))*p(k,l,m)
       enddo

      ! right hand side of Allen-Cahn equation (see Eq.(9.29))
 	   dpi = dpi - 2.*mij(i,j)/real(no_of_grain(l,m))*(ppp - 8./pi*sqrt(p(i,l,m)*p(j,l,m))*e(i,j))
   enddo
	   ! time integration (calculate the phase-field variable at time t+dt)
     pp(i,l,m) = p(i,l,m)+dpi*dt
  enddo

 enddo
enddo

! adjust phase-field variables within 0~1
do m=1,grid_y
  do l=1,grid_x
    phi=0.
    do i=1,num_grains
      if(pp(i,l,m).le.0.) pp(i,l,m)=0.
      if(pp(i,l,m).ge.1.) pp(i,l,m)=1.
      phi=phi+pp(i,l,m)
    enddo
		a=1.
    if(phi.ne.1.) a=1./phi
    do i=1,num_grains
      p(i,l,m) =a*pp(i,l,m)
      pp(i,l,m)=0.
    enddo
  enddo
enddo

! output of result
if(mod(step,step_out).eq.0)then
	iout = iout+1
  write(*,*) 'output = ',iout
  call output(iout,grid_x,grid_y,num_grains,p,no_of_grain,grain_id)
endif

enddo ! for main iteration loop

end program mpf

!===================================
! file output of results
! VTK file can be visualized by ParaView software.
! ParaView can be downloaded at https://www.paraview.org/
!===================================
subroutine output(iout,grid_x,grid_y,num_grains,p,no_of_grain,grain_id)
implicit none

character(len=20) filename
integer, dimension(num_grains,grid_x,grid_y):: grain_id
integer, dimension(grid_x,grid_y):: no_of_grain
double precision, dimension(num_grains,grid_x,grid_y):: p

integer:: iout,grid_x,grid_y,num_grains,l,m,i,j,ori
double precision:: phi,gb

write(filename,'(a,i3.3,a)') 'mpf_result',iout,'.vtk'
open(1,file=filename)
write(1,'(a)') '# vtk DataFile Version 3.0'
write(1,'(a)') 'output.vtk'
write(1,'(a)') 'ASCII'
write(1,'(a)') 'DATASET STRUCTURED_POINTS'
write(1,'(a,3i5)') 'DIMENSIONS', grid_x,grid_y,1
write(1,'(a,3f4.1)') 'ORIGIN ', 0.0, 0.0, 0.0
write(1,'(a,3i2)') 'ASPECT_RATIO', 1, 1, 1
write(1,'(a,1i11)') 'POINT_DATA', grid_x*grid_y*1
write(1,'(a)') 'SCALARS grain_number float'
write(1,'(a)') 'LOOKUP_TABLE default'
do m=1,grid_y
 do l=1,grid_x
  phi = 0.
  do j=1,no_of_grain(l,m)
	  i = grain_id(j,l,m)
    if(p(i,l,m).gt.phi) then
      phi=p(i,l,m)
      ori=grain_id(j,l,m)
    endif
  enddo
  write(1,*) ori
 enddo
enddo
write(1,'(a)') 'SCALARS phase_field_1 float'
write(1,'(a)') 'LOOKUP_TABLE default'
do m=1,grid_y
do l=1,grid_x
	write(1,*) p(1,l,m)
enddo
enddo
write(1,'(a)') 'SCALARS grain_boundary float'
write(1,'(a)') 'LOOKUP_TABLE default'
do m=1,grid_y
do l=1,grid_x
 gb = 0.
 do j=1,num_grains
	gb = gb + p(j,l,m)*p(j,l,m)
 enddo
 write(1,*) gb
enddo
enddo

close(1)

return
end
