!
! Purpose   :: Phase-field model for martensitic transformation in 2D
! Reference :: Y. U. Wang and A. G. Khachaturyan, Acta Materialia, 45 (1999), p. 759.
!              * Please refere the following textbook for equations numbered in this source code
!              T. Takaki and A. Yamanaka, Phase-field method, Yokendo-Ltd., (2012), Sections 7 & 8 (in Japanese)
! Programmer:: Akinori Yamanaka (Tokyo Univ. Agri. Tech.)
! Date      :: Jan., 2019
!
! << Please note that programmer does not take any responsibility or liability for any damage or loss caused by using this program >>
!

!================================
! Main program
!================================
program martensite
 implicit none

!================================
! Setting of model parameters
!================================
      integer, parameter:: ig = 7
      integer, parameter:: grid = 128               !! number of computational grids (must be 2^ig)
      integer, parameter:: ndm = grid-1
      integer, parameter:: half_grid = grid/2
      integer, parameter:: step_end = 10000         !! number of time steps
      integer, parameter:: step_out = 500           !! interval of output
      double precision, parameter:: length = 256.e-09    !! length of computationa domain [m]
      double precision, parameter:: dx = length/grid     !! spacing of computational grids
      double precision, parameter:: dt = 0.02            !! time increment [dimensionless time]
      double precision, parameter:: temp = 355.          !! temperature [K]
      double precision, parameter:: pi = 3.141592
      double precision, parameter:: rr = 8.3145          !! gas constant
      double precision, parameter:: mobility = 1.0       !! mobility of phase-field
      double precision, parameter:: vm0 = 7.0e-6         !! molar volume [m3/mol]
      double precision, parameter:: aa0 = (3.5e+08*(405.-temp)/405.)*vm0/rr/temp     !! magnitude of driving force [dimensionless]
      double precision, parameter:: aa = 0.15            !! constant A for chemical free energy
      double precision, parameter:: bb = 3.*aa+12.       !! B
      double precision, parameter:: cc = 2.*aa+12.       !! C
      double precision, parameter:: grad_energ_coeff=(1.0e-15)/rr/temp/dx/dx  !! gradient energy coefficient [dimensionless]
      double precision, parameter:: stress_unit=(1.0e+11)*vm0/rr/temp   !! stress is regularized by this number
      double precision, parameter:: c11=1.4*stress_unit                 !! elastic constant: C11=140GPa
      double precision, parameter:: c22=c11
      double precision, parameter:: c33=c11
      double precision, parameter:: c44=0.28*stress_unit                !! C44=(C11-C12)/2 because of isotropic material
      double precision, parameter:: c55=c44
      double precision, parameter:: c66=c44
      double precision, parameter:: c12=0.84*stress_unit                !! C12=84GPa
      double precision, parameter:: c21=c12
      double precision, parameter:: c13=c12
      double precision, parameter:: c31=c12
      double precision, parameter:: c23=c12
      double precision, parameter:: c32=c12
      double precision, parameter:: ram0=c12                       !! Lame's constant ¥lambda
      double precision, parameter:: mu0=c44                        !! Lame's constant ¥mu
      double precision, parameter:: nu0=ram0/2./(ram0+mu0)         !! Poisson's ratio
      double precision, parameter:: munu0=2.*mu0*(1.-nu0)
      double precision, parameter:: sig22_a=0.                     !! applied stress along y direction
!      double precision, parameter::  sig22_a=2.e+06*vm0/rr/temp
      double precision, parameter:: ep11_a=-c12/(c11-c12)/(c11+2.*c12)*sig22_a      !! applied strain due to the applied stress
      double precision, parameter:: ep22_a=(c11+c12)/(c11-c12)/(c11+2.*c12)*sig22_a

!==========================================================
! Declaration of arrays and variables for phase field model
!==========================================================
      double precision, dimension(0:grid,0:grid):: p1,p2     !! phase-field variables 1 and 2
      double precision:: chem_pot1,elast_pot1,grad_pot1  !! chemical, elastic, and gradient potential for p1
      double precision:: chem_pot2,elast_pot2,grad_pot2
      double precision:: dpdt1,dpdt2                     !! right hand side of Allen-Chan equation for p1 and p2
      double precision, dimension(0:grid,0:grid):: eigen11_r,eigen22_r    !! eigen strain in real space
      double precision, dimension(0:grid,0:grid):: eigen12_r,eigen21_r
      double precision, dimension(0:grid,0:grid):: eigen11_f_real,eigen11_f_imag !! real and imagenary part of eigen strain in fourier space
      double precision, dimension(0:grid,0:grid):: eigen22_f_real,eigen22_f_imag
      double precision, dimension(0:grid,0:grid):: eigen12_f_real,eigen12_f_imag
      double precision, dimension(0:grid,0:grid):: eigen21_f_real,eigen21_f_imag
      double precision, dimension(0:grid,0:grid):: inhomo_strain11,inhomo_strain22   !! inhomogeneous strains
      double precision, dimension(0:grid,0:grid):: inhomo_strain12,inhomo_strain21
      double precision, dimension(0:grid,0:grid):: stress11,stress22                 !! Cauchy stress
      double precision, dimension(0:grid,0:grid):: stress33,stress12
      double precision, dimension(0:grid,0:grid):: elastic_strain11,elastic_strain22 !! elastic strain
      double precision, dimension(0:grid,0:grid):: elastic_strain12,elastic_strain21
      double precision, dimension(4,4):: eigen0_1,eigen0_2                       !! eigen strain for each variants
      double precision:: homo_strain11,homo_strain22,homo_strain12,homo_strain21 !! homogeneous strain
      double precision:: el_strain11,el_strain22
      double precision:: eigen11_p1,eigen22_p1,eigen11_p2,eigen22_p2
      double precision:: sigma11_r,sigma22_r,sigma11_i,sigma22_i,sigma12_r,sigma12_i
      double precision:: sum11,sum22,sum12,sum21
      double precision:: s1,s2,rand

!==========================================================
! Declaration of arrays and variables for Fourier transform
!==========================================================
      double precision, dimension(0:grid,0:grid):: xi,xr       !! imagenary and real part of variable to be transformed
      double precision, dimension(0:grid,0:grid):: xif,xrf     !! imagenary and real part of variable to be inverse transformed
      double precision, dimension(0:grid):: s,c              !! sin and cos
      double precision:: kxx,kyy,kxy,nnn,k_mag             !! wave vectors in Fourier space
      integer, dimension(0:grid):: ik(0:grid)

      integer:: i,j,k,l,ii,jj,kk,iii,jjj,m,n,ip,im,jp,jm
      integer:: nstep,iout

! make sin and cos table for FFT
      call sincos(pi,ig,grid,half_grid,c,s,ik)  !! pi, grid, half_grid→c,s,ik

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
         p1(i,j)=rand()       ! rand() generates random numbers ranging from 0 to 1
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
       enddo
      enddo

      call fft(-1,ig,grid,half_grid,ndm,c,s,ik,xr,xi)   !! forward FFT

      do i=0,ndm
       do j=0,ndm
       eigen11_f_real(i,j)=xr(i,j)  ! real part of eigen strain in Fourier space
       eigen11_f_imag(i,j)=xi(i,j)
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
       enddo
      enddo

      call fft(-1,ig,grid,half_grid,ndm,c,s,ik,xr,xi)

      do i=0,ndm
       do j=0,ndm
       eigen22_f_real(i,j)=xr(i,j)
       eigen22_f_imag(i,j)=xi(i,j)
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
       enddo
      enddo

      call fft(-1,ig,grid,half_grid,ndm,c,s,ik,xr,xi)

      do i=0,ndm
       do j=0,ndm
       eigen12_f_real(i,j)=xr(i,j)
       eigen12_f_imag(i,j)=xi(i,j)
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
        nnn=dsqrt(k_mag)         !! magnutude of wave vector |k|
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
       enddo
      enddo

      call fft(1,ig,grid,half_grid,ndm,c,s,ik,xr,xi)  !! inverse FFT

      inhomo_strain11 = xr    !! inhomogeneous strain in real space

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
       enddo
      enddo

      call fft(1,ig,grid,half_grid,ndm,c,s,ik,xr,xi)

      inhomo_strain22 =xr

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
       enddo
      enddo

      call fft(1,ig,grid,half_grid,ndm,c,s,ik,xr,xi)

      inhomo_strain12 = xr
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
       if(i.eq.ndm) ip=0      !! periodic boundary condition is assumed
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
end program martensite


!===============
! subroutines
!===============
      subroutine sincos(pi,ig,grid,half_grid,c,s,ik)
      implicit none

      integer:: it,it1,it2,mc,mn,ik(0:grid),grid,half_grid,ig
      double precision:: q,pi
      double precision:: c(0:grid),s(0:grid)

      q=2.*pi/grid
      do it=0,half_grid-1
       c(it)=dcos(q*it)
       s(it)=dsin(q*it)
      enddo

      ik(0)=0
      mn=half_grid
      mc=1
      do it1=1,ig
       do it2=0,mc-1
       ik(it2+mc)=ik(it2)+mn
       enddo
      mn=mn/2
      mc=2*mc
      enddo

      return
      end subroutine sincos

!============================
! Fast Fourier transformation
!============================
      subroutine fft(qs,ig,grid,half_grid,ndm,c,s,ik,xr,xi)
      implicit none

      integer:: qs,ig,grid,half_grid,ndm
      integer:: i,ic,ir,j,ik(0:grid)
      double precision:: xrf(0:grid),xif(0:grid)
      double precision:: xr(0:grid,0:grid),xi(0:grid,0:grid)
      double precision:: c(0:grid),s(0:grid)

      do ir=0,ndm
       do ic=0,ndm
        xrf(ic)=xr(ir,ic)
        xif(ic)=xi(ir,ic)
       enddo
       call fft1(qs,grid,ig,half_grid,xrf,xif,c,s)
       do ic=0,ndm
        xr(ir,ic)=xrf(ik(ic))
        xi(ir,ic)=xif(ik(ic))
       enddo
      enddo

      do ic=0,ndm
       do ir=0,ndm
        xrf(ir)=xr(ir,ic)
        xif(ir)=xi(ir,ic)
       enddo
       call fft1(qs,grid,ig,half_grid,xrf,xif,c,s)
       do ir=0,ndm
        xr(ir,ic)=xrf(ik(ir))
        xi(ir,ic)=xif(ik(ir))
       enddo
      enddo

      if(qs.gt.0) goto 100  ! if forward FFT, return

      ! if inverse FFT, the resutls must be divided by total number of grids
      do i=0,ndm
       do j=0,ndm
       xr(i,j)=xr(i,j)/grid/grid
       xi(i,j)=xi(i,j)/grid/grid
       enddo
      enddo

  100 continue
      end subroutine fft

!=======================
! 1-dimensional FFT
!=======================
      subroutine fft1(qs,grid,ig,half_grid,xrf,xif,c,s)
      implicit none

      integer:: half_grid,qs,ig,grid
      integer:: ix,ka,kb,l2,lf,mf,n2,nf
      double precision:: tj,tr,c(0:grid),s(0:grid)
      double precision:: xrf(0:grid),xif(0:grid)

      l2=1
      do lf=1,ig
      n2=half_grid/l2
       do mf=1,l2
        do nf=0,n2-1
        ix=nf*l2
        ka=nf+2*n2*(mf-1)
        kb=ka+n2
        tr=xrf(ka)-xrf(kb)
        tj=xif(ka)-xif(kb)
        xrf(ka)=xrf(ka)+xrf(kb)
        xif(ka)=xif(ka)+xif(kb)
        xrf(kb)=tr*c(ix)+tj*qs*s(ix)
        xif(kb)=tj*c(ix)-tr*qs*s(ix)
        enddo
       enddo
      l2=l2*2
      enddo

      return
      end subroutine fft1

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
