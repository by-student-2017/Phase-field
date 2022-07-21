!----------------------------------------------------------------------
!
!     Phase-field calculation program for binary alloy
!     For free dendrite analysis
!
!     00/11/14 MACHIKO ODE
!     This is a shareware program
!
!----------------------------------------------------------------------
      
      module common_v
      
      integer::l=0,i=1,j=1,Nx,Ny,Nxc,Nyc
      integer::triangle=15,modsave=2000,lsave=0
      real(8)::M,ep,ep2,W
      real(8)::dx,dy,dt,dt1,Ds,Dl,Dl2,Ds2
      real(8)::Cle,Cse,Tmelt,T,Vm,co
      real(8)::v,yk,sigma,me,ke,beta,comnoise
      
      real(8),allocatable::phi(:,:,:), com(:,:,:)
      real(8),allocatable::fcl(:,:),fcc(:,:),Cl(:,:),Cs(:,:)
      
      real(8),parameter::R=8.314   ! R = Gas constant
      
      end module common_v
      
!======================================================================
!     main program
!======================================================================
      
      program main
      
      use common_v
      implicit none
      
      real(8)::pd,phix,phiy,phixx,phiyy,phixy
      real(8)::e1,e2,e3,e4,e5,th,eta,deta
      real(8)::p,hp,pg,gd,gg,fp,dc,c,dphi
      real(8)::fccw,fcce,fccn,fccs,fccl
      real(8)::xj1,xj2,xj3,xj4
      real(8)::d1,d2,d3,d4,d5
      real(8)::a,aa,bb,cc
      
!----------------------------------------------------------------------
!     read calculation condition
!----------------------------------------------------------------------
      
      call cal_cond
      
!----------------------------------------------------------------------
!     set the initial phase-condition
!----------------------------------------------------------------------
      
      call init_cond
      
      allocate( Cs(0:Nx+1,0:Ny+1) ); allocate( Cl(0:Nx+1,0:Ny+1) )
      allocate( fcl(0:Nx+1,0:Ny+1) ); allocate( fcc(0:Nx+1,0:Ny+1) )
      
!----------------------------------------------------------------------
!     set the phase-field parameters (& timestep)
!----------------------------------------------------------------------
      ! ep = epsilon = Interface energy anisotropy, lambda = 3.0*dx
      ep=sqrt(18.0/2.2*dx*sigma); ep2=ep*ep
      W=2.2*sigma/dx   ! W = Coefficient of excess free energy at the interface
      call mobility
      a=Cle/Cse*(1.0-Cse)/(1.0-Cle)   ! a = Isochemical potential conditions
      
      write(6,*) "set all the calculation conditions"
      write(6,*) "now calculating ..."
      
!======================================================================
!     calculate the governing equations
!======================================================================
      
  500 l=l+1   ! number of cycle
      
!----------------------------------------------------------------------
!     calculate Cs & Cl
!----------------------------------------------------------------------
      
      do i=0, Nxc+1
        do j=0, Nyc+1
          
          p=phi(0,i,j)
          hp=p**3*( 6.0*p*p - 15.0*p + 10.0 )   ! h(phi) = hp
          
          if( phi(0,i,j) .lt. 0.001) then
            Cl(i,j)=com(0,i,j)
            Cs(i,j)=Cl(i,j)/(a+(1.0-a)*Cl(i,j))
            
          else if( phi(0,i,j) .gt. 0.999) then
            Cs(i,j)=com(0,i,j)
            Cl(i,j)=a*Cs(i,j)/(1.0+(a-1.0)*Cs(i,j))
            
          else
            aa=hp*(1.0-a)
            bb=hp + (1.0-hp)*a + com(0,i,j)*(1.0-a)
            cc=com(0,i,j)
            
            Cs(i,j)=( bb - sqrt(bb*bb-4.0*aa*cc) )/(2.0*aa)
            Cl(i,j)=( cc - hp*Cs(i,j) )/(1.0-hp)
          end if
          ! Vm = Molar volume, T = temperature
          fcl(i,j)=(R*T/Vm)*log( Cl(i,j)/(1.0 - Cl(i,j)) )
          fccl=(R*T/Vm)*(1/( Cl(i,j)*(1.0 - Cl(i,j)) ))        ! eq.8 for i=L
          fccs=(R*T/Vm)*(1/( Cs(i,j)*(1.0 - Cs(i,j)) ))        ! eq.8 for i=S
          fcc(i,j)=(fccl*fccs)/( (1.0 - hp)*fccs + hp*fccl )   ! eq.8
          
        end do
      end do
      
!----------------------------------------------------------------------
!     governing equations
!----------------------------------------------------------------------
      
      do i=1, Nxc
        do j=1, Nyc
          
          p=phi(0,i,j); c=com(0,i,j)
          
!-----------------------------------------------------------------------
!     time saving
!-----------------------------------------------------------------------
          
          pd=( phi(0,i+1,j) + phi(0,i-1,j) + phi(0,i,j+1) + phi(0,i,j-1) )/4.0
          
          if( pd .le. 1.0e-5 ) then
            
            dphi=0.0
            dc=Dl*(( com(0,i,j+1) + com(0,i,j-1) - 2.0*c )/(dy*dy) &
     &            +( com(0,i+1,j) + com(0,i-1,j) - 2.0*c )/(dx*dx))
            ! Dl = diffusion coefficient of Liquid
          else if( pd .ge. (1.0-1.0e-5) ) then
            
            dphi=0.0
            dc=Ds*(( com(0,i,j+1) + com(0,i,j-1) - 2.0*c )/(dy*dy) &
     &            +( com(0,i+1,j) + com(0,i-1,j) - 2.0*c )/(dx*dx))
            ! Ds = diffusion coefficient of Solid
!-----------------------------------------------------------------------
!     non-time saving, eq.5 and eq.6
!-----------------------------------------------------------------------
          
          else
            
            pg=30*p*p*(1.0 - p)*(1.0 - p)
            gd=2.0*p*(1.0 - p)*(1.0 - 2.0*p)
            
            gg=pg*log( (1.0 - Cse)/(1.0 - Cle)*(1.0 - Cl(i,j))/(1.0 - Cs(i,j)) )
            
            fp=(R*T/Vm)*gg-W*gd   ! d(eq.1)/d(phi)
            
            phix=( phi(0,i-1,j) - phi(0,i+1,j) )/(2.0*dx)
            phiy=( phi(0,i,j-1) - phi(0,i,j+1) )/(2.0*dy)
            
            phixx=( phi(0,i-1,j) + phi(0,i+1,j) - 2.0*p )/(dx*dx)
            phiyy=( phi(0,i,j-1) + phi(0,i,j+1) - 2.0*p )/(dy*dy)
            phixy=( phi(0,i+1,j+1) + phi(0,i-1,j-1) - phi(0,i-1,j+1) - phi(0,i+1,j-1) )/(2.0*dx*2.0*dy)
            th=atan( phiy/(phix + 1.0e-20) )
            
            
            eta=1.0+v*cos(yk*th)   ! related with eq.4
            deta=v*yk*sin(yk*th)   ! -d(eta)/d(th)
            
            e1=ep2*eta*eta*( phixx + phiyy )
            e2=ep2*eta*(-deta)*( sin(2.0*th) * (phiyy - phixx) + 2.0*cos(2.0*th)*phixy )
            e3=0.5*ep2
            e4=deta*deta+eta*( -v*yk*yk*cos(yk*th) )
            e5=2.0*sin(2.0*th)*phixy - phixx - phiyy - cos(2.0*th)*(phiyy - phixx)
            
            ! M = mobility
            dphi=M*( e1 + e2 - e3*e4*e5 + fp )   ! eq.5
            
            d1=Ds; d2=Ds; d3=Ds; d4=Ds; d5=Ds
            
            if( p .le. 0.9 ) d1=Dl
            if( phi(0,i-1,j) .le. 0.9) d2=Dl; if( phi(0,i+1,j) .le. 0.9) d3=Dl
            if( phi(0,i,j+1) .le. 0.9) d4=Dl; if( phi(0,i,j-1) .le. 0.9) d5=Dl
            
            fccw=2.0*d1/fcc(i,j)*d2/fcc(i-1,j)/( d1/fcc(i,j) + d2/fcc(i-1,j) )
            fcce=2.0*d1/fcc(i,j)*d3/fcc(i+1,j)/( d1/fcc(i,j) + d3/fcc(i+1,j) )
            fccs=2.0*d1/fcc(i,j)*d4/fcc(i,j+1)/( d1/fcc(i,j) + d4/fcc(i,j+1) )
            fccn=2.0*d1/fcc(i,j)*d5/fcc(i,j-1)/( d1/fcc(i,j) + d5/fcc(i,j-1) )
            
            xj1=( -fcl(i,j) + fcl(i-1,j) )/dx*fccw
            xj2=( -fcl(i,j) + fcl(i+1,j) )/dx*fcce
            xj3=( -fcl(i,j) + fcl(i,j+1) )/dy*fccs
            xj4=( -fcl(i,j) + fcl(i,j-1) )/dy*fccn
            
            dc=( xj1 + xj2 )/dx + ( xj3 + xj4 )/dy   ! eq.6
            
          end if
          
          phi(1,i,j)=p+dphi*dt; com(1,i,j)=c+dc*dt
          
        end do
      end do
      
!-----------------------------------------------------------------------
!     end governing quation calculations
!-----------------------------------------------------------------------
!     boundary conditions
!-----------------------------------------------------------------------
      
      do i=0, Nx+1
        phi(1,i,0)=phi(1,i,1); phi(1,i,Ny+1)=phi(1,i,Ny)
        com(1,i,0)=com(1,i,1); com(1,i,Ny+1)=com(1,i,Ny)
      end do
      
      do j=0, Ny+1
        phi(1,0,j)=phi(1,1,j); phi(1,Nx+1,j)=phi(1,Nx,j)
        com(1,0,j)=com(1,1,j); com(1,Nx+1,j)=com(1,Nx,j)
      end do
      
!-----------------------------------------------------------------------
!     renewal of phase & concentration fields
!-----------------------------------------------------------------------
      
      do i=0, Nx+1
        do j=0, Ny+1
          phi(0,i,j)=phi(1,i,j); com(0,i,j)=com(1,i,j)
        end do
      end do
      
!-----------------------------------------------------------------------
!     noise
!-----------------------------------------------------------------------
      
      call noise
      
!-----------------------------------------------------------------------
!     area set for time saving
!-----------------------------------------------------------------------
      
      if( l .ge. 100 ) call areaset
      
!-----------------------------------------------------------------------
!     out put
!-----------------------------------------------------------------------
      
      if( mod(l,modsave) .eq. 0 ) call outsave 
      
!-----------------------------------------------------------------------
!     end condition
!-----------------------------------------------------------------------
      
      if( phi(0,1,Ny-10) .le. 0.5) goto 500
      write(6,*) "calculation has finished", l, "cycles"
      
      end
      
!-----------------------------------------------------------------------
!     subroutine
!========================================================================
!     read calculation condition
!-----------------------------------------------------------------------
      
      subroutine cal_cond
      
      use common_v
      implicit none
      real(8)::data(20)
      open(1,file="parameters.txt")
      read(1,'(11x,e9.3)') (data(i),i=1,19)
      close(1)
      !write(*,'(e9.3)') (data(i),i=1,19)   ! check input parameters
      modsave=int(data(1))  ! outout at every modsave
      triangle=int(data(2)) ! initial triangle
      ! e.g., Al-1.96 mol%Cu alloy
      Nx=int(data(3))  ! Number of mesh in x direction
      Ny=int(data(4))  ! Number of mesh in y direction
      dx=data(5)     ! mesh size, 1.0e-8 = 0.01 micro meter
      dy=data(6)     ! mesh size, 1.0e-8 = 0.01 micro meter
      Dl=data(7)     ! Di(i=L), diffusion coefficient of Liquid
      Ds=data(8)     ! Di(i=S), diffusion coefficient of Solid
      co=data(9)     ! initial solute content
      T =data(10)    ! initial temperature [K]
      Tmelt=data(11) ! melting point [K]
      Vm=data(12)    ! moler volume
      me=data(13)    ! liquidus slope
      ke=data(14)    ! partition coefficient
      beta=data(15)  ! kinetic coefficient
      v =data(16)    ! anisotropy, (1.0+v*cos(yk*theta))
      yk=data(17)    ! anisotropy, (1.0+v*cos(yk*theta))
      sigma=data(18) ! interface energy
      comnoise=data(19) ! noise
      
      Nyc=Ny; Nyc=Nx; Dl2=2.0*Dl/dx; Ds2=2.0*Ds/dx
      Cle=( Tmelt - T )/me; Cse=Cle*ke
      
      end
      
!-----------------------------------------------------------------------
!     init_cond     Default value setting
!-----------------------------------------------------------------------
      
      subroutine init_cond
      use common_v
      implicit none
      
      allocate( phi(0:1,0:Nx+1, 0:Ny+1) )
      allocate( com(0:1,0:Nx+1, 0:Ny+1) )
      
      do i=1, Nx
        do j=1, Ny
          
          phi(0,i,j)=0.0; phi(1,i,j)=0.0
          com(0,i,j)=co ; com(1,i,j)=co
          
        end do
      end do
      ! Set the initial solid phase of the triangle near the origin.
      do i=1, triangle+5
        do j=1, triangle+5
          
          if( j .lt. (-i+triangle) ) then
            phi(0,i,j)=1.0; phi(1,i,j)=1.0
            com(0,i,j)=Cse; com(1,i,j)=Cse
            
          end if
          
          if( j .eq. (-i+triangle) ) then
            phi(0,i,j)=0.5; phi(1,i,j)=0.5
          end if
          
        end do
      end do
      
      do i=0, Nx+1
        phi(0,i,0)=phi(0,i,1); phi(0,i,Ny+1)=phi(0,i,Ny)
        com(0,i,0)=com(0,i,1); com(0,i,Ny+1)=com(0,i,Ny)
      end do
      
      do j=0, Ny+1
        phi(0,0,j)=phi(0,1,j); phi(0,Nx+1,j)=phi(0,Nx,j)
        com(0,0,j)=com(0,1,j); com(0,Nx+1,j)=com(0,Nx,j)
      end do
      
      return
      end
      
!-----------------------------------------------------------------------
!     pf mobility
!-----------------------------------------------------------------------
      subroutine mobility
      use common_v
      implicit none
      
      real(8)::p1,p2,hp1,hp2,fun1,fun2,zeta,fccle,fccse,alpha
      integer::ip1
      
      ep2=ep*ep; zeta=0.0
      
      fccle=(R*T/Vm)*(1/( Cle*(1.0 - Cle) ))
      fccse=(R*T/Vm)*(1/( Cse*(1.0 - Cse) ))
      
      do ip1=1, 998 !do p1=0.001, 0.998, 0.001 !Deleted feature: Loop variable at (1) must be integer
        p1=float(ip1)/1000.0
        p2=p1+0.001
        hp1=p1**3*( 6.0*p1*p1 - 15.0*p1 + 10.0 )   ! h(phi)
        hp2=p2**3*( 6.0*p2*p2 - 15.0*p2 + 10.0 )   ! h(phi)
        fun1=hp1*(1.0 - hp1)/( (1.0 - hp1)*fccse + hp1*fccle )*(1/( p1*(1.0 - p1) ))  ! zeta function of eq.14
        fun2=hp2*(1.0 - hp2)/( (1.0 - hp2)*fccse + hp1*fccle )*(1/( p2*(1.0 - p2) ))  ! zeta function of eq.14
        zeta=zeta+(fun1 + fun2)/2.0*0.001   ! d(phi) = 0.001
        ! integral for zeta function of eq.14
      end do
      
      alpha=(R*T/Vm)*( (1.0 - ke)/me )*beta
      M=1.0/( (ep*ep/sigma)*(alpha + ( ep/(dl*sqrt(2.0*W)) )*fccse*fccle*((Cle - Cse)**2)*zeta) )
      ! eq.13 and eq.14
      
      dt=dx**2/(5.0*M*ep**2)
      dt1=dx**2/(5.0*dl)
      dt=dmin1(dt,dt1)
      
      return
      end
      
!-----------------------------------------------------------------------
!     comnoise, Setting noise fluctuation for density.
!-----------------------------------------------------------------------
      
      subroutine noise
      use common_v
      implicit none
      real(8)::countnoise, comtot, comdam, cnoise
      integer::lam=12869,c=6925,mmu=32768,x=19724
      
      countnoise=0.0; comtot=0.0
      
      do i=1, Nxc
        do j=1, Nyc
          
          if( phi(0,i,j) .gt. 0.01 .and. phi(0,i,j) .le. 0.5 ) then
            
            comdam=com(0,i,j)
            
            x=mod( (x*lam+c), mmu )                 ! random seed
            cnoise=( real(x)/mmu - 0.5 )*comnoise   ! set random value
            
            com(0,i,j)=comdam*( 1.0 + cnoise )
            comtot=comtot+comdam*cnoise
            countnoise=countnoise+1.0
            
          end if
          
        end do
      end do
      
      do i=1, Nxc
        do j=1, Nyc
          
          if( phi(0,i,j) .gt. 0.01 .and. phi(0,i,j) .le. 0.5 ) then
            
            com(0,i,j)=com(0,i,j)-(comtot/countnoise)   ! Guarantees solute preservation
          
          end if
          
        end do
      end do
      
      do i=0, Nx+1
        
        com(0,i,0)=com(0,i,1); com(0,i,Ny+1)=com(0,i,Ny)
        
      end do
      
      do j=0, Ny+1
        
        com(0,0,j)=com(0,1,j); com(0,Nx+1,j)=com(0,Nx,j)
        
      end do
      
      return
      end
      
!-----------------------------------------------------------------------
!     area set
!-----------------------------------------------------------------------
      
      subroutine areaset
      use common_v
      implicit none
      
      do i=1, Nx
        if( abs(com(0,i,1)/co - 1.0 ) .gt. 1.0e-5 ) Nxc=i
      end do
      
      do j=1, Ny
        if( abs(com(0,1,j)/co - 1.0 ) .gt. 1.0e-5 ) Nyc=j
      end do
      
      Nxc=Nxc+10; Nyc=Nyc+10
      if( Nxc .gt. Nx ) Nxc=Nx; if( Nyc .gt. Ny ) Nyc=Ny
      
      return
      
      end
      
!-----------------------------------------------------------------------
!     outsave
!-----------------------------------------------------------------------
      
      subroutine outsave
      use common_v
      implicit none
      
      character*3::out_num
      character*15::fpout
      integer::one,ten,hand
      
      lsave=lsave+1
      one=mod(lsave,10)
      ten=mod( int(real(lsave)/10.0),10 )
      hand=mod( int(real(lsave)/100.0),10 )
      
      one=48+one; ten=48+ten; hand=48+hand;
      out_num=char(hand)//char(ten)//char(one)
      
      fpout="output"//out_num//".dat"
      open(14, file=fpout, err=1000)
      write(6,*) fpout
      write(14,410) Nxc,Nyc,l
      
      do i=1, Nxc
        do j=1, Nyc
          write(14,300) i,j,phi(0,i,j),com(0,i,j)  ! phi = phase-field, com = concentration
        end do
        write(14,*) ""
      end do
      
      close(14)
      
      return
      
  410 format("# ",i4,i4,i10,/)
  300 format(i5,i5,e12.5,e12.5)
 1000 write(6,*) "Error in file open"
      
      end