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
      
      integer::l=0,i=1,j=1,m,n,mc,nc
      integer::triangle=15,modsave=2000,lsave=0
      real(8)::xm,ep,ep2,w
      real(8)::dx,dy,dt,dt1,ds,dl,dl2,ds2
      real(8)::cle,cse,tmpmelt,tmp,vm,co
      real(8)::v,yk,sigma,xme,ke,beta,comnoise
      
      real(8),allocatable::phi(:,:,:), com(:,:,:)
      real(8),allocatable::fcl(:,:),fcc(:,:),cl(:,:),cs(:,:)
      
      real(8),parameter::r=8.314
      
      end module common_v
      
!======================================================================
!     main program
!======================================================================
      
      program main
      
      use common_v
      implicit none
      
      real(8)::pd,phix,phiy,phixx,phiyy,phixy
      real(8)::e1,e2,e3,e4,e5,th,eta,deta
      real(8)::p,pp,pg,gd,gg,fp,dc,c,dphi
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
      
      allocate( cs(0:m+1,0:n+1) ); allocate( cl(0:m+1,0:n+1) )
      allocate( fcl(0:m+1,0:n+1) ); allocate( fcc(0:m+1,0:n+1) )
      
!----------------------------------------------------------------------
!     set the phase-field parameters (& timestep)
!----------------------------------------------------------------------
      
      ep=sqrt(18.0/2.2*dx*sigma); ep2=ep*ep
      w=2.2*sigma/dx
      call mobility
      a=cle/cse*(1.0-cse)/(1.0-cle)
      
      write(6,*) "set all the calculation conditions"
      write(6,*) "now calculating ..."
      
!======================================================================
!     calculate the governing equations
!======================================================================
      
  500 l=l+1
      
!----------------------------------------------------------------------
!     calculate cs & cl
!----------------------------------------------------------------------
      
      do i=0, mc+1
        do j=0, nc+1
          
          p=phi(0,i,j)
          pp=p**3*(10.0-15.0*p+6.0*p*p)
          
          if( phi(0,i,j) .lt. 0.001) then
            cl(i,j)=com(0,i,j)
            cs(i,j)=cl(i,j)/(a+(1.0-a)*cl(i,j))
            
          else if( phi(0,i,j) .gt. 0.999) then
            cs(i,j)=com(0,i,j)
            cl(i,j)=a*cs(i,j)/(1.0+(a-1.0)*cs(i,j))
            
          else
            aa=pp*(1.0-a)
            bb=pp + (1.0-pp)*a + com(0,i,j)*(1.0-a)
            cc=com(0,i,j)
            
            cs(i,j)=( bb - sqrt(bb*bb-4.0*aa*cc) )/(2.0*aa)
            cl(i,j)=( cc - pp*cs(i,j) )/(1.0-pp)
          end if
          
          fcl(i,j)=r*tmp/vm*log( cl(i,j)/(1.0-cl(i,j)) )
          fccl=r*tmp/vm/( cl(i,j)/(1.0-cl(i,j)) )
          fccs=r*tmp/vm/( cs(i,j)/(1.0-cs(i,j)) )
          fcc(i,j)=fccl*fccs/( (1.0-pp)*fccs + pp*fccl )
          
        end do
      end do
      
!----------------------------------------------------------------------
!     governing equations
!----------------------------------------------------------------------
      
      do i=1, mc
        do j=1, nc
          
          p=phi(0,i,j); !=com(0,i,j)
          
!-----------------------------------------------------------------------
!     time saving
!-----------------------------------------------------------------------
          
          pd=( phi(0,i+1,j) + phi(0,i-1,j) + phi(0,i,j+1) + phi(0,i,j-1) )/4.0
          
          if( pd .le. 1.0e-5 ) then
            
            dphi=0.0
            dc=dl*( (com(0,i,j+1) + com(0,i,j-1) - 2.0*c)/(dy*dy) &
     &             +(com(0,i+1,j) + com(0,i-1,j) - 2.0*c)/(dx*dx) )
            
          else if( pd .ge. (1.0-1.0e-5) ) then
            
            dphi=0.0
            dc=ds*( (com(0,i,j+1) + com(0,i,j-1) - 2.0*c)/(dy*dy) &
     &             +(com(0,i+1,j) + com(0,i-1,j) - 2.0*c)/(dx*dx) )
          
!-----------------------------------------------------------------------
!     non-time saving
!-----------------------------------------------------------------------
          
          else
            
            pg=30*p*p*(1.0-p)*(1.0-p)
            gd=2.0*p*(1.0-p)*(1.0-2.0*p)
            
            gg=pg*log( (1.0-cse)/(1.0-cle)*(1.0-cl(i,j))/(1.0-cs(i,j)) )
            
            fp=r*tmp/vm*gg-w*gd
            
            phix=( phi(0,i-1,j) - phi(0,i+1,j) )/(2.0*dx)
            phiy=( phi(0,i,j-1) - phi(0,i,j+1) )/(2.0*dy)
            
            phixx=( phi(0,i-1,j) + phi(0,i+1,j) - 2.0*p )/(dx*dx)
            phiyy=( phi(0,i,j-1) + phi(0,i,j+1) - 2.0*p )/(dy*dy)
            phixy=( phi(0,i+1,j+1) + phi(0,i-1,j-1) - phi(0,i-1,j+1) - phi(0,i+1,j-1) )/(2.0*dx*2.0*dy)
            th=atan( phiy/(phix+1.0e-20) )
            
            
            eta=1.0+v*cos(yk*th)
            deta=v*yk*sin(yk*th)
            
            e1=ep2*eta*eta*( phixx + phiyy )
            e2=ep2*eta*(-deta)*( sin(2.0*th) * (phiyy-phixx) + 2.0*cos(2.0*th)* phixy )
            e3=0.5*ep2
            e4=deta*deta+eta*( -v*yk*yk*cos(yk*th) )
            e5=2.0*sin(2.0*th)*phixy - phixx - phiyy -cos(2.0*th)*( phiyy - phixx )
            
            
            dphi=xm*( e1 + e2 - e3*e4*e5 + fp )
            
            d1=ds; d2=ds; d3=ds; d4=ds; d5=ds
            
            if( p .le. 0.9 ) d1=dl
            if( phi(0,i-1,j) .le. 0.9) d2=dl; if( phi(0,i+1,j) .le. 0.9) d3=dl
            if( phi(0,i,j+1) .le. 0.9) d4=dl; if( phi(0,i,j-1) .le. 0.9) d5=dl
            
            fccw=2.0*d1/fcc(i,j)*d2/fcc(i-1,j)/( d1/fcc(i,j) + d2/fcc(i-1,j) )
            fcce=2.0*d1/fcc(i,j)*d3/fcc(i+1,j)/( d1/fcc(i,j) + d3/fcc(i+1,j) )
            fccs=2.0*d1/fcc(i,j)*d4/fcc(i,j+1)/( d1/fcc(i,j) + d4/fcc(i,j+1) )
            fccn=2.0*d1/fcc(i,j)*d5/fcc(i,j-1)/( d1/fcc(i,j) + d5/fcc(i,j-1) )
            
            xj1=( -fcl(i,j) + fcl(i-1,j) )/dx*fccw
            xj2=( -fcl(i,j) + fcl(i+1,j) )/dx*fcce
            xj3=( -fcl(i,j) + fcl(i,j+1) )/dy*fccs
            xj4=( -fcl(i,j) + fcl(i,j-1) )/dy*fccn
            
            dc=( xj1 + xj2 )/dx + ( xj3 + xj4 )/dy
            
          end if
          
          phi(1,i,j)=p+dphi*dt; com(1,i,j)=c+dc*dt
          
        end do
      end do
      
!-----------------------------------------------------------------------
!     end governing quation calculations
!-----------------------------------------------------------------------
!     boundary conditions
!-----------------------------------------------------------------------
      
      do i=0, m+1
        phi(1,i,0)=phi(1,i,1); phi(1,i,n+1)=phi(1,i,n)
        com(1,i,0)=com(1,i,1); com(1,i,n+1)=com(1,i,n)
      end do
      
      do j=0, n+1
        phi(1,0,j)=phi(1,1,j); phi(1,m+1,j)=phi(1,m,j)
        com(1,0,j)=com(1,1,j); com(1,m+1,j)=com(1,m,j)
      end do
      
!-----------------------------------------------------------------------
!     renewal of phase & concentration fields
!-----------------------------------------------------------------------
      
      do i=0, m+1
        do j=0, n+1
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
      
      if( phi(0,1,n-10) .le. 0.5) goto 500
      write(6,*) "calculation has finished"
      
      end
      
!-----------------------------------------------------------------------
!     subroutine
!========================================================================
!     read calculation condition
!-----------------------------------------------------------------------
      
      subroutine cal_cond
      
      use common_v
      implicit none
      
      m=750         !(x-direction mesh number)
      n=750         !(y-direction mesh number)
      dx=1.0e-8     !mesh size
      dy=1.0e-8     !mesh size
      dl=3.0e-9     !dl
      ds=3.0e-13    !ds
      co=0.0196     !initial solute content
      tmp=900       !initial temperature
      tmpmelt=933.3 !melting point
      vm=10.547e-6  !moler volume
      xme=640.0     !liquidus slope
      ke=0.14       !partition coefficient
      beta=0.0      !kinetic coefficient
      v=0.03        !anisotropy
      yk=4.0        !anisotropy
      sigma=0.093   !interface energy
      comnoise=0.0  !noise
      
      nc=n; mc=m; dl2=2.0*dl/dx; ds2=2.0*ds/dx
      cle=( tmpmelt - tmp)/xme; cse=cle*ke
      
      end
      
!-----------------------------------------------------------------------
!     init_cond
!-----------------------------------------------------------------------
      
      subroutine init_cond
      use common_v
      implicit none
      
      allocate( phi(0:1,0:m+1, 0:n+1) )
      allocate( com(0:1,0:m+1, 0:n+1) )
      
      do i=1, m
        do j=1, n
          
          phi(0,i,j)=0.0; phi(1,i,j)=0.0
          com(0,i,j)=co ; com(1,i,j)=co
          
        end do
      end do
      
      do i=1, triangle+5
        do j=1, triangle+5
          
          if( j .lt. (-i+triangle) ) then
            phi(0,i,j)=1.0; phi(1,i,j)=1.0
            com(0,i,j)=cse; com(1,i,j)=cse
            
          end if
          
          if( j .eq. (-i+triangle) ) then
            phi(0,i,j)=0.5; phi(1,i,j)=0.5
          end if
          
        end do
      end do
      
      do i=0, m+1
        phi(0,i,0)=phi(0,i,1); phi(0,i,n+1)=phi(0,i,n)
        com(0,i,0)=com(0,i,1); com(0,i,n+1)=com(0,i,n)
      end do
      
      do j=0, n+1
        phi(0,0,j)=phi(0,1,j); phi(0,m+1,j)=phi(0,m,j)
        com(0,0,j)=com(0,1,j); com(0,m+1,j)=com(0,m,j)
      end do
      
      return
      end
      
!-----------------------------------------------------------------------
!     pf mobility
!-----------------------------------------------------------------------
      subroutine mobility
      use common_v
      implicit none
      
      real(8)::p1,p2,pp1,pp2,fun1,fun2,zeta,fccle,fccse,alpha
      integer::ip1
      
      ep2=ep*ep; zeta=0.0
      
      fccle=r*tmp/vm/( cle*(1.0 - cle) )
      fccse=r*tmp/vm/( cse*(1.0 - cse) )
      
      do ip1=1, 998 !do p1=0.001, 0.998, 0.001 !Deleted feature: Loop variable at (1) must be integer
        p1=float(ip1)/1000.0
        p2=p1+0.001
        pp1=p1**3*(10.0 - 15.0*p1 + 6.0*p1*p1 )
        pp2=p2**3*(10.0 - 15.0*p2 + 6.0*p2*p2 )
        fun1=pp1*(1.0 - pp1)/( (1.0 - pp1)*fccse + pp1*fccle )/( p1*(1.0 - p1) )
        fun2=pp2*(1.0 - pp2)/( (1.0 - pp2)*fccse + pp1*fccle )/( p2*(1.0 - p2) )
        zeta=zeta+(fun1 + fun2)*0.001/2.0
        
      end do
      
      alpha=beta*r*tmp*(1.0 - ke)/(vm*xme)
      xm=1.0/( ep*ep/sigma*(alpha + ep/(dl*sqrt(2.0*w))*zeta*fccse*fccle*(cle-cse)**2) )
      
      
      dt=dx**2/(5.0*xm*ep**2)
      dt1=dx**2/(5.0*dl)
      dt=dmin1(dt,dt1)
      
      return
      end
      
!-----------------------------------------------------------------------
!     comnoise
!-----------------------------------------------------------------------
      
      subroutine noise
      use common_v
      implicit none
      real(8)::countnoise, comtot, comdam, cnoise
      integer::lam=12869,c=6925,mmu=32768,x=19724
      
      countnoise=0.0; comtot=0.0
      
      do i=1, mc
        do j=1, nc
          
          if( phi(0,i,j) .gt. 0.01 .and. phi(0,i,j) .le. 0.5 ) then
            
            comdam=com(0,i,j)
            
            x=mod( (x*lam+c), mmu )
            cnoise=( real(x)/mmu - 0.5 )*comnoise
            
            com(0,i,j)=comdam*( 1.0 + cnoise )
            comtot=comtot+comdam*cnoise
            countnoise=countnoise+1.0
            
          end if
          
        end do
      end do
      
      do i=1, mc
        do j=1, nc
          
          if( phi(0,i,j) .gt. 0.01 .and. phi(0,i,j) .le. 0.5 ) then
            
            com(0,i,j)=com(0,i,j)-(comtot/countnoise)
          
          end if
          
        end do
      end do
      
      do i=0, m+1
        
        com(0,i,0)=com(0,i,1); com(0,i,n+1)=com(0,i,n)
        
      end do
      
      do j=0, n+1
        
        com(0,0,j)=com(0,1,j); com(0,m+1,j)=com(0,m,j)
        
      end do
      
      return
      end
      
!-----------------------------------------------------------------------
!     area set
!-----------------------------------------------------------------------
      
      subroutine areaset
      use common_v
      implicit none
      
      do j=1, n
        if( abs(com(0,1,j)/co - 1.0 ) .gt. 1.0e-5 ) nc=j
      end do
      
      do i=1, m
        if( abs(com(0,i,1)/co - 1.0 ) .gt. 1.0e-5 ) mc=i
      end do
      
      nc=nc+10; mc=mc+10
      if( nc .gt. n ) nc=n; if( mc .gt. m ) mc=m
      
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
      
      write(14,410) mc,nc,l
      
      do i=1, mc
        do j=1, nc
          write(14,300) i,j,phi(0,i,j),com(0,i,j)
        end do
      end do
      
      close(14)
      
      return
      
  410 format(i4,i4,i10)
  300 format(i5,i5,e12.5,e12.5)
 1000 write(6,*) "Error in file open"
      
      end