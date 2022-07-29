      implicit double precision (a-h,o-z)
!-----------------------------------------------------------------------
      !parameter(npo=4)
      !parameter(mlax=51,nlax=mlax+1)
      !parameter(mlay=51,nlay=mlay+1)
      allocatable p(:,:,:), pp(:,:,:)
      !dimension p(npo,0:nlax,0:nlay),pp(npo,0:nlax,0:nlay)
      allocatable aij(:,:), wij(:,:), tij(:,:), eij(:,:)
      !dimension aij(npo,npo),wij(npo,npo),tij(npo,npo),eij(npo,npo)
      allocatable mfij(:,:,:), nfij(:,:)
      !dimension mfij(npo,mlax,mlay),nfij(mlax,mlay)
      allocatable ls(:), ms(:)
      !dimension ls(npo),ms(npo)
      character*6 fname
      character(len=20) filename
      real data(20)
      integer iout
      iout  = -1
!-----------------------------------------------------------------------
      open(1,file="parameters.txt")
      read(1,'(11x,e10.3)') (data(i),i=1,11)
      close(1)
!
      npo   = int(data(1))
      dx    = data(2)
      dy    = data(3)
      mlax  = int(data(4))
      nlax=mlax+1
      mlay  = int(data(5))
      nlay=mlay+1
      !
      allocate (  p(npo,0:nlax,0:nlay) )
      allocate ( pp(npo,0:nlax,0:nlay) )
      !
      allocate (aij(npo,npo))
      allocate (wij(npo,npo))
      allocate (tij(npo,npo))
      allocate (eij(npo,npo))
      !
      allocate ( mfij(npo,mlax,mlay) )
      allocate ( nfij(mlax,mlay) )
      !
      allocate (ls(npo))
      allocate (ms(npo))
!
      gamma = 1.
      delta = 7.*dx
      pree  = data(6)
      activ = data(7)
      tempe = data(8)
      boltz = 1.38e-23
      amobi = pree*exp(-activ/(boltz*tempe))
      eee   = data(9)
      pi    = 3.141592654
!
      aaa   = 2./pi*sqrt(2.*delta*gamma)
      www   = 4.*gamma/delta
      pmobi = amobi*pi*pi/(8.*delta)
!
      dt    = dx*dx/(5.*pmobi*aaa*aaa)
      write(*,*) 'Time increment =',dt,'[s]'
      nend  = int(data(10))
      nout  = int(data(11))
!
      do i=1,npo
      do j=1,npo
         wij(i,j)=www
         aij(i,j)=aaa
         tij(i,j)=pmobi
         eij(i,j)=0.
         if(i.eq.j) then
            wij(i,j)=0.
            aij(i,j)=0.
            tij(i,j)=0.
         endif
         if(i.eq.npo.or.j.eq.npo) eij(i,j)=eee
         if(i.gt.j) eij(i,j)=-eij(i,j)
      enddo
      enddo
!
!
      ls(1)=10
      ms(1)=10
      ls(2)=40
      ms(2)=20
      ls(3)=25
      ms(3)=40
!
      r0=3.*dx
      do m=1,mlay
      do l=1,mlax
         phi=0.
         do i=1,npo-1
            r=sqrt((real(l-ls(i))*dx)**2+(real(m-ms(i))*dy)**2)
            r=r-r0
            a=sqrt(2.*www)/aaa*r
            p(i,l,m)=0.5*(1.-sin(a))
            if(a.ge. pi/2.) p(i,l,m)=0.
            if(a.le.-pi/2.) p(i,l,m)=1.
            phi=phi+p(i,l,m)
         enddo
         p(npo,l,m)=1.-phi
      enddo
      enddo
!
      do m=1,mlay
      do l=1,mlax
         do i=1,npo
            pp(i,l,m)=0.
         enddo
      enddo
      enddo
!
      do m=0,nlay
      do i=1,npo
         p(i,0,m)   =p(i,mlax,m)
         p(i,nlax,m)=p(i,1,m)
      enddo
      enddo
      do l=0,nlax
      do i=1,npo
         p(i,l,0)   =p(i,l,mlay)
         p(i,l,nlay)=p(i,l,1)
      enddo
      enddo
!
      nstep=0
      do nstep=1,nend
      write(*,*) 'Step No. =',nstep
!
      do l=1,mlax
      do m=1,mlay
         n=0
         do i=1,npo
            if(p(i,l  ,m  ).gt.0..or.                         &
     &         (p(i,l  ,m  ).eq.0..and.p(i,l+1,m  ).gt.0..or. &
     &                                 p(i,l-1,m  ).gt.0..or. &
     &                                 p(i,l  ,m+1).gt.0..or. &
     &                                 p(i,l  ,m-1).gt.0.)   ) then
            n=n+1
            mfij(n,l,m)=i
            endif
         enddo
         nfij(l,m)=n
      enddo
      enddo
!
      do m=1,mlay
      do l=1,mlax
         do n1=1,nfij(l,m)
            i=mfij(n1,l,m)
            dpi=0.
         do n2=1,nfij(l,m)
            j=mfij(n2,l,m)
            ppp=0.
         do n3=1,nfij(l,m)
            k=mfij(n3,l,m)
            ppp=ppp+ &
     &        (aij(i,k)*aij(i,k)-aij(j,k)*aij(j,k))/2.       &
     &        *((p(k,l+1,m)-2.*p(k,l,m)+p(k,l-1,m))/(dx*dx)  &
     &         +(p(k,l,m+1)-2.*p(k,l,m)+p(k,l,m-1))/(dy*dy)) &
     &         +(wij(i,k)-wij(j,k))*p(k,l,m)
         enddo
            pee=p(i,l,m)*p(j,l,m)
            dpi=dpi-2.*tij(i,j)/real(nfij(l,m))     &
     &              *(ppp-8./pi*sqrt(pee)*eij(i,j))
         enddo
            pp(i,l,m)=p(i,l,m)+dpi*dt
         enddo
      enddo
      enddo
!
      do m=1,mlay
      do l=1,mlax
         phi=0.
         do i=1,npo
            if(pp(i,l,m).le.0.) pp(i,l,m)=0.
            if(pp(i,l,m).ge.1.) pp(i,l,m)=1.
            phi=phi+pp(i,l,m)
         enddo
         a=1.
         if(phi.ne.1.) a=1./phi
         do i=1,npo
            p(i,l,m) =a*pp(i,l,m)
            pp(i,l,m)=0.
         enddo
      enddo
      enddo
!
      do m=0,nlay
      do i=1,npo
         p(i,0,m)   =p(i,mlax,m)
         p(i,nlax,m)=p(i,1,m)
      enddo
      enddo
      do l=0,nlax
      do i=1,npo
         p(i,l,0)   =p(i,l,mlay)
         p(i,l,nlay)=p(i,l,1)
      enddo
      enddo
!
!
!-----------------------------------------------------------------------
      nout=50
      if(nstep.eq.nstep/nout*nout) then
      iout = iout + 1

      write(filename,'(a,i3.3,a)') 'mpf_result',iout,'.vtk'
!      open(100,file='out'//iout//'.vtk')
      open(100,file=filename)
      write(100,'(a)') '# vtk DataFile Version 3.0'
      write(100,'(a)') 'output.vtk'
      write(100,'(a)') 'ASCII'
      write(100,'(a)') 'DATASET STRUCTURED_POINTS'
      write(100,'(a,3i5)') 'DIMENSIONS', mlax,mlay,1
      write(100,'(a,3f4.1)') 'ORIGIN ', 0.0, 0.0, 0.0
      write(100,'(a,3i2)') 'ASPECT_RATIO', 1, 1, 1
      write(100,'(a,1i11)') 'POINT_DATA', mlax*mlay*1
      write(100,'(a)') 'SCALARS phase_field float'
      write(100,'(a)') 'LOOKUP_TABLE default'
      do m=1,mlay
      do l=1,mlax
         pgb=0.
         do i=1,npo
            pgb=pgb+p(i,l,m)*p(i,l,m)
         enddo
         write(100,'(1f10.6)') pgb
!         write(100,'(4e16.7)') real(l-1)*dx,real(m-1)*dy,pgb,p(npo,l,m)
      enddo
      enddo
      !
      write(100,'(a)') 'SCALARS phase_field_1 float'
      write(100,'(a)') 'LOOKUP_TABLE default'
      do m=1,mlay
      do l=1,mlax
         write(100,'(1f10.6)') p(1,l,m)
      enddo
      enddo
      !
      write(100,'(a)') 'SCALARS phase_field_2 float'
      write(100,'(a)') 'LOOKUP_TABLE default'
      do m=1,mlay
      do l=1,mlax
         write(100,'(1f10.6)') p(2,l,m)
      enddo
      enddo
      !
      write(100,'(a)') 'SCALARS phase_field_3 float'
      write(100,'(a)') 'LOOKUP_TABLE default'
      do m=1,mlay
      do l=1,mlax
         write(100,'(1f10.6)') p(3,l,m)
      enddo
      enddo
      !
      write(100,'(a)') 'SCALARS phase_field_4 float'
      write(100,'(a)') 'LOOKUP_TABLE default'
      do m=1,mlay
      do l=1,mlax
         write(100,'(1f10.6)') p(4,l,m)
      enddo
      enddo
      !
      close(100)
      endif
!-----------------------------------------------------------------------
!
      enddo
!
      deallocate (p)
      deallocate (pp)
      !
      deallocate (aij)
      deallocate (wij)
      deallocate (tij)
      deallocate (eij)
      !
      deallocate (mfij)
      deallocate (nfij)
      !
      deallocate (ls)
      deallocate (ms)
!
      stop
      end
