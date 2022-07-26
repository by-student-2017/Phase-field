      implicit double precision (a-h,o-z)
c
      parameter(npo=4)
      parameter(mlax=51,nlax=mlax+1)
      parameter(mlay=51,nlay=mlay+1)
      dimension p(npo,0:nlax,0:nlay),pp(npo,0:nlax,0:nlay)
      dimension aij(npo,npo),wij(npo,npo),tij(npo,npo),eij(npo,npo)
      dimension mfij(npo,mlax,mlay),nfij(mlax,mlay)
      dimension ls(npo),ms(npo)
      character*6 fname
c
      dx    = 0.5e-6
      dy    = dx
      gamma = 1.
      delta = 7.*dx
      pree  = 6.2e-6
      activ = 2.08e-19
      tempe = 800
      boltz = 1.38e-23
      amobi = pree*exp(-activ/(boltz*tempe))
      eee   = 1.e+6
      pi    = 3.141592654
c
      aaa   = 2./pi*sqrt(2.*delta*gamma)
	www   = 4.*gamma/delta
      pmobi = amobi*pi*pi/(8.*delta)
c
      dt    = dx*dx/(5.*pmobi*aaa*aaa)
      write(*,*) 'Time increment =',dt,'[s]'
      nend  = 500
c
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
c
c
      ls(1)=10
      ms(1)=10
      ls(2)=40
      ms(2)=20
      ls(3)=25
      ms(3)=40
c
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
c
      do m=1,mlay
      do l=1,mlax
         do i=1,npo
            pp(i,l,m)=0.
         enddo
      enddo
      enddo
c
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
c
      nstep=0
      do nstep=1,nend
      write(*,*) 'Step No. =',nstep
c
      do l=1,mlax
      do m=1,mlay
         n=0
         do i=1,npo
            if(p(i,l  ,m  ).gt.0..or.
     &         (p(i,l  ,m  ).eq.0..and.p(i,l+1,m  ).gt.0..or.
     &                                 p(i,l-1,m  ).gt.0..or.
     &                                 p(i,l  ,m+1).gt.0..or.
     &                                 p(i,l  ,m-1).gt.0.)   ) then
            n=n+1
            mfij(n,l,m)=i
            endif
         enddo
         nfij(l,m)=n
      enddo
      enddo
c
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
            ppp=ppp+
     &        (aij(i,k)*aij(i,k)-aij(j,k)*aij(j,k))/2.
     &        *((p(k,l+1,m)-2.*p(k,l,m)+p(k,l-1,m))/(dx*dx)
     &         +(p(k,l,m+1)-2.*p(k,l,m)+p(k,l,m-1))/(dy*dy))
     &         +(wij(i,k)-wij(j,k))*p(k,l,m)
         enddo
            pee=p(i,l,m)*p(j,l,m)
            dpi=dpi-2.*tij(i,j)/real(nfij(l,m))
     &              *(ppp-8./pi*sqrt(pee)*eij(i,j))
         enddo
            pp(i,l,m)=p(i,l,m)+dpi*dt
         enddo
      enddo
      enddo
c
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
c
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
c
c
c^^^^^^^可視化ソフトVisible2D用の出力：noutステップ毎に結果を出力^^^^^^^^^^
      nout=50
      if(nstep.eq.nstep/nout*nout) then
c
      i100000=nstep/100000
      i10000=(nstep-i100000*100000)/10000
      i1000=(nstep-i100000*100000-i10000*10000)/1000
      i100=(nstep-i100000*100000-i10000*10000-i1000*1000)/100
      i10=(nstep-i100000*100000-i10000*10000-i1000*1000-i100*100)/10
      i1=(nstep-i100000*100000-i10000*10000-i1000*1000-i100*100-i10*10)
      fname=char(48+i100000)//char(48+i10000)//char(48+i1000)
     &      //char(48+i100)//char(48+i10)//char(48+i1)
c
      open(100,file='out'//fname//'.dat')
      write(100,*)
      write(100,*)
      write(100,'(2i10)') mlax-1,mlay-1
      write(100,'(A)') 'x y gb phase-4'
      do m=1,mlay
      do l=1,mlax
         pgb=0.
         do i=1,npo
            pgb=pgb+p(i,l,m)*p(i,l,m)
         enddo
         write(100,'(4e16.7)') real(l-1)*dx,real(m-1)*dy,pgb,p(npo,l,m)
      enddo
      enddo
      close(100)
      endif
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
      enddo
c
      stop
	end
