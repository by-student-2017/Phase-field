      implicit double precision (a-h,o-z)
      parameter(mlax=101,nlax=mlax+1)
      parameter(mlay=101,nlay=mlay+1)
      dimension p(0:nlax,0:nlay),pp(0:nlax,0:nlay),ppp(0:nlax,0:nlay)
      character*6 fname
c
      dx    = 0.5e-6
      dy    = dx
      gamma = 1.
      delta = 4.*dx
      ram   = 0.1
      bbb=2.*log((1.+(1.-2.*ram))/(1.-(1.-2.*ram)))/2
c
      aaa   = sqrt(3.*delta*gamma/bbb)
      www   = 6.*gamma*bbb/delta
      pmobi = 1.
      beta  = 0.
c
      dt    = dx*dx*dx*dx/(40.*pmobi*aaa*aaa)
      write(*,*) 'Time increment =',dt,'[s]'
      nend  = 100000
c
      r0=real(mlax/4)*dx
      l0=mlax/2
      m0=mlay/2
      do m=1,mlay
      do l=1,mlax
         ll=l-l0
         mm=m-m0
         if(ll.ge.0.and.ll.ge.abs(mm)) then
          r0=3./4.*dx*real(mlax-1)
          p(l,m)=0.5*(1.-tanh(sqrt(2.*www)/(2.*aaa)*( dx*real(l-1)-r0)))
         endif
         if(ll.le.0.and.abs(ll).ge.abs(mm)) then
          r0=1./4.*dx*real(mlax-1)
          p(l,m)=0.5*(1.-tanh(sqrt(2.*www)/(2.*aaa)*(-dx*real(l  )+r0)))
         endif
         if(mm.ge.0.and.mm.ge.abs(ll)) then
          r0=3./4.*dx*real(mlay-1)
          p(l,m)=0.5*(1.-tanh(sqrt(2.*www)/(2.*aaa)*( dy*real(m-1)-r0)))
         endif
         if(mm.le.0.and.abs(mm).ge.abs(ll)) then
          r0=1./4.*dx*real(mlay-1)
          p(l,m)=0.5*(1.-tanh(sqrt(2.*www)/(2.*aaa)*(-dy*real(m  )+r0)))
         endif
      enddo
      enddo
c
      do m=0,nlay
         p(0,m)   =p(mlax,m)
         p(nlax,m)=p(1,m)
      enddo
      do l=0,nlax
         p(l,0)   =p(l,mlay)
         p(l,nlay)=p(l,1)
      enddo
c
      do nstep=1,nend
      write(*,*) 'Step No. =',nstep
c
      do m=1,mlay
      do l=1,mlax
         dpx=(p(l+1,m  )-p(l-1,m  ))/(2.*dx)
         dpy=(p(l  ,m+1)-p(l  ,m+1))/(2.*dy)
         df =www*p(l,m)*p(l,m)*(1.-p(l,m))*(1.-p(l,m))
     &      +aaa*aaa/2.*(dpx*dpx+dpy*dpy)
         rpx =(p(l-1,m)-2.*p(l,m)+p(l+1,m))/(dx*dx)
         rpy =(p(l,m-1)-2.*p(l,m)+p(l,m+1))/(dy*dy)
         dpi1=-aaa*aaa*(rpx+rpy)
         dpi2=-4.*www*p(l,m)*(1.-p(l,m))*(p(l,m)-0.5+beta)
         dpi =dpi1+dpi2
         ppp(l,m)=dpi
      enddo
      enddo
c
      do m=0,nlay
         ppp(0,m)   =ppp(mlax,m)
         ppp(nlax,m)=ppp(1,m)
      enddo
      do l=0,nlax
         ppp(l,0)   =ppp(l,mlay)
         ppp(l,nlay)=ppp(l,1)
      enddo
c
c
      do m=1,mlay
      do l=1,mlax
         rpx =(ppp(l-1,m)-2.*ppp(l,m)+ppp(l+1,m))/(dx*dx)
         rpy =(ppp(l,m-1)-2.*ppp(l,m)+ppp(l,m+1))/(dy*dy)
         dpi =rpx+rpy
         pp(l,m)=p(l,m)+pmobi*dpi*dt
      enddo
      enddo
c
      do m=1,mlay
      do l=1,mlax
         p(l,m)=pp(l,m)
      enddo
      enddo
      do m=0,nlay
         p(0,m)   =p(mlax,m)
         p(nlax,m)=p(1,m)
      enddo
      do l=0,nlax
         p(l,0)   =p(l,mlay)
         p(l,nlay)=p(l,1)
      enddo
c
c^^^^^^^可視化ソフトVisible2D用の出力：noutステップ毎に結果を出力^^^^^^^^^^
      nout=10000
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
      write(100,'(A)') 'x y phase'
      do m=1,mlay
      do l=1,mlax
         write(100,'(3e16.7)') real(l-1)*dx,real(m-1)*dy,p(l,m)
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
