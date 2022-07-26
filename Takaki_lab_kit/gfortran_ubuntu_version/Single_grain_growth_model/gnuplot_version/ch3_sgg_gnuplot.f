      implicit double precision (a-h,o-z)
      parameter(mlax=51,nlax=mlax+1)
      parameter(mlay=51,nlay=mlay+1)
      dimension p(0:nlax,0:nlay),pp(0:nlax,0:nlay)
      character*6 fname
c
      dx    = 0.5e-6
      dy    = dx
      gamma = 1.
      delta = 4.*dx
      pree  = 6.2e-6
      activ = 2.08e-19
      tempe = 800
      boltz = 1.38e-23
      amobi = pree*exp(-activ/(boltz*tempe))
      ram   = 0.1
      bbb=2.*log((1.+(1.-2.*ram))/(1.-(1.-2.*ram)))/2
c
      aaa   = sqrt(3.*delta*gamma/bbb)
      www   = 6.*gamma*bbb/delta
      pmobi = amobi*sqrt(2.*www)/(6.*aaa)
      beta  = 0.5
c
      dt    = dx*dx/(5.*pmobi*aaa*aaa)
      write(*,*) 'Time increment =',dt,'[s]'
      nend  = 100
c
      r0=real(mlax/2)*dx
      do m=1,mlay
      do l=1,mlax
         x=real(l-1)*dx
         y=real(m-1)*dy
         r=sqrt(x*x+y*y)-r0
         p(l,m)=0.5*(1.-tanh(sqrt(2.*www)/(2.*aaa)*r))
      enddo
      enddo
c
      do m=0,nlay
         p(0,m)    =p(1,m)
         p(nlax,m)=p(mlax,m)
      enddo
      do l=0,nlax
         p(l,0)    =p(l,1)
         p(l,nlay)=p(l,mlay)
      enddo
c
      do nstep=1,nend
      write(*,*) 'Step No. =',nstep
c
      do m=1,mlay
      do l=1,mlax
         rpx =(p(l-1,m)-2.*p(l,m)+p(l+1,m))/(dx*dx)
         rpy =(p(l,m-1)-2.*p(l,m)+p(l,m+1))/(dy*dy)
         dpi1=aaa*aaa*(rpx+rpy)
         dpi2=4.*www*p(l,m)*(1.-p(l,m))*(p(l,m)-0.5+beta)
         dpi =dpi1+dpi2
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
         p(0,m)    =p(1,m)
         p(nlax,m)=p(mlax,m)
      enddo
      do l=0,nlax
         p(l,0)    =p(l,1)
         p(l,nlay)=p(l,mlay)
      enddo
c
c^^^^^^^可視化ソフトVisible2D用の出力：noutステップ毎に結果を出力^^^^^^^^^^
      nout=10
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
      write(100,*) ""
      do l=1,mlax
         write(100,'(3e16.7)') real(l-1)*dx,real(m-1)*dy,p(l,m)
      enddo
      enddo
      close(100)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
      endif
c
      enddo
c
      stop
      end
