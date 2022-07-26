      implicit double precision (a-h,o-z)
      parameter(mlax=101,nlax=mlax+1)
      parameter(mlay=101,nlay=mlay+1)
      dimension p(0:nlax,0:nlay),pp(0:nlax,0:nlay)
      dimension t(0:nlax,0:nlay),tt(0:nlax,0:nlay)
      character*6 fname
c
      dx    = 20.e-9     !x方向格子サイズ[m]
      dy    = dx         !x方向格子サイズ[m]
      gamma = 0.37       !界面エネルギー[J/m^2]
      delta = 4.*dx      !界面幅[m]
      ram   = 0.1        !界面入り口のφの値
      bbb   = 2.*log((1.+(1.-2.*ram))/(1.-(1.-2.*ram)))/2
c
c    << pure Ni >>
      cndct = 84.01      !熱伝導率[W/mK]
      speht = 5.42e+06   !比熱[J/Km^3]
      rlate = 2.350e+09  !潜熱[J/m^3]
      tmelt = 1728.      !融点[K]
      skine = 2.         !界面カイネティック係数[m/Ks]
      t_ini = 1511.2     !初期温度[K]
c
      aaa   = sqrt(3.*delta*gamma/bbb)
      www   = 6.*gamma*bbb/delta
      pmobi = bbb*tmelt*skine/(3.*delta*rlate)
c
      dtp   = dx*dx/(5.*pmobi*aaa*aaa)
      dtt   = dx*dx/(5.*cndct/speht)
      dt    = min(dtp,dtt)
      write(*,*) 'Time increment =',dt,'[s]'
      nend  = 1000
c
      r0=2.*dx
      do m=1,mlay
      do l=1,mlax
         x=real(l-1)*dx-real(mlax/2)*dx
         y=real(m-1)*dy-real(mlay/2)*dy
         r=sqrt(x*x+y*y)-r0
         p(l,m)=0.5*(1.-tanh(sqrt(2.*www)/(2.*aaa)*r))
         t(l,m)=t_ini
      enddo
      enddo
c
      do m=0,nlay
         p(0,m)   =p(mlax,m)
         p(nlax,m)=p(1,m)
         t(0,m)   =t(mlax,m)
         t(nlax,m)=t(1,m)
      enddo
      do l=0,nlax
         p(l,0)   =p(l,mlay)
         p(l,nlay)=p(l,1)
         t(l,0)   =t(l,mlay)
         t(l,nlay)=t(l,1)
      enddo
c
      do nstep=1,nend
      write(*,*) 'Step No. =',nstep
c
      do m=1,mlay
      do l=1,mlax
c
        bet =-15./(2.*www)*rlate*(t(l,m)-tmelt)/tmelt*p(l,m)*(1.-p(l,m))
        rpx =(p(l-1,m)-2.*p(l,m)+p(l+1,m))/(dx*dx)
        rpy =(p(l,m-1)-2.*p(l,m)+p(l,m+1))/(dy*dy)
        dpi1=aaa*aaa*(rpx+rpy)
        dpi2=4.*www*p(l,m)*(1.-p(l,m))*(p(l,m)-0.5+bet)
        dpi =dpi1+dpi2
        dpt =pmobi*dpi
        pp(l,m)=p(l,m)+dpt*dt
c
        rtx=(t(l-1,m  )-2.*t(l,m)+t(l+1,m  ))/(dx*dx)
        rty=(t(l  ,m-1)-2.*t(l,m)+t(l  ,m+1))/(dy*dy)
        pd =30.*p(l,m)*p(l,m)*(1.-p(l,m))*(1.-p(l,m))
        dtt =(cndct*(rtx+rty)+rlate*pd*dpt)/speht
        tt(l,m)=t(l,m)+dtt*dt
c
      enddo
      enddo
c
      do m=1,mlay
      do l=1,mlax
         p(l,m)=pp(l,m)
         t(l,m)=tt(l,m)
      enddo
      enddo
      do m=0,nlay
         p(0,m)   =p(mlax,m)
         p(nlax,m)=p(1,m)
         t(0,m)   =t(mlax,m)
         t(nlax,m)=t(1,m)
      enddo
      do l=0,nlax
         p(l,0)   =p(l,mlay)
         p(l,nlay)=p(l,1)
         t(l,0)   =t(l,mlay)
         t(l,nlay)=t(l,1)
      enddo
c
c^^^^^^^可視化ソフトVisible2D用の出力：noutステップ毎に結果を出力^^^^^^^^^^
      nout=100
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
      write(100,'(A)') 'x y phase temperature'
      do m=1,mlay
      do l=1,mlax
         write(100,'(4e16.7)') real(l-1)*dx,real(m-1)*dy,p(l,m),t(l,m)
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
