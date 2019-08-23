subroutine grid1(ip,n)

  !    source quantities are are calculated: n_i
  !    right now only ion quantitities are calculated...

  use gem_com
  use gem_equil
  implicit none
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: enerb,vxdum,dum,xdot,ydot
  INTEGER :: m,n,i,j,k,l,ns,ip
  real :: wx0,wx1,wy0,wy1,wz0,wz1,vte
  real :: sz,wght,r,th,cost,sint,b,qr,dv
  real :: xt,yt,rhog,pidum,vpar,xs,dely,vfac
  real :: lbfs(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  real :: lbfr(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: myden(0:imx,0:jmx,0:1),myjpar(0:imx,0:jmx,0:1)
  real :: mydene(0:imx,0:jmx,0:1),myupar(0:imx,0:jmx,0:1)
  real :: mydti(0:imx,0:jmx,0:1),mydte(0:imx,0:jmx,0:1)
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
  real :: grp,gxdgyp,rhox(4),rhoy(4)


  rho=0.
  jion = 0.
  mydte = 0.
  ns=1
  if(idg.eq.1)write(*,*)'enter ion grid1',mm(1)
  do ns = 1,nsm
     den(ns,:,:,:)=0.
     jpar(ns,:,:,:)=0.
     myden = 0.
     myjpar = 0.
     mydti = 0.
grid11_start_tm=grid11_start_tm+MPI_WTIME()
!$acc parallel loop gang vector private(rhox,rhoy)
     do m=1,mm(ns)
        dv=float(lr(1))*(dx*dy*dz)
        r=x3(ns,m)-0.5*lx+lr0

        k = int(z3(ns,m)/delz)
        wz0 = ((k+1)*delz-z3(ns,m))/delz
        wz1 = 1-wz0
        th = wz0*thfnz(k)+wz1*thfnz(k+1)

        i = int((r-rin)/dr)
        wx0 = (rin+(i+1)*dr-r)/dr
        wx1 = 1.-wx0
        k = int((th+pi)/dth)
        wz0 = (-pi+(k+1)*dth-th)/dth
        wz1 = 1.-wz0
        dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
             +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1)
        dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
             +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1)
        grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
             +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1)
        bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
             +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
        radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
             +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1)
        dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
             +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1)
        qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
             +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1)
        grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
             +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1)
        gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
             +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1)
        fp = wx0*f(i)+wx1*f(i+1)
        jfnp = wz0*jfn(k)+wz1*jfn(k+1)
        psipp = wx0*psip(i)+wx1*psip(i+1)
        !         b=1.-lr0/br0*cost
        b=1.-tor+tor*bfldp

        rhog=sqrt(2.*b*mu(ns,m)*mims(ns))/(q(ns)*b)*iflr

        rhox(1) = rhog*(1-tor)+rhog*grp*tor
        rhoy(1) = rhog*gxdgyp/grp*tor
        rhox(2) = -rhox(1)
        rhoy(2) = -rhoy(1)
        rhox(3) = 0
        rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
        rhox(4) = 0
        rhoy(4) = -rhoy(3)

        vfac=0.5*(mims(ns)*u3(ns,m)**2 + 2.*mu(ns,m)*b )
        wght=w3(ns,m)/dv

        vpar = u3(ns,m) !linearly correct

        !    now do 1,2,4 point average, where lr is the no. of points...
!$acc loop seq
        do l=1,lr(1)
           xs=x3(ns,m)+rhox(l) !rwx(1,l)*rhog
           yt=y3(ns,m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
           xt=mod(xs+800.*lx,lx)
           yt=mod(yt+800.*ly,ly)
           xt = min(xt,lx-1.e-8)
           yt = min(yt,ly-1.e-8)

           i=int(xt/dx+0.5)
           j=int(yt/dy+0.5)
           k=int(z3(ns,m)/dz+0.5)-gclr*kcnt


!$acc atomic update
           myden(i,j,k) = myden(i,j,k) + wght
!$acc end atomic
!$acc atomic update
           myjpar(i,j,k) = myjpar(i,j,k)+wght*vpar
!$acc end atomic
!$acc atomic update
           mydti(i,j,k) = mydti(i,j,k)+wght*vfac
!$acc end atomic

        enddo
     enddo
!$acc end parallel
!$acc wait
grid11_end_tm=grid11_end_tm+MPI_WTIME()

     if(idg.eq.1)write(*,*)myid,'pass ion grid1'
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     !   enforce periodicity
     call enforce(myden(:,:,:))
     call enforce(myjpar)
     call enforce(mydti)
     !      call filter(myden(:,:,:))
     !      call filter(myjpar(:,:,:))

     do i=0,im
        do j=0,jm
           do k=0,mykm
              den(ns,i,j,k)=q(ns)*myden(i,j,k)/n0/jac(i,k)*cn0s(ns)
              jpar(ns,i,j,k) = q(ns)*myjpar(i,j,k)/n0/jac(i,k)*cn0s(ns)
              mydti(i,j,k) = mydti(i,j,k)/n0/jac(i,k)*cn0s(ns)
           enddo
        enddo
     enddo

     do i = 0,im
        do j = 0,jm
           do k = 0,1
              mydti(i,j,k) = (mydti(i,j,k)-gt0i(i)*den(ns,i,j,k))/gn0s(1,i)
           end do
        end do
     end do
     call MPI_ALLREDUCE(den(ns,0:im,0:jm,0:1),  &
          dti(ns,0:im,0:jm,0:1),             &
          (imx+1)*(jmx+1)*2,MPI_REAL8,       &
          MPI_SUM,GRID_COMM,ierr)


     do i=0,im
        do j=0,jm
           do k=0,mykm
              rho(i,j,k)=rho(i,j,k)+den(ns,i,j,k)
              jion(i,j,k) = jion(i,j,k)+jpar(ns,i,j,k)
           enddo
        enddo
     enddo
  end do

end subroutine grid1


