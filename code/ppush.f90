subroutine ppush(n,ns)

  use gem_com
  use gem_equil
  implicit none
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,dum1,bstar
  INTEGER :: m,i,j,k,l,n,ns
  INTEGER :: np_old,np_new
  real :: rhog,vfac,kap,vpar,pidum,kaptp,kapnp,xnp
  real :: b,th,r,enerb,cost,sint,qr,laps,sz,ter
  real :: xt,xs,yt,xdot,ydot,zdot,pzdot,edot,pzd0,vp0
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
  real :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp,vncp,vparspp,psip2p,bdcrvbp,curvbzp,dipdrp
  integer :: mynopi
  

  pidum = 1./(pi*2)**1.5*vwidth**3
  mynopi = 0
  nopi(ns) = 0
ppush1_start_tm=ppush1_start_tm+MPI_WTIME()
!$acc parallel loop gang vector private(rhox,rhoy)
  do m=1,mm(ns)
     r=x2(m,ns)-0.5*lx+lr0
     k = int(z2(m,ns)/delz)
     wz0 = ((k+1)*delz-z2(m,ns))/delz
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

     curvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
          +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1)
     bdcrvbp = wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
          +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1)
     grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
          +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1)

     fp = wx0*f(i)+wx1*f(i+1)
     jfnp = wz0*jfn(k)+wz1*jfn(k+1)
     psipp = wx0*psip(i)+wx1*psip(i+1)
     psp = wx0*psi(i)+wx1*psi(i+1)
     ter = wx0*t0s(ns,i)+wx1*t0s(ns,i+1)
     kaptp = wx0*capts(ns,i)+wx1*capts(ns,i+1)
     kapnp = wx0*capns(ns,i)+wx1*capns(ns,i+1)
     xnp = wx0*xn0s(ns,i)+wx1*xn0s(ns,i+1)
     vncp = wx0*phincp(i)+wx1*phincp(i+1)
     vparspp = wx0*vparsp(ns,i)+wx1*vparsp(ns,i+1)
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)
     b=1.-tor+tor*bfldp
     pzp = mims(ns)*u2(m,ns)/b-q(ns)*psp/br0

     rhog=sqrt(2.*b*mu(m,ns)*mims(ns))/(q(ns)*b)*iflr

     rhox(1) = rhog*(1-tor)+rhog*grp*tor
     rhoy(1) = rhog*gxdgyp/grp*tor
     rhox(2) = -rhox(1)
     rhoy(2) = -rhoy(1)
     rhox(3) = 0
     rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
     rhox(4) = 0
     rhoy(4) = -rhoy(3)
     !    calculate avg. e-field...
     !    do 1,2,4 point average, where lr is the no. of points...

     phip=0.
     exp1=0.
     eyp=0.
     ezp=0.
     delbxp=0.
     delbyp=0.
     dpdzp = 0.
     dadzp = 0.
     aparp = 0.

     !  4 pt. avg. done explicitly for vectorization...
!$acc loop seq
     do l=1,lr(1)
        !
        xs=x2(m,ns)+rhox(l) !rwx(1,l)*rhog
        yt=y2(m,ns)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
        !
        !   particle can go out of bounds during gyroavg...
        xt=mod(xs+800.*lx,lx)
        yt=mod(yt+800.*ly,ly)
        xt = min(xt,lx-1.0e-8)
        yt = min(yt,ly-1.0e-8)

        include "ppushngp.h"
     enddo
     exp1 = exp1/4.
     eyp = eyp/4.
     ezp = ezp/4.
     delbxp = delbxp/4.
     delbyp = delbyp/4.
     dpdzp = dpdzp/4.
     dadzp = dadzp/4.
     aparp = aparp/4.
     !
     vfac = 0.5*(mims(ns)*u2(m,ns)**2 + 2.*mu(m,ns)*b)
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw

     vpar = u2(m,ns)-q(ns)/mims(ns)*aparp*nonlin(ns)*0.
     bstar = b*(1+mims(ns)*vpar/(q(ns)*b)*bdcrvbp)
     enerb=(mu(m,ns)+mims(ns)*vpar*vpar/b)/q(ns)*b/bstar*tor

     kap = kapnp - (1.5-vfac/ter)*kaptp-vpar*mims(ns)/ter*vparspp*vparsw
     dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
     vxdum = (eyp/b+vpar/b*delbxp)*dum1
     xdot = vxdum*nonlin(ns) -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
     ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlin(ns) &
          +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
          (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0   &
          +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
          -mims(ns)*vpar**2/(q(ns)*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
          -dipdrp/radiusp*mims(ns)*vpar**2/(q(ns)*bstar*b)*grcgtp*lr0/q0*qhatp
     zdot =  vpar*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
          +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*mims(ns)*vpar**2/(q(ns)*bstar*b)*q0*br0*grcgtp/jfnp

     pzd0 = tor*(-mu(m,ns)/mims(ns)/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
          +mu(m,ns)*vpar/(q(ns)*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
     pzdot = pzd0 + (q(ns)/mims(ns)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
          +q(ns)/mims(ns)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*ipara

     edot = q(ns)*(xdot*exp1+(ydot-vp0)*eyp+zdot*ezp)                      &
          +q(ns)*pzdot*aparp*tor     &
          +q(ns)*vpar*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)    &
          -q(ns)*vpar*delbxp*vp0

     x3(m,ns) = x2(m,ns) + 0.5*dt*xdot
     y3(m,ns) = y2(m,ns) + 0.5*dt*ydot
     z3(m,ns) = z2(m,ns) + 0.5*dt*zdot
     u3(m,ns) = u2(m,ns) + 0.5*dt*pzdot

     dum = 1-w2(m,ns)*nonlin(ns)*0.
     if(ildu.eq.1)dum = (tgis(ns)/ter)**1.5*exp(vfac*(1/tgis(ns)-1./ter))
     vxdum = (eyp/b+vpar/b*delbxp)*dum1
     !         vxdum = eyp+vpar/b*delbxp
     w3(m,ns)=w2(m,ns) + 0.5*dt*(vxdum*kap + edot/ter)*dum*xnp

     !         if(x3(m,ns)>lx .or. x3(m,ns)<0.)w3(m,ns) = 0.


     !if(itube/=1) then
     !if(abs(pzp-pzi(m,ns))>pzcrit(ns).or.abs(vfac-eki(m,ns))>0.2*eki(m,ns))then
     !   mynopi = mynopi+1
     !   x3(m,ns) = xii(m,ns)
     !   z3(m,ns) = z0i(m,ns)
     !   r = x3(m,ns)-lx/2+lr0
     !   k = int(z3(m,ns)/delz)
     !   wz0 = ((k+1)*delz-z3(m,ns))/delz
     !   wz1 = 1-wz0
     !   th = wz0*thfnz(k)+wz1*thfnz(k+1)

     !   i = int((r-rin)/dr)
     !   wx0 = (rin+(i+1)*dr-r)/dr
     !   wx1 = 1.-wx0
     !   k = int((th+pi)/dth)
     !   wz0 = (-pi+(k+1)*dth-th)/dth
     !   wz1 = 1.-wz0
     !   b = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
     !        +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
     !   u3(m,ns) = u0i(m,ns)
     !   u2(m,ns) = u3(m,ns)
     !   w3(m,ns) = 0.
     !   w2(m,ns) = 0.
     !   x2(m,ns) = x3(m,ns)
     !   z2(m,ns) = z3(m,ns)
     !end if
     !end if

     laps=anint((z3(m,ns)/lz)-.5)*(1-peritr)
     r=x3(m,ns)-0.5*lx+lr0
     i = int((r-rin)/dr)
    ! i = min(i,nr-1)
     i = max(i,0)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     y3(m,ns)=mod(y3(m,ns)-laps*2*pi*qr*lr0/q0*sign(1.0,q0)+8000.*ly,ly)
     if(x3(m,ns)>lx.and.iperidf==0)then
        x3(m,ns) = lx-1.e-8
        z3(m,ns)=lz-z3(m,ns)
        x2(m,ns) = x3(m,ns)
        z2(m,ns) = z3(m,ns)
        w2(m,ns) = 0.
        w3(m,ns) = 0.
     end if
     if(x3(m,ns)<0..and.iperidf==0)then
        x3(m,ns) = 1.e-8
        z3(m,ns)=lz-z3(m,ns)
        x2(m,ns) = x3(m,ns)
        z2(m,ns) = z3(m,ns)
        w2(m,ns) = 0.
        w3(m,ns) = 0.
     end if
     z3(m,ns)=mod(z3(m,ns)+8.*lz,lz)
     x3(m,ns)=mod(x3(m,ns)+800.*lx,lx)
     x3(m,ns) = min(x3(m,ns),lx-1.0e-8)
     y3(m,ns) = min(y3(m,ns),ly-1.0e-8)
     z3(m,ns) = min(z3(m,ns),lz-1.0e-8)

  enddo
!$acc end parallel
ppush1_end_tm=ppush1_end_tm+MPI_WTIME()
!!$acc update host(w3,x2,x3,y3,z2,z3,u3,w2)
  call MPI_ALLREDUCE(mynopi,nopi(ns),1,MPI_integer, &
       MPI_SUM, MPI_COMM_WORLD,ierr)

  np_old=mm(ns)
  !$acc update host(z3)
  call init_pmove(z3(:,ns),np_old,lz,ierr)

  call test_pmove(x2(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(x3(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(y2(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(y3(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(z2(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(z3(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(u2(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(u3(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(w2(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(w3(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(mu(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(xii(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(z0i(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(pzi(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(eki(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call test_pmove(u0i(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  !$acc update host(x2,x3,y2,y3,z2,z3,u2,u3,w2,w3)
  

  !!$acc update host(x2) async(1)
  !call pmove(x2(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(x2) 
  !!$acc update host(x3) async(1)
  !call pmove(x3(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(x3)
  !!$acc update host(y2) async(1)
  !call pmove(y2(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(y2)
  !!$acc update host(y3) async(1)
  !call pmove(y3(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(y3)
  !!$acc update host(z2) async(1)
  !call pmove(z2(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(z2)
  !call pmove(z3(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(z3)
  !!$acc update host(u2) async(1)
  !call pmove(u2(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(u2)
  !!$acc update host(u3) async(1)
  !call pmove(u3(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(u3)
  !!$acc update host(w2) async(1)
  !call pmove(w2(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(w2)
  !!$acc update host(w3) async(1)
  !call pmove(w3(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(w3)
  !call pmove(mu(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(mu)
  !call pmove(xii(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(xii)
  !call pmove(z0i(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(z0i)
  !call pmove(pzi(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(pzi)
  !call pmove(eki(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(eki)
  !call pmove(u0i(:,ns),np_old,np_new,ierr)
  !if (ierr.ne.0) call ppexit
  !!$acc update device(u0i)
  !!$acc wait(1)
  !!!$acc update device(w3,x2,x3,y2,y3,z2,z3,u2,u3,w2,mu,xii,z0i,pzi,eki,u0i)

  call end_pmove(ierr)
  mm(ns)=np_new

  !      return
end subroutine ppush

