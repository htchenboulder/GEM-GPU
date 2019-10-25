subroutine cpush(n,ns)

  use gem_com
  use gem_equil
  implicit none
  INTEGER :: n
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: wx0,wx1,wy0,wy1,wz0,wz1,w3old
  INTEGER :: m,i,j,k,ns,l
  INTEGER :: np_old,np_new
  real :: vfac,vpar,vxdum,dum,xdot,ydot,zdot,pzdot,edot,pzd0,vp0
  real :: rhog,xt,yt,zt,kap,xs,pidum,dum1,kaptp,kapnp,xnp
  real :: b,th,r,enerb,cost,sint,qr,laps,sz,ter,bstar
  real :: myke,mypfl_es(1:nsubd),mypfl_em(1:nsubd),myavewi
  real :: myefl_es(1:nsubd),myefl_em(1:nsubd),mynos
  real :: ketemp,pfltemp
  real :: efltemp,nostemp
  real :: sbuf(10),rbuf(10)
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
  real :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp,vncp,vparspp,psip2p,bdcrvbp,curvbzp,dipdrp

  sbuf(1:10) = 0.
  rbuf(1:10) = 0.
  myavewi = 0.
  myke=0.
  mypfl_es=0.
  mypfl_em=0.
  myefl_es=0.
  myefl_em=0.
  mynos=0.
  ketemp=0.
  pfltemp=0.
  efltemp=0.
  nostemp=0.
  pidum = 1./(pi*2)**1.5*vwidth**3
cpush1_start_tm=cpush1_start_tm+MPI_WTIME()
!$acc parallel loop gang vector private(rhox,rhoy)
  do m=1,mm(ns)
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
     b=1.-tor+tor*bfldp
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)
     pzp = mims(ns)*u3(ns,m)/b-q(ns)*psp/br0
     vncp = wx0*phincp(i)+wx1*phincp(i+1)
     vparspp = wx0*vparsp(ns,i)+wx1*vparsp(ns,i+1)

     rhog=sqrt(2.*b*mu(ns,m)*mims(ns))/(q(ns)*b)*iflr

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
     delbxp = 0.
     delbyp = 0.
     dpdzp = 0.
     dadzp = 0.
     aparp = 0.

     !  4 pt. avg. written out explicitly for vectorization...
!$acc loop seq
     do l=1,lr(1)
        xs=x3(ns,m)+rhox(l) !rwx(1,l)*rhog
        yt=y3(ns,m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
        !   BOUNDARY
        xt=mod(xs+800.*lx,lx)
        yt=mod(yt+800.*ly,ly)
        xt = min(xt,lx-1.0e-8)
        yt = min(yt,ly-1.0e-8)

        include "cpushngp.h"
     enddo

     exp1=exp1/4.
     eyp=eyp/4.
     ezp=ezp/4.
     delbxp=delbxp/4.
     delbyp=delbyp/4.
     dpdzp = dpdzp/4.
     dadzp = dadzp/4.
     aparp = aparp/4.


     vfac = 0.5*(mims(ns)*u3(ns,m)**2 + 2.*mu(ns,m)*b)
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw

     vpar = u3(ns,m)-q(ns)/mims(ns)*aparp*nonlin(ns)*0.
     bstar = b*(1+mims(ns)*vpar/(q(ns)*b)*bdcrvbp)
     enerb=(mu(ns,m)+mims(ns)*vpar*vpar/b)/q(ns)*b/bstar*tor

     kap = kapnp - (1.5-vfac/ter)*kaptp-vpar*mims(ns)/ter*vparspp*vparsw
     dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
     vxdum = (eyp/b+vpar/b*delbxp)*dum1
     xdot = vxdum*nonlin(ns) -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
     ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlin(ns)     &
          +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
          (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
          +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
          -mims(ns)*vpar**2/(q(ns)*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
          -dipdrp/radiusp*mims(ns)*vpar**2/(q(ns)*bstar*b)*grcgtp*lr0/q0*qhatp
     zdot =  vpar*b/bstar*(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
          +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*mims(ns)*vpar**2/(q(ns)*bstar*b)*q0*br0*grcgtp/jfnp

     pzd0 = tor*(-mu(ns,m)/mims(ns)/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
          +mu(ns,m)*vpar/(q(ns)*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
     pzdot = pzd0 + (q(ns)/mims(ns)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
          +q(ns)/mims(ns)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*ipara

     edot = q(ns)*(xdot*exp1+(ydot-vp0)*eyp+zdot*ezp)                      &
          +q(ns)*pzdot*aparp*tor     &
          +q(ns)*vpar*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)   &
          -q(ns)*vpar*delbxp*vp0

     x3(ns,m) = x2(ns,m) + dt*xdot
     y3(ns,m) = y2(ns,m) + dt*ydot
     z3(ns,m) = z2(ns,m) + dt*zdot
     u3(ns,m) = u2(ns,m) + dt*pzdot

     dum = 1-w3(ns,m)*nonlin(ns)*0.
     if(ildu.eq.1)dum = (tgis(ns)/ter)**1.5*exp(vfac*(1/tgis(ns)-1./ter))
     !         vxdum = eyp+vpar/b*delbxp
     vxdum = (eyp/b+vpar/b*delbxp)*dum1
     w3old = w3(ns,m)
     w3(ns,m) = w2(ns,m) + dt*(vxdum*kap+edot/ter)*dum*xnp

     if(abs(w3(ns,m)).gt.1.0.and.nonlin(ns)==1)then
        w3(ns,m) = 0.
        w2(ns,m) = 0.
     end if

     laps=anint((z3(ns,m)/lz)-.5)*(1-peritr)
     r=x3(ns,m)-0.5*lx+lr0
     i = int((r-rin)/dr)
     i = min(i,nr-1)
     i = max(i,0)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     y3(ns,m)=mod(y3(ns,m)-laps*2*pi*qr*lr0/q0*sign(1.0,q0)+8000.*ly,ly)
     if(x3(ns,m)>lx.and.iperidf==0)then
        x3(ns,m) = lx-1.e-8
        z3(ns,m)=lz-z3(ns,m)
        x2(ns,m) = x3(ns,m)
        z2(ns,m) = z3(ns,m)
        w2(ns,m) = 0.
        w3(ns,m) = 0.
     end if
     if(x3(ns,m)<0..and.iperidf==0)then
        x3(ns,m) = 1.e-8
        z3(ns,m)=lz-z3(ns,m)
        x2(ns,m) = x3(ns,m)
        z2(ns,m) = z3(ns,m)
        w2(ns,m) = 0.
        w3(ns,m) = 0.
     end if
     z3(ns,m)=mod(z3(ns,m)+8.*lz,lz)
     x3(ns,m)=mod(x3(ns,m)+800.*lx,lx)
     x3(ns,m) = min(x3(ns,m),lx-1.0e-8)
     y3(ns,m) = min(y3(ns,m),ly-1.0e-8)
     z3(ns,m) = min(z3(ns,m),lz-1.0e-8)

     !     particle diagnostics done here because info is available...
     k = int(x3(ns,m)/(lx/nsubd))
     k = min(k,nsubd-1)
     k = k+1
     mypfl_es(k)=mypfl_es(k) + w3old*(eyp)
     mypfl_em(k)=mypfl_em(k) + w3old*(vpar*delbxp/b)
     myefl_es(k)=myefl_es(k) + vfac*w3old*(eyp)
     myefl_em(k)=myefl_em(k) + vfac*w3old*(vpar*delbxp/b)
     myke=myke + vfac*w3(ns,m)
     mynos=mynos + w3(ns,m)
     myavewi = myavewi+abs(w3(ns,m))

     !     xn+1 becomes xn...
     u2(ns,m)=u3(ns,m)
     x2(ns,m)=x3(ns,m)
     y2(ns,m)=y3(ns,m)
     z2(ns,m)=z3(ns,m)
     w2(ns,m)=w3(ns,m)

     !     100     continue
  enddo
!$acc end parallel
!$acc wait
cpush1_end_tm=cpush1_end_tm+MPI_WTIME()
!$acc update host(z3,x2,x3,y2,y3,z2,u2,u3,w2,w3)

  sbuf(1)=myke
  sbuf(2)=myefl_es(nsubd/2)
  sbuf(3)=mypfl_es(nsubd/2)
  sbuf(4)=mynos
  sbuf(5)=myavewi
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,           &
       MPI_COMM_WORLD,ierr)

  ketemp=rbuf(1)
  efltemp=rbuf(2)
  pfltemp=rbuf(3)
  nostemp=rbuf(4)
  !if(myid ==0 )write(*,*) 'before cpush tmm', real(tmm(1)),real(tmm(1),8)
  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  avewi(ns,n) = (rbuf(5)/ real(tmm(1)))/real(ntube) 
  nos(1,n)=(nostemp/ real(tmm(1)))/real(ntube) 
  ke(1,n)=(ketemp/( 2.*real(tmm(1))*mims(ns) ))/real(ntube)
  !if(myid ==0 )write(*,*) 'after cpush tmm', real(tmm(1)),real(tmm(1),8)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  sbuf(1:nsubd) = myefl_es(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd
     efl_es(ns,k,n)=(rbuf(k)/ real(tmm(1)))/real(ntube) *totvol/vol(k)*cn0s(ns)
  end do

  sbuf(1:nsubd) = myefl_em(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd
     efl_em(ns,k,n)=(rbuf(k)/ real(tmm(1)))/real(ntube) *totvol/vol(k)*cn0s(ns)
  end do

  sbuf(1:nsubd) = mypfl_es(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd
     pfl_es(ns,k,n)=(rbuf(k)/ real(tmm(1)))/real(ntube) *totvol/vol(k)*cn0s(ns)
  end do

  sbuf(1:nsubd) = mypfl_em(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd
     pfl_em(ns,k,n)=(rbuf(k)/ real(tmm(1)))/real(ntube) *totvol/vol(k)*cn0s(ns)
  end do

  !      pfl(1,n)=pfltemp/( real(tmm(1)) )
  !      efl(1,n)=mims(ns)/tets(1)*efltemp/( real(tmm(1)) )

  np_old=mm(ns)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call init_pmove(z3(ns,:),np_old,lz,ierr)
  
  !$acc update device(z3) 
  call pmove(x2(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(x2)
  call pmove(x3(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(x3)
  call pmove(y2(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(y2)
  call pmove(y3(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(y3)
  call pmove(z2(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(z2)
  call pmove(z3(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(z3)
  call pmove(u2(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(u2)
  call pmove(u3(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(u3)
  call pmove(w2(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(w2)
  call pmove(w3(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(w3)
  call pmove(mu(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(mu)
  call pmove(xii(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(xii)
  call pmove(z0i(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(z0i)
  call pmove(pzi(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(pzi)
  call pmove(eki(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(eki)
  call pmove(u0i(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(u0i)
  !$acc wait
  call end_pmove(ierr)
  mm(ns)=np_new
  !     write(*,*)MyId,mm(ns)

  !      return
end subroutine cpush

