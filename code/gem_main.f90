!   global variables...
program gem_main

  use gem_com
  use gem_equil
  use gem_fft_wrapper

  implicit none
  integer :: n,i,j,k,ip

  call initialize
  ! use the following two lines for r-theta contour plot
  if(iCrs_Sec==1)then
     call pol2d
     goto 100
  end if
  if(iget.eq.0)call loader_wrapper
  if(iget.eq.1)then
     call restart(1,0)
  end if

  if(isft==1)then
          !!!without wtime (htc)
     call ftcamp
     goto 100
  end if

  starttm=MPI_WTIME()
  do  timestep=ncurr,nm
     tcurr = tcurr+dt
     accumulate_start_tm=accumulate_start_tm+MPI_WTIME()
     call accumulate(timestep-1,0)
     accumulate_end_tm=accumulate_end_tm+MPI_WTIME()
     poisson_start_tm=poisson_start_tm+MPI_WTIME()
     call poisson(timestep-1,0)
     poisson_end_tm=poisson_end_tm+MPI_WTIME()
     field_start_tm=field_start_tm+MPI_WTIME()
     call field(timestep-1,0)
     field_end_tm=field_end_tm+MPI_WTIME()
     diagnose_start_tm=diagnose_start_tm+MPI_WTIME()
     call diagnose(timestep-1)
     diagnose_end_tm=diagnose_end_tm+MPI_WTIME()
     reporter_start_tm=reporter_start_tm+MPI_WTIME()
     call reporter(timestep-1)
     reporter_end_tm=reporter_end_tm+MPI_WTIME()
     push_start_tm=ppush_start_tm+MPI_WTIME()
     call push_wrapper(timestep,1)
     push_end_tm=push_end_tm+MPI_WTIME()
     accumulate_start_tm=accumulate_start_tm+MPI_WTIME()
     call accumulate(timestep,1)
     accumulate_end_tm=accumulate_end_tm+MPI_WTIME()
     poisson_start_tm=poisson_start_tm+MPI_WTIME()
     call poisson(timestep,1)
     poisson_end_tm=poisson_end_tm+MPI_WTIME()
     field_start_tm=field_start_tm+MPI_WTIME()
     call field(timestep,1)
     field_end_tm=field_end_tm+MPI_WTIME()
     push_start_tm=push_start_tm+MPI_WTIME()
     call push_wrapper(timestep,0)
     push_end_tm=push_end_tm+MPI_WTIME()
     
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  end do
  lasttm=MPI_WTIME()
  !ftcamp_start_tm=ftcamp_start_tm+MPI_WTIME()
  !call ftcamp
  !ftcamp_end_tm=ftcamp_end_tm+MPI_WTIME()
  tottm=lasttm-starttm
  accumulate_tot_tm=accumulate_end_tm-accumulate_start_tm
  poisson_tot_tm=poisson_end_tm-poisson_start_tm
  field_tot_tm=field_end_tm-field_start_tm
  diagnose_tot_tm=diagnose_end_tm-diagnose_start_tm
  reporter_tot_tm=reporter_end_tm-reporter_start_tm
  push_tot_tm=push_end_tm-push_start_tm
  ppush_tot_tm=ppush_end_tm-ppush_start_tm
  cpush_tot_tm=cpush_end_tm-cpush_start_tm
  grid1_tot_tm=grid1_end_tm-grid1_start_tm
  ppush1_tot_tm=ppush1_end_tm-ppush1_start_tm
  cpush1_tot_tm=cpush1_end_tm-cpush1_start_tm
  grid11_tot_tm=grid11_end_tm-grid11_start_tm
  poisson0_tot_tm=poisson0_end_tm-poisson0_start_tm
  exact_total_tm=tottm-poisson0_tot_tm
  if(myid==0)then
  write(*,*),'exact_total_tm',exact_total_tm,'tottm',tottm,'accumulate_tot_tm',accumulate_tot_tm,'poisson_tot_tm',poisson_tot_tm,'field_tot_tm',field_tot_tm,&
             'diagnose_tot_tm',diagnose_tot_tm,'reporter_tot_tm',reporter_tot_tm,'push_tot_tm',push_tot_tm, 'ppush_tot_tm',ppush_tot_tm, &
             'cpush_tot_tm',cpush_tot_tm,'gird1_tot_tm',grid1_tot_tm,'ppush1_tot_tm',ppush1_tot_tm,&
              'cpush1_tot_tm',cpush1_tot_tm,'grid11_tot_tm',grid11_tot_tm,'poisson0_tot_tm',poisson0_tot_tm
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

100 call MPI_FINALIZE(ierr)

end program gem_main

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine hybinit

  use gem_com
  use gem_equil
  implicit none
  GCLR = int(MyId/ntube)
  GLST = numprocs/ntube-1
  TCLR = modulo(MyId,ntube)
  TLST = ntube-1

  !***MPI variables

  mykm = 1
  rngbr = GCLR+1
  if (GCLR.eq.GLST)rngbr = 0

  lngbr = GCLR-1
  if (GCLR.eq.0) lngbr = GLST

  idnxt = TCLR+1
  if (TCLR.eq.TLST)idnxt = 0

  idprv = TCLR-1
  if (TCLR.eq.0) idprv = TLST
  !     
  !      return
end subroutine hybinit

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine init

  use gem_com
  use gem_equil
  implicit none
  character(len=62) dumchar
  INTEGER :: i,j,k,m,n,ns,idum,i1,j1,k1
  INTEGER ::mm1
  INTEGER ::mm2,lr1
  real :: mims1,tets1,q1,kappan,kappat,r,qr,th,cost,dum,zdum
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp
  real :: grp,gxdgyp,jacp,jfnp,gn0ep,gt0ep,gt0ip,grdgtp,gthp
  real,DIMENSION (1:5) :: gn0sp
  real :: wx0,wx1,wz0,wz1,b

  IU=cmplx(0.,1.)
  pi=4.0*atan(1.0)
  pi2 = pi*2.

  open(115,file='gem.in')
  read(115,*) dumchar
  read(115,*) itube,iperi,iperidf,ibunit,mimp,mcmp,chgi,chgc
  read(115,*) dumchar
  read(115,*) Rovera,elon0,selon0,tria0,stria0,rmaj0p,q0,shat0,teti,tcti,rhoia
  read(115,*) dumchar
  read(115,*) Rovlni,Rovlti,Rovlne,Rovlte,Rovlnc,Rovltc,ncne
  read(115,*) dumchar
  read(115,*) imx,jmx,kmx,mmx,mmxe,nmx,nsmx,modemx,ntube,icrs_sec,ipg,isphi
  read(115,*) dumchar
  read(115,*) lxa,lymult,delra,delri,delre,delrn,nrst,eprs,nlow,jcnt
  im=imx;jm=jmx;km=kmx
  read(115,*) dumchar
  read(115,*) dt,nm,nsm,xshape,yshape,zshape
  read(115,*) dumchar
  read(115,*) iput,iget,ision,isiap,peritr,llk,mlk,onemd,izonal,adiabatic_electron,ineq0,iflut ! jycheng
  read(115,*) dumchar
  read(115,*) nplot,xnplt,modem,nzcrt,npze,npzi,npzc,npzb
  read(115,*) dumchar
  read(115,*) isft,mynf,frmax,ifskp,iphbf,iapbf,idpbf
  read(115,*) dumchar
  read(115,*) cut,amp,tor,ishift,fradi,kxcut,kycut,bcut
  read(115,*) dumchar
  read(115,*) r0a,width,vpp,vt0,yd0
  read(115,*) dumchar
  read(115,*) c4,ifluid,isg,amie,nuacs,rneui,vexbsw,vparsw,mach,gamma_E
  read(115,*) dumchar
  read(115,*) betai,nonlin(1),nonlin(2),nonline,ipara,vwidth,vwidthe,vcut,isuni,idg

  !$acc update device(nonlin)
  if(isft==1) iget = 1
  nfreq = kmx*mynf
  call new_gem_com()
  !      rin = r0-lx/2
  !      rout = r0+lx/2
  rina = r0a-lxa/2.
  routa = r0a+lxa/2.

  call new_equil()
  !!$acc update device(jfn,radius,t0s,tgis,dbdr,gr,bdcrvb,qhat,f,grdgt,capts,grcgt,dydr,gxdgy,phincp,sf,capns,bfld,xn0s,vparsp,thfnz,psip,psip2,curvbz,dipdr,dbdth)
  lx=lxa*a
  ly=2.*pi*r0/q0abs/lymult
  !beta has been used all the time, so the definition is fixed, declared in gem_com.f. betai is declared in equil.f
  !in gem.in betai is intended to be GYRO's beta_e (within a factor of 2). It is defined to be (mu0 ne Te)/Bunit**2
  beta=betai 
  br0 = rmaj0
  lr0 = r0
  qp = q0p
  lz = pi2*q0abs*rmaj0
  delz = lz/ntheta
  rneu = nuacs
  if(onemd==1)then
     kxcut=pi2*shat0*pi2/ly*2.1
     kycut=pi2/ly*1.1
  end if
  achii=0.
  achie=0.
  addi=0.

  if(myid.eq.master)then
     open(19, file='eqdat', status='unknown')
     write(19,*) 'r0a,lxa,beta,a,rmaj0=', r0a,lxa,beta,a,rmaj0
     write(19,*)'elon0,tria0,rmaj0p,stria0,selon0,q0,shat=', elon0,tria0,rmaj0p,stria0,selon0,q0,r0*q0p/q0
     write(19,*)'xn0s(1,nr/2),xn0s(2,nr/2),t0s(1,nr/2),t0s(2,nr/2)=',xn0s(1,nr/2),xn0s(2,nr/2),t0s(1,nr/2),t0s(2,nr/2)
     write(19,*)'capns(1,nr/2),capns(2,nr/2),capts(1,nr/2),capts(2,nr/2)=', capns(1,nr/2),capns(2,nr/2),capts(1,nr/2),capts(2,nr/2)
     write(19,*)'capne(nr/2),capte(nr/2)=', capne(nr/2),capte(nr/2)
     write(19,*)'t0i0,t0e0,ni0,ne0=', t0i(nr/2),t0e(nr/2),xn0i(nr/2),xn0e(nr/2)
     write(19,*)'a/cs=', a/sqrt(t0e(nr/2)/mimp)
     write(19,*)'i,   sf,   ni,   ne,   nc,   ti,   tc,   capni,   capnc, captc'
     do i = 0,nr
        write(19,99)i,sf(i),xn0i(i),xn0e(i),t0i(i),t0e(i),capni(i),capne(i),capti(i),capte(i)
     end do
     write(19,99)i,cn0e,cn0s(1),cn0s(2),cn0s(3)
     close(19)
  end if
98 format(7(1x,e14.7))
99 format(1x,i3,2x,15(1x,e10.3))
  ! equilibrium data for transport coefficients calculation by mablab with flux data
  if(myid.eq.master)then
     open(912, file='eqflux', status='unknown')
     do j = 1, nsubd
        i = int(nr*(2*j-1)/(2*nsubd))
        write(912,9991) xn0e(i)*cn0e,capne(i),xn0c(i)*cn0c,capnc(i),t0i(i),capti(i)
     enddo
     close(912)
  endif
9991 format(6(2x,e12.4)) 
  !
  !      write(*,*)'br0,lr0,q0,qp = ', br0,lr0,q0,qp

  if(isuni.eq.0)vwidthe=vwidth
  dte = dt
  iadi = 0
  if(isg.gt.0.)fradi = isg
  if(ifluid.eq.0)then
     iadi = 1
     fradi = cn0e*xn0e(nr/2)/t0e(nr/2)
  end if

  !     begin reading species info, ns=1,nsm...
  if(nsm.le.0) write(*,*)'invalid nsm',nsm
  read(115,*) dumchar
  ns = 1
  read(115,*) dumchar
  read(115,*) mm1,mm2,mims(3),q(3),lr1
  mims(1)=mimp;mims(2)=mcmp;q(1)=chgi;q(2)=chgc
  !$acc update device(mims,q)
  read(115,*) dumchar
  read(115,*)iflr,iorb
  !!!mm1 in the new version of gem is the particle number per ntube
  tmm(1)=mm1
  tmm(2)=mm2
  mm(:)=int(mm1/numprocs)*ntube
  mme = int(mm2/numprocs)
  !!!mm1 in the new version of gem is the particle number per ntube
  !if (MyId.eq.Last) mm(ns)=int(mm1-int(Last,8)*int(mm(ns),8),4)
  if (MyId.eq.Last) mm(ns)=int(mm1-int(mm1/numprocs)*Last)*ntube
  
       !write(*,*)'in init  ',Myid,mm1-int8(Last)*int8(mm(ns)),mm(ns)
  tets(1)=1
  lr(1)=lr1
  lr(2)=lr1
  lr(3)=lr1
  !$acc update device(lr)
  pzcrite = abs(psi(nr)-psi(0))/br0/npze
  encrit = 1.
  pzcrit(1) = q(1)*abs(psi(nr)-psi(0))/br0/npzi
  pzcrit(2) = q(2)*abs(psi(nr)-psi(0))/br0/npzc
  pzcrit(3) = q(3)*abs(psi(nr)-psi(0))/br0/npzb

  kapn(ns)=kappan
  kapt(ns)=kappat

  emass = 1./amie
  qel = -1.

  mbeam = 2
  qbeam = 1

  if(iget.eq.1) amp=0.
  !     totvol is the square for now...
  dx=lx/real(im)
  dy=ly/real(jm)
  dz=lz/real(km)
  !      totvol=lx*ly*lz

  e0=lr0/q0/br0
  !     
  do i=0,nxpp
     xg(i)=i*dx !dx*(tclr*nxpp+i)
  enddo
  do j=0,jm
     yg(j)=dy*real(j)
  enddo
  kcnt=1

  do k=0,mykm
     n=GCLR*kcnt+k 
     zg(k)=dz*real(n)
  enddo

  !      jcnt = 3  !jmx/ntube                                                                                                                                                                         
  mstart = 0
  ntor0 = mstart+1
  do m = 0,jcnt-1
     isgnft(m) = 1
     j1 = mstart+int((real(m)+1.0)/2)
     jft(m) = j1
     if(m==0)then
        isgnft(m) = 1
        jft(m) = 0
     end if
     if(m>0.and.modulo(m,2)==0)then
        isgnft(m) = -1
        jft(m) = jmx-j1
     end if
  end do

  !     set arrays for mode history plots...
  lmode(1)=0
  lmode(2)=0
  lmode(3)=0
  lmode(4)=0
  mmode(1)=1 !int(.33*ly/2/pi)-1
  mmode(2)=2 !int(.33*ly/2/pi)
  mmode(3)=3 !int(.33*ly/2/pi)+1
  mmode(4)=4 !int(.33*ly/2/pi)+2
  nmode(1)=km/2
  nmode(2)=km/2
  nmode(3)=km/2
  nmode(4)=km/2

  !     initialize bfld   

  zfnth(0) = 0.
  do j = 1,ntheta
     zfnth(j) = zfnth(j-1)+dth*q0*br0*(1./jfn(j-1)+1./jfn(j))/2
  end do
  if(q0<0)zfnth = zfnth+lz

  thfnz(0) = -pi
  thfnz(ntheta) = pi
  if(q0<0.)then
     thfnz(0) = pi
     thfnz(ntheta) = -pi
  end if

  if(q0>0)then
     k = 0
     do j = 1,ntheta-1
        zdum = j*lz/ntheta
        do i = k,ntheta-1
           if(zfnth(i)<=zdum.and.zfnth(i+1)>zdum)then
              k = i
              dum = (zdum-zfnth(i))*dth/(zfnth(i+1)-zfnth(i))
              thfnz(j) = i*dth-pi+dum
              go to 127
           end if
        end do
127     continue
     end do
  end if

  if(q0<0)then
     k = 0
     do j = 1,ntheta-1
        zdum = lz-j*lz/ntheta
        do i = k,ntheta-1
           if(zfnth(i)>=zdum.and.zfnth(i+1)<zdum)then
              k = i
              dum = (zdum-zfnth(i))*dth/(zfnth(i+1)-zfnth(i))
              thfnz(ntheta-j) = i*dth-pi+dum
              go to 128
           end if
        end do
128     continue
     end do
  end if

  do i1 = 0,nxpp
     r = xg(i1)-0.5*lx+lr0
     do k1 = 0,mykm
        k = int(zg(k1)/delz)
        k = min(k,ntheta-1)
        wz0 = ((k+1)*delz-zg(k1))/delz
        wz1 = 1-wz0
        th = wz0*thfnz(k)+wz1*thfnz(k+1)
        i = int((r-rin)/dr)
        i = min(i,nr-1)
        wx0 = (rin+(i+1)*dr-r)/dr
        wx1 = 1.-wx0
        k = int((th+pi)/dth)
        k = min(k,ntheta-1)
        wz0 = (-pi+(k+1)*dth-th)/dth
        wz1 = 1.-wz0
        dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
             +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
        dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
             +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
        grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
             +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 

        grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
             +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 
        gthp = wx0*wz0*gth(i,k)+wx0*wz1*gth(i,k+1) &
             +wx1*wz0*gth(i+1,k)+wx1*wz1*gth(i+1,k+1) 

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
        jacp = wx0*wz0*jacob(i,k)+wx0*wz1*jacob(i,k+1) &
             +wx1*wz0*jacob(i+1,k)+wx1*wz1*jacob(i+1,k+1) 
        gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
             +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
        fp = wx0*f(i)+wx1*f(i+1)        
        jfnp = wz0*jfn(k)+wz1*jfn(k+1)
        psipp = wx0*psip(i)+wx1*psip(i+1)        
        gn0ep = wx0*xn0e(i)+wx1*xn0e(i+1)        
        gt0ep = wx0*t0e(i)+wx1*t0e(i+1)        
        do ns = 1, nsm
           gn0sp(ns) = wx0*xn0s(ns,i)+wx1*xn0s(ns,i+1)
           gn0s(ns,i1) = gn0sp(ns)
        enddo
        gt0ip = wx0*t0s(1,i)+wx1*t0s(1,i+1)        
        b=1.-tor+tor*bfldp
        cfx(i1,k1) = br0/b**3*fp/radiusp*dbdtp*grcgtp
        cfy(i1,k1) = br0/b**3*fp/radiusp* &
             (dydrp*dbdtp-lr0/q0*qhatp*dbdrp)*grcgtp
        bdgxcgy(i1,k1) = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp             
        bmag(i1,k1) = b
        jac(i1,k1) = jacp*jfnp
        bdgrzn(i1,k1) = q0*br0/radiusp/b*psipp*grcgtp/jfnp
        gn0e(i1) = gn0ep
        gt0e(i1) = gt0ep
        gt0i(i1) = gt0ip

        ggxdgy(i1,k1) = gxdgyp
        ggy2(i1,k1) = dydrp**2*grp**2 + (r0/q0*qhatp*gthp)**2 + 2*dydrp*r0/q0*qhatp*grdgtp
        ggx(i1,k1) = grp
     end do
  end do

  iseed = -(1777+myid*13)
  idum = ran2(iseed)
  phi = 0.
  apar = 0.
  dene = 0.
  upar = 0.
  upa0 = 0.
  camp = 0.
  delbx=0.
  dpdz(:,:,:)=0.
  delby(:,:,:)=0.
  dadz(:,:,:)=0.0
  !$acc update device(phi)
  !if ES only
  !$acc update device(apar,dpdz,delby,delbx,dadz)

  if(myid.eq.master)then
     write(*,*)zfnth(ntheta),thfnz(ntheta/2),thfnz(ntheta/2+1)
     if(myid.eq.master)open(9, file='plot', &
          status='unknown',position='append')
     write(9,*)'dt,beta= ',dt, beta
     write(9,*)'amp,vpp,yd0 = ',amp,vpp,yd0
     write(9,*)'peritr,ifluid= ',peritr,ifluid
     write(9,*)'tor,nonlin= ',tor,nonlin(1),nonlin(2)
     write(9,*)'isuni= ',isuni, 'amie= ',amie
     write(9,*)'kxcut,kycut,bcut= ',kxcut,kycut,bcut
     write(9,*)'fradi,isg= ',fradi,isg
     write(9,*)'llk,mlk,onemd =',llk,mlk,onemd
     write(9,*)'vwidth (i,e), vcut= ', vwidth,vwidthe,vcut
     write(9,*)'rneu= ', rneu
     write(9,*)'V-ExB switch= ', vexbsw
     write(9,*)'V-parallel switch= ', vparsw
     write(9,*)'mm1= ',mm1
     write(9,*)'pzcrite,encrit = ',pzcrite,encrit
     write(9,*) 'lxa,lymult,delra,r0a,rina,routa=',lxa,lymult,delra,r0a,rina,routa
     write(9,*) 'a,r0,rmaj0,q0,lx,ly,lz=',a,r0,rmaj0,q0,lx,ly,lz
     write(9,*) 't0,kyrhoi_local=',t0i(nr/2),2*pi*sqrt(mims(1))*sqrt(t0i(nr/2))/ly
     close(9)
  end if

  if(myid.eq.master)then
     write(*,*)'dt,beta= ',dt, beta
     write(*,*)'amp,vpp,yd0 = ',amp,vpp,yd0
     write(*,*)'peritr,ifluid= ',peritr,ifluid
     write(*,*)'tor,nonlin= ',tor,nonlin(1),nonlin(2)
     write(*,*)'isuni= ',isuni, 'amie= ',amie
     write(*,*)'kxcut,kycut,bcut= ',kxcut,kycut,bcut
     write(*,*)'fradi,isg= ',fradi,isg
     write(*,*)'llk,mlk,onemd =',llk,mlk,onemd
     write(*,*)'vwidth (i,e), vcut= ', vwidth,vwidthe,vcut
     write(*,*)'rneu= ', rneu
     write(*,*)'V-ExB switch= ', vexbsw
     write(*,*)'V-parallel switch= ', vparsw
     write(*,*)'mm1= ',mm1
     write(*,*)'pzcrite,encrit = ',pzcrite,encrit
     write(*,*)'nue0 = ',nue0(1),nue0(nr/2),nue0(nr-1)
     write(*,*)'xn0e(1),xnir0 = ',xn0e(1),xnir0
     write(*,*)'frequ, eru = ', frequ, eru
     write(*,*) 'lxa,lymult,delra,r0a,rina,routa=',lxa,lymult,delra,r0a,rina,routa
     write(*,*) 'a,r0,rmaj0,q0,lx,ly,lz=',a,r0,rmaj0,q0,lx,ly,lz
     write(*,*) 't0,kyrhoi_local=',t0i(nr/2),2*pi*sqrt(mims(1))*sqrt(t0i(nr/2))/ly
     write(*,*) 'coefu = ', xu**2*frequ
     write(*,*) 'ktheta*rhos = ',2*pi*sqrt(mims(1))*sqrt(t0e(nr/2))/ly
     write(*,*) 'cs/a, q0, q0p, s^hat = ',sqrt(t0e(nr/2)/2.)/a, q0, q0p, q0p/q0*r0
     write(*,*) 'rho* = rhos/a = ', sqrt(mims(1))*sqrt(t0e(nr/2))/a
     write(*,*) 'f0p,psip(nr/2),Bunit,candyf0p = ',f0p,psip(nr/2),bunit,candyf0p
     write(*,*) 'lxa min = ', ly*q0/(2*pi*r0*q0p)/a
     write(*,*) 't0i(nr/2)= ', t0i(nr/2)
     write(*,*) 'Gyrokrs = ', 2*pi*sqrt(mims(1))*sqrt(t0e(nr/2))/ly/bunit
  end if
  close(115)
!$acc update device(jfn,radius,t0s,tgis,dbdr,gr,bdcrvb,qhat,f,grdgt,capts,grcgt,dydr,gxdgy,phincp,sf,capns,bfld,xn0s,vparsp,thfnz,psip,psip2,curvbz,dipdr,dbdth)
  !      return
end subroutine init

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine grad(ip)

  !  currently set up for periodic in x,y,z

  use gem_com
  use gem_equil
  implicit none
  INTEGER :: i,j,k,ip
  real :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)
  real :: tmp(0:imx,0:jmx,0:1)

  !!$acc update host(ex,ey,ez)
  if(idg==1)write(*,*)'enter grad'
  call gradu(phi(:,:,:),-ex,-ey)
  !!$acc update device(ex,ey,ez)
  !!$acc update host(dadz,dpdz,apar,delbx,delby)
  delbx = 0.
  delby = 0.
  if(ifluid.eq.1)then
     call gradu(apar(:,:,:),ux,uy)
     delbx(:,:,:) = uy(:,:,:)
     delby(:,:,:) = -ux(:,:,:)
  end if
  !!$acc update device(dadz,dpdz,apar,delbx,delby)

  !      return
end subroutine grad

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine modes2(u,modehis,n)

  !     calculate mode histories. calculate modes for u at timestep n
  !     and store in modehis(mode,n).

  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  !     
  real :: u(0:imx,0:jmx,0:1)
  COMPLEX :: modehis(modemx,0:nmx)
  COMPLEX :: modebuf
  INTEGER :: n,i,j,k,l,m,ii

  INTEGER :: mode,jj,thek,oproc,ifirst

  !     

  if(n.eq.0) return

  do mode=1,modem
     oproc=int(nmode(mode)/kcnt*ntube)

     if (MyId.eq.oproc) then
        thek=0
        do j=0,jm-1
           do i=0,im-1
              tmpx(i)=u(i,j,thek)
           enddo

           !     FT in x....
           call ccfft('x',1,imx,1.0,tmpx,coefx,workx,0)
           ii=lmode(mode) !+1
           if(lmode(mode).lt.0) write(*,*) 'lmode < 0, error'
           tmpy(j)=tmpx(ii)/real(im)
        enddo

        !     FT in y....
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        jj=mmode(mode)  !+1
        if(mmode(mode).lt.0) write(*,*) 'mmode < 0, error'
        modebuf=tmpy(jj)/real(jm)

     endif

     call MPI_BCAST(modebuf,1,MPI_DOUBLE_COMPLEX,oproc, &
          MPI_COMM_WORLD,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     !         write(*,*)myid,modebuf
     modehis(mode,n)=modebuf
  enddo
  !     
  !      return
end subroutine modes2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine restart(iflag,n)

  use gem_com
  use gem_equil
  implicit none
  INTEGER :: m,ns,iflag,n,i,j,k,ip
  character(len=70) fname
  character(len=5) holdmyid

  write(holdmyid,'(I5.5)') MyId
  fname=directory//'dump_'//holdmyid//'.b'
  if(iflag.eq.1) then
     open(139+MyId,file=fname,form='unformatted',status='old')
     read(139+MyId)ncurr
     read(139+MyId)tcurr,rmpp,rmaa
     do ns = 1,nsm
        read(139+MyId)mm(ns)
        do m=1,mm(ns)
           read(139+MyId) mu(m,ns)
           read(139+MyId) x2(m,ns),y2(m,ns),z2(m,ns),u2(m,ns),w2(m,ns)
           read(139+MyId) xii(m,ns),z0i(m,ns),pzi(m,ns),eki(m,ns),u0i(m,ns)
           w2(m,ns)=w2(m,ns)/cut
           x3(m,ns)=x2(m,ns)
           y3(m,ns)=y2(m,ns)
           z3(m,ns)=z2(m,ns)
           u3(m,ns)=u2(m,ns)
           w3(m,ns)=w2(m,ns)
        enddo
     end do
120  continue

     read(139+MyId)mme
     do  m=1,mme
        read(139+MyId) x2e(m),y2e(m),z2e(m),u2e(m),w2e(m),mue2(m),ipass(m)
        read(139+MyId) xie(m),z0e(m),pze(m),eke(m),mue(m),u0e(m)
        w2e(m)=w2e(m)/cut
        x3e(m)=x2e(m)
        y3e(m)=y2e(m)
        z3e(m)=z2e(m)
        u3e(m)=u2e(m)
        w3e(m)=w2e(m)
        mue3(m)=mue2(m)
     end do

     do i = 0,6
        do j = 0,50000
           read(139+myid)camp(i,j)
        end do
     end do
     close(139+MyId)
  endif

  if(iflag.eq.2) then
     open(139+MyId,file=fname,form='unformatted',status='unknown')
     write(139+MyId)n+1
     write(139+MyId)tcurr-dt,rmpp,rmaa
     do ns = 1,nsm
        write(139+MyId)mm(ns)
        do m=1,mm(ns)
           write(139+MyId) mu(m,ns)
           write(139+MyId) x2(m,ns),y2(m,ns),z2(m,ns),u2(m,ns),w2(m,ns)
           write(139+MyId) xii(m,ns),z0i(m,ns),pzi(m,ns),eki(m,ns),u0i(m,ns)
        enddo
     end do

     write(139+MyId)mme
     do  m=1,mme
        write(139+MyId) x2e(m),y2e(m),z2e(m),u2e(m),w2e(m),mue2(m),ipass(m)
        write(139+MyId) xie(m),z0e(m),pze(m),eke(m),mue(m),u0e(m)
     end do

     do i = 0,6
        do j = 0,50000
           write(139+myid)camp(i,j)
        end do
     end do
     close(139+MyId)
  endif

end subroutine restart
!-----------------------------------------------------------------------

!        Normal distribution random no. generator, stand. dev. = 1.
!        Version 2 does it Willy's way...

subroutine parperp(vpar,vperp2,m,pi,cnt,MyId)
  use gem_com, only : iseed2
  INTERFACE
     real function ran(i)
       integer :: i
     end function ran
  END INTERFACE
        

  !INTERFACE
     !real function revers(num,n)
       !integer :: num,n
     !end function revers
  !END INTERFACE

  real :: vpar,vperp2,r1,r2,t,pi
  INTEGER :: m,iflag,cnt,MyId
  real :: c0,c1,c2
  real :: d1,d2,d3
  data c0,c1,c2/2.515517,0.802853,0.010328/
  data d1,d2,d3/1.432788,0.189269,0.001308/


  !r1=revers(m+MyId*cnt,7)
  !r2=revers(m+MyId*cnt,11)

  !!!random generator for large particle number myid*cnt>2*10^31(htc)
  r1=ran(iseed2)
  r2=ran(iseed2)
  


  !.....quiet start---see denavit pf '71(?) & abramowitz hand book
  !.....fibonacci start---see denavit comm. pla. phy. & con. fus. '81
  ! warning: we have g1=1 in the x-direction. This surpresses all odd
  !          modes in the x-direction!!!

  iflag=1
  if(r1.le.0.5) go to 110
  r1=1.-r1
  iflag=-1
110 continue
  if(r1.ge.1.e-6) then
     t=sqrt(log(1.0/(r1*r1)))
  else
     t=5.0
     !write(*,*)'parperp2 warning  m= ',m
  endif

  vpar=t-(c0+c1*t+c2*t**2)/(1.+d1*t+d2*t**2+d3*t**3)
  vpar=vpar*iflag

  vperp2=-2.0*log(r2)

  !        return
end subroutine parperp

!---------------------------------------------------------------------- 
!    calculate weights and delta j for proper periodicity 
!    in the toroidal direction, weights and delta j are
!    dependant only on minor radius, hence x, hence
!    weight is a vector in 0:imx

subroutine weight

  use gem_com
  use gem_equil
  implicit none
  INTEGER :: i,j
  real :: dely,aweight,r,qr,wx0,wx1

  !         peritr = 0
  if (GCLR.eq.Master.and.peritr.eq.0) then
     do i=0,im
        r=real(i)*dx-0.5*lx+lr0
        j = int((r-rin)/dr)
        j = min(j,nr-1)
        wx0 = (rin+(j+1)*dr-r)/dr
        wx1 = 1.-wx0
        qr = wx0*sf(j)+wx1*sf(j+1)
        dely=modulo(2.*pi*lr0/q0*qr*sign(1.0,q0)+800.0*ly,ly)
        deljp(i)=int(dely/dy)
        deljm(i)=0
        aweight=modulo(dely,dy)
        weightp(i)=1.-aweight/dy
        weightm(i)=1.
     enddo
  elseif (GCLR.eq.GLST.and.peritr.eq.0) then
     do i=0,im
        r=real(i)*dx-0.5*lx+lr0
        j = int((r-rin)/dr)
        j = min(j,nr-1)
        wx0 = (rin+(j+1)*dr-r)/dr
        wx1 = 1.-wx0
        qr = wx0*sf(j)+wx1*sf(j+1)
        dely=modulo(2.*pi*lr0/q0*qr*sign(1.0,q0)+800.*ly,ly)
        deljm(i)=int(dely/dy)
        deljp(i)=0
        aweight=modulo(dely,dy)
        weightm(i)=1.-aweight/dy
        weightp(i)=1.
     enddo
  else
     do i=0,im
        deljp(i)=0
        deljm(i)=0
        weightp(i)=1.
        weightm(i)=1.
     enddo
  endif

  do j=0,jm
     do i=0,im
        jpl(i,j)=modulo(j+deljp(i)+8*jm,jm)
        jpn(i,j)=modulo(j+deljp(i)+1+8*jm,jm)
        jmi(i,j)=modulo(j-deljm(i)+8*jm,jm)
        jmn(i,j)=modulo(j-deljm(i)-1+8*jm,jm)
        weightpn(i)=1.-weightp(i)
        weightmn(i)=1.-weightm(i)
     enddo
  enddo

  !      return
end subroutine weight

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine eqmo(ip)

  use gem_com
  use gem_equil

  implicit none
  real :: lbfr(0:imx,0:jmx)
  real :: lbfs(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  integer :: i,j,k,ip
  real :: eta

  ez(:,:,:) = 0.
  dpdz = 0.
  do i = 0,im
     do j = 0,jm
        do k = 1,mykm-1
           ez(i,j,k) = (phi(i,j,k-1)-phi(i,j,k+1))/(2.*dz)
        end do
     end do
  end do

  rbfs = phi(:,:,0)

  call MPI_SENDRECV(rbfs(0,0),(imx+1)*(jmx+1), &
       MPI_REAL8,rngbr,204,                    &
       lbfr(0,0),(imx+1)*(jmx+1),              &
       MPI_REAL8,lngbr,204,                    &
       TUBE_COMM,stat,ierr)

  lbfs=phi(:,:,1)

  call MPI_SENDRECV(lbfs(0,0),(imx+1)*(jmx+1), &
       MPI_REAL8,lngbr,205,                    &
       rbfr(0,0),(imx+1)*(jmx+1),              &
       MPI_REAL8,rngbr,205,                    &
       TUBE_COMM,stat,ierr)                    
  call MPI_BARRIER(TUBE_COMM,ierr)             
  do i=0,im
     do j=0,jm
        ez(i,j,0)=(weightp(i)*lbfr(i,jpl(i,j)) &
             +weightpn(i)*lbfr(i,jpn(i,j))        &
             -phi(i,j,1))/(2.*dz)

        ez(i,j,1)=( phi(i,j,0)-(weightm(i)*rbfr(i,jmi(i,j)) &
             +weightmn(i)*rbfr(i,jmn(i,j))))/(2.*dz)

     enddo
  enddo

  dpdz(:,:,:) = -ez(:,:,:)

end subroutine eqmo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine spec(n)
  use gem_com
  use gem_equil
  implicit none
  integer :: i,j,k,l,m,n
  real :: pf,efe,efi,pfi,efc,pfc,efb,pfb,pflxgb,eflxgb,x
  real :: pf_em,efe_em,efi_em,pfi_em,efc_em,pfc_em,efb_em,pfb_em
  real :: tdum,v(0:imx)

!  call zon(upart+amie*upa0t,v)
!  call joule

  eflxgb = xn0e(nr2)*cn0e*t0e(nr2)*sqrt(t0e(nr2)/mimp)*rhoia**2
  pflxgb = eflxgb/t0e(nr2)
  i = tcurr-dt
  pf = 0.
  efe = 0.
  efi = 0.
  pfi = 0.
  efc = 0.
  pfc = 0.
  efb = 0.
  pfb = 0.

  pf_em = 0.
  efe_em = 0.
  efi_em = 0.
  pfi_em = 0.
  efc_em = 0.
  pfc_em = 0.
  efb_em = 0.
  pfb_em = 0.
  k = 2
  x = real(nsubd)/real(nsubd-2*k)
  do j = 1+k,nsubd-k
     pf = pf+pfle_es(j,n)*vol(j)/totvol*x
     efe = efe+efle_es(j,n)*vol(j)/totvol*x
     efi = efi+efl_es(1,j,n)*vol(j)/totvol*x
     pfi = pfi+pfl_es(1,j,n)*vol(j)/totvol*x
     efc = efc+efl_es(2,j,n)*vol(j)/totvol*x
     pfc = pfc+pfl_es(2,j,n)*vol(j)/totvol*x
     efb = efb+efl_es(3,j,n)*vol(j)/totvol*x
     pfb = pfb+pfl_es(3,j,n)*vol(j)/totvol*x
  end do

  do j = 1+k,nsubd-k
     pf_em = pf_em+pfle_em(j,n)*vol(j)/totvol*x
     efe_em = efe_em+efle_em(j,n)*vol(j)/totvol*x
     efi_em = efi_em+efl_em(1,j,n)*vol(j)/totvol*x
     pfi_em = pfi_em+pfl_em(1,j,n)*vol(j)/totvol*x
     efc_em = efc_em+efl_em(2,j,n)*vol(j)/totvol*x
     pfc_em = pfc_em+pfl_em(2,j,n)*vol(j)/totvol*x
     efb_em = efb_em+efl_em(3,j,n)*vol(j)/totvol*x
     pfb_em = pfb_em+pfl_em(3,j,n)*vol(j)/totvol*x
  end do

  tdum = tcurr-dt
  if(myid.eq.master)then
     open(9, file='plot', status='unknown',position='append')
     open(11, file='flux', status='unknown',position='append')
     open(17, file='yyre', status='unknown',position='append')

     write(*,10)i,rmsphi(n),rmsapa(n),pf,efe,pfi,efi,avewi(1,n),&
          avewe(n),yyre(1,0),yyim(1,0),yyamp(1,0)

10   format(1x,i6,12(2x,e10.3))
11   format(6x,5(2x,e12.5))
12   format(1x,i6,12(2x,e12.5))
13   format(1x,i6,4(2x,e10.3),2x,i7,2x,i7)
15   format(1x,i6,8(2x,e10.3))

     write(9,10)i,rmsphi(n),rmsapa(n),pf,efe,pfi,efi,avewi(1,n),&
          avewe(n),yyre(1,0),yyim(1,0),yyamp(1,0)

     write(11,12)i,pf/pflxgb,pfi/pflxgb,pfc/pflxgb,efe/eflxgb,efi/eflxgb,&
          efc/eflxgb,pf_em/pflxgb,pfi_em/pflxgb,pfc_em/pflxgb,efe_em/eflxgb,&
          efi_em/eflxgb,efc_em/eflxgb

     write(17,12)i,yyre(1,0),yyre(1,1),yyre(1,2),yyre(1,3),yyre(1,4)
     close(9)
     close(11)
     close(17)
  end if

  return
  if(gclr==kmx/2 .and. tclr==0)then
     open(22, file='yyre2', status='unknown',position='append')
     open(23, file='mdhis', status='unknown',position='append')
     open(24, file='mdhisa', status='unknown',position='append')
     open(25, file='stress', status='unknown',position='append')

     write(23,16)tdum,mdhis(0),mdhis(1),mdhis(2),mdhis(3),mdhis(4),&
          mdhisa(0),mdhisa(1),mdhisa(2),mdhisa(3),mdhisa(4)
     write(24,16)tdum,mdhisb(0),mdhisb(1),mdhisc(0),mdhisc(1),&
          mdhisd(0),mdhisd(1)
     write(25,17)tdum,(v(i),i = 0,imx-1)

     do  i = 0,6
        write(22,14)tdum,i,real(phihis(i,0)),(real(phihis(i,j)), &
             aimag(phihis(i,j)), j = 1,jcnt-2,2)
        write(22,14)tdum,i,real(aparhis(i,0)),(real(aparhis(i,j)), &
             aimag(aparhis(i,j)), j = 1,jcnt-2,2)
     end do
     close(22)
     close(23)
     close(24)
     close(25)
14   format(1x,f10.1,1x,i2,10(2x,e12.5))
16   format(1x,f10.1,1x,10(2x,e12.5))
17   format(1x,f10.1,256(2x,e12.5))
  end if

end subroutine spec

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!!new fast random generator (htc)
real function ran(idum)
        !integer, parameter :: k8b=selected_int_kind(18)
        integer,intent(inout) :: idum
        integer,parameter :: IA=16807, IM=2147483647,IQ=127773,IR=2836
        real, save :: am
        integer, save :: ix=-1, iy=-1,k
        
        if (idum <= 0 .or. iy<0) then
                am=nearest(1.0,-1.0)/IM
                iy=ior(ieor(888889999,abs(idum)),1)
                ix=ieor(777755555,abs(idum))
                idum=abs(idum)+1
        end if
        ix=ieor(ix,ishft(ix,13))
        ix=ieor(ix,ishft(ix,-17))
        ix=ieor(ix,ishft(ix,5))
        k=iy/IQ
        iy=IA*(iy-k*IQ)-IR*k
        if (iy < 0) iy=iy+IM
        ran=am*ior(iand(IM,ieor(ix,iy)),1)

end function ran


real function ran2(idum)

  parameter( IM1=2147483563,  &
       IM2=2147483399, &
       AM=1.0/IM1,&
       IMM1=IM1-1,&
       IA1=40014,&
       IA2=40692,&
       IQ1=53668,&
       IQ2=52774,&
       IR1=12211,&
       IR2=3791,&
       NTAB=32,&
       NDIV=1+IMM1/NTAB,&
       EPS=1.2e-7,&
       RNMX=1.0-EPS &
       )

  integer :: j,k,idum2=123456789,iy=0,iv(0:NTAB-1)
  real :: temp

  save idum2, iy,iv
  !      write(*,*)'idum2,iy  ',idum2,iy
  if(idum.le.0)then
     if(-idum.lt.1)then
        idum=1
     else
        idum = -idum
     end if
     idum2 = idum
     do j = NTAB+7,0,-1
        k = idum/IQ1
        idum = IA1*(idum-k*IQ1)-k*IR1
        if(idum.lt.0)idum = idum+IM1
        if(j.lt.NTAB)iv(j) = idum
     end do
     iy = iv(0)
  end if

  k = idum/IQ1
  idum = IA1*(idum-k*IQ1)-k*IR1
  if(idum.lt.0)idum = idum+IM1
  k = idum2/IQ2
  idum2 = IA2*(idum2-k*IQ2)-k*IR2
  if(idum2.lt.0)idum2 = idum2+IM2
  j = iy/NDIV
  iy = iv(j)-idum2
  iv(j) = idum
  if(iy<1)iy = iy+IMM1
  temp = AM*iy
  if(temp>RNMX)then
     ran2 = RNMX
  else
     ran2 = temp
  end if
  return
end function ran2

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine loadi(ns)

  use gem_com
  use gem_equil
  implicit none
  INTEGER :: i,j,k,m,idum,ns,m1
  INTEGER :: np_old,np_new
  real :: vpar,vperp2,r,qr,th,b,cost,ter
  real :: avgv,myavgv,avgw,myavgw
  real :: dumx,dumy,dumz,jacp
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,psp
  real :: grp,gxdgyp,zoldp
  real :: wx0,wx1,wz0,wz1

  !!!new definition of tmm(1) (htc)
  cnt=int(tmm(1)/numprocs)*ntube

  myavgv=0.
  avgv=0.
  avgw = 0.
  myavgw = 0.

  !      ns = 1
  m = 0
  do j = 1,500000000

     !     load a slab of ions...

     dumx=lx*ran2(iseed) !revers(MyId*cnt+j,2) !ran2(iseed)
     dumy=ly*ran2(iseed) !revers(MyId*cnt+j,3) !ran2(iseed)
     dumz=lz*ran2(iseed) !revers(MyId*cnt+j,5) !ran2(iseed)
     dumz = min(dumz,lz-1.e-8)
     r = lr0+dumx-0.5*lx
     th = (dumz-lz/2)/(q0*br0)
     i = int((r-rin)/dr)
     k = int((pi+th)/dth)
     jacp = jacob(i,k)
     if(ran2(iseed)<(0.5*jacp/jacmax))then
        m = m+1
        if(m>mm(ns))goto 170
        x2(m,ns)=min(dumx,lx-1e-8)
        y2(m,ns)=min(dumy,ly-1e-8)

        k = int((th+pi)/dth)
        wz0 = (-pi+(k+1)*dth-th)/dth
        wz1 = 1-wz0
        z2(m,ns) = wz0*zfnth(k)+wz1*zfnth(k+1)
        z2(m,ns)=min(z2(m,ns),lz-1e-8)

        call parperp(vpar,vperp2,m,pi,cnt,MyId)
        !   normalizations will be done in following loop...

        r=x2(m,ns)-0.5*lx+lr0
        cost=cos(th)
        i = int((r-rin)/dr)
        wx0 = (rin+(i+1)*dr-r)/dr
        wx1 = 1.-wx0
        k = int((th+pi)/dth)
        wz0 = (-pi+(k+1)*dth-th)/dth
        wz1 = 1.-wz0
        bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
             +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
        psp = wx0*psi(i)+wx1*psi(i+1)
        ter = wx0*t0s(ns,i)+wx1*t0s(ns,i+1)
        if(ildu==1)ter = tgis(ns)
        b=1.-tor+tor*bfldp

        u2(m,ns)=vpar/sqrt(mims(ns)/ter)
        mu(m,ns)=0.5*vperp2/b*ter
        eki(m,ns) = mu(m,ns)*b+0.5*mims(ns)*u2(m,ns)**2
        pzi(m,ns) = mims(ns)*u2(m,ns)/b-q(ns)*psp/br0
        z0i(m,ns) = z2(m,ns)
        xii(m,ns) = x2(m,ns)
        u0i(m,ns) = u2(m,ns)
        myavgv=myavgv+u2(m,ns)

        !    LINEAR: perturb w(m,ns) to get linear growth...
        !         w2(m,ns)=2.*amp*(revers(MyId*cnt+m,13)-0.5) !(ran2(iseed) - 0.5 )
        w2(m,ns)=2.*amp*sin(pi2/ly*y2(m,ns))*exp(-(z2(m,ns)-lz/2)**2/(lz/8)**2)*exp(-(x2(m,ns)-0.4*lx)**2/(lx/8)**2)
        if(ns==2)w2(m,ns) = 0.
        myavgw=myavgw+w2(m,ns)
     end if
  enddo
170 continue
  myavgw = myavgw/mm(ns)
  !      write(*,*)'myid ', myid,x2(10),y2(20),u2(20),w2(20)
  !    subtract off avg. u...
  call MPI_ALLREDUCE(myavgv,avgv,1, &
       MPI_REAL8, &
       MPI_SUM,MPI_COMM_WORLD,ierr)
  if(idg.eq.1)write(*,*)'all reduce'
  !change int to real (htc)
  !avgv=avgv/real(tmm(1))
  avgv=(avgv/real(tmm(1)))/real(ntube)
  
  do m=1,mm(ns)
     u2(m,ns)=u2(m,ns)-avgv
     x3(m,ns)=x2(m,ns)
     y3(m,ns)=y2(m,ns)
     z3(m,ns)=z2(m,ns)
     u3(m,ns)=u2(m,ns)
     !         w2(m,ns) = w2(m,ns)-myavgw
     w3(m,ns)=w2(m,ns)
  enddo

  np_old=mm(ns)
  call init_pmove(z2(:,ns),np_old,lz,ierr)
  if(idg.eq.1)write(*,*)'pass init_pmove'
  !
  call pmove(x2(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(x2)
  call pmove(x3(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(x3)
  call pmove(y2(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(y2(:,ns))
  call pmove(y3(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(y3(:,ns))
  call pmove(z2(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(z2(:,ns))
  call pmove(z3(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(z3(:,ns))
  call pmove(u2(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(u2(:,ns))
  call pmove(u3(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(u3(:,ns))
  call pmove(w2(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(w2(:,ns))
  call pmove(w3(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(w3(:,ns))
  call pmove(mu(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(mu(:,ns))

  call pmove(xii(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(xii(:,ns))
  call pmove(z0i(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(z0i(:,ns))
  call pmove(pzi(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(pzi(:,ns))
  call pmove(eki(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(eki(:,ns))
  call pmove(u0i(:,ns),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  !$acc update device(u0i(:,ns))

  !     
  call end_pmove(ierr)
  mm(ns)=np_new

  !      return
end subroutine loadi
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine enforce(u)
  use gem_com
  use gem_equil
  implicit none
  INTEGER :: m,i,j,k,ns,jj
  real :: u(0:imx,0:jmx,0:1)      
  real :: lbfs(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  real :: lbfr(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: dum,dum1,dely,th,wy1,ydum

  do j=0,jm-1
     do k=0,mykm
        u(0,j,k) = u(0,j,k)+u(im,j,k)
     enddo
  enddo

  do i=0,im-1 
     do k=0,mykm
        u(i,0,k) = u(i,0,k)+u(i,jm,k)
        u(i,jm,k) = u(i,0,k)
     enddo
  enddo

  rbfs=u(:,:,mykm)
  call MPI_SENDRECV(rbfs(0,0),(imx+1)*(jmx+1), &
       MPI_REAL8, &
       rngbr,101, &
       lbfr(0,0),(imx+1)*(jmx+1), &
       MPI_REAL8, &
       lngbr,101, &
       tube_comm,stat,ierr) 
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  lbfs=u(:,:,0)
  call MPI_SENDRECV(lbfs(0,0),(imx+1)*(jmx+1), &
       MPI_REAL8, &
       lngbr,102, &
       rbfr(0,0),(imx+1)*(jmx+1), &
       MPI_REAL8, &
       rngbr,102, &
       tube_comm,stat,ierr) 
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  do i=0,im
     do j=0,jm
        u(i,j,0)=u(i,j,0)  &
             +weightp(i)*lbfr(i,jpl(i,j))  &
             +weightpn(i)*lbfr(i,jpn(i,j)) 
     enddo
  enddo
  do i=0,im
     do j=0,jm
        u(i,j,mykm)=u(i,j,mykm)           &
             +weightm(i)*rbfr(i,jmi(i,j)) &
             +weightmn(i)*rbfr(i,jmn(i,j))
     enddo
  enddo

  call enfxy(u)
  !      return
end subroutine enforce
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine enfxy(u)
  use gem_com
  use gem_equil
  implicit none
  real, intent(inout) :: u(0:imx,0:jmx,0:1)
  integer :: i,j,k,l,m,n,jj
  real :: ydum,th,dely,wy1

  !    periodic bc in y...
  do k=0,mykm
     do i=0,im-1
        u(i,jm,k)=u(i,0,k)
     enddo
  enddo

  !   bc for x
  do k=0,mykm
     do j=0,jm
        u(im,j,k)=u(0,j,k)
     enddo
  enddo
  
  !      return
end subroutine enfxy
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine gradu(u,ux,uy)
  use gem_com
  use gem_equil
  implicit none
  real, intent(in) :: u(0:imx,0:jmx,0:1)
  real, intent(out) :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)
  integer :: i,j,k,l,m,n,jj,ju,jl
  real :: ydum,th,dely,wy1,ul


  do j=0,jm-1
     ju = j+1
     jl = j-1
     if(j.eq.0)jl = jm-1
     do i=0,im-1
        do k=0,mykm
           uy(i,j,k)=(u(i,ju,k)-u(i,jl,k))/(2.*dy)
        enddo
     enddo
  enddo
  
  do i=1,im-1
     do j=0,jm-1
        do k=0,mykm
           ux(i,j,k)=(u(i+1,j,k)-u(i-1,j,k))/(2.*dx)
        enddo
     enddo
  enddo

  ! do boundary i=0
  do j=0,jm-1
     do k=0,mykm
        ul=u(im-1,j,k)
        ux(0,j,k)=(u(1,j,k)-ul)/(2.*dx)
     enddo
  enddo
  call enfxy(ux)
  call enfxy(uy)
  !      return
end subroutine gradu

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine gradx(u,ux)
  use gem_com
  use gem_equil
  implicit none
  real :: u(0:imx,0:jmx,0:1)
  real :: ux(0:imx,0:jmx,0:1)
  integer :: i,j,k,l,m,n,jj,ju,jl
  real :: ydum,th,dely,wy1,ul

  do i=1,im-1
     do j=0,jm-1
        do k=0,mykm
           ux(i,j,k)=(u(i+1,j,k)-u(i-1,j,k))/(2.*dx)
        enddo
     enddo
  enddo

  ! do boundary i=0
  do j=0,jm-1
     do k=0,mykm
        ul=u(im-1,j,k)
        ux(0,j,k)=(u(1,j,k)-ul)/(2.*dx)
     enddo
  enddo

  call enfxy(ux)

end subroutine gradx

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine grady(u,uy)
  use gem_com
  use gem_equil
  implicit none
  real :: u(0:imx,0:jmx,0:1)
  real :: uy(0:imx,0:jmx,0:1)
  integer :: i,j,k,l,m,n,jj,ju,jl
  real :: ydum,th,dely,wy1,ul

  do j=0,jm-1
     ju = j+1
     jl = j-1
     if(j.eq.0)jl = jm-1
     do i=0,im-1
        do k=0,mykm
           uy(i,j,k)=(u(i,ju,k)-u(i,jl,k))/(2.*dy)
        enddo
     enddo
  enddo

  call enfxy(uy)

end subroutine grady

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine enfz(u)
  use gem_com
  use gem_equil
  implicit none
  real :: u(0:imx,0:jmx,0:1)
  real :: lbfr(0:imx,0:jmx)
  real :: lbfs(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  integer :: i,j

  do i = 0,im
     do j = 0,jm
        rbfs(i,j)=u(i,j,mykm)
     end do
  end do
  call MPI_SENDRECV(rbfs(0,0),(imx+1)*(jmx+1), &
       MPI_REAL8,rngbr,204, &
       lbfr(0,0),(imx+1)*(jmx+1), &
       MPI_REAL8,lngbr,204,       &
       tube_comm,stat,ierr)
  do i = 0,im
     do j = 0,jm
        lbfs(i,j)=u(i,j,0)
     end do
  end do
  call MPI_SENDRECV(lbfs(0,0),(imx+1)*(jmx+1),  &
       MPI_REAL8,lngbr,205,  &
       rbfr(0,0),(imx+1)*(jmx+1), &
       MPI_REAL8,rngbr,205,  &
       tube_comm,stat,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  do i=0,im
     do j=0,jm
        u(i,j,0)=(weightp(i)*lbfr(i,jpl(i,j))  &
             +weightpn(i)*lbfr(i,jpn(i,j))  &
             +u(i,j,0) )/2.  
        u(i,j,mykm)=(weightm(i)*rbfr(i,jmi(i,j)) &  
             +weightmn(i)*rbfr(i,jmn(i,j))  &
             +u(i,j,mykm) )/2.
     enddo
  enddo

  !      return
end subroutine enfz
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine initialize
  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  real :: dum,dum1,dum2,jacp,xndum,r,wx0,wx1
  !        complex(8),dimension(0:1) :: x,y
  real,dimension(0:1) :: x,y
  integer :: n,i,j,k,ip

  call ppinit(MyId,numprocs,ntube,TUBE_COMM,GRID_COMM)
  ! write(*,*)'ppinit  ',myid,numprocs,ntube,TUBE_COMM,GRID_COMM

  !     program begins....

  !     reset timestep counter.
  Last=numprocs-1
  timestep=0
  tcurr = 0.

  !     read input data and initialize...

  call hybinit
  do i=0,Last
     if (MyId.eq.i) call init
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  enddo
  if(idg.eq.1)write(*,*)'past init'

  dum = 0.
  dum1 = 0.
  do i = 0,im-1
     r = xg(i)-0.5*lx+lr0
     j = int((r-rin)/dr)
     j = min(j,nr-1)
     wx0 = (rin+(j+1)*dr-r)/dr
     wx1 = 1.-wx0
     xndum = wx0*xn0e(j)+wx1*xn0e(j+1)
     dum = dum+(jac(i,0)+jac(i+1,0)+jac(i,1)+jac(i+1,1))/4
     dum1 = dum1+xndum*(jac(i,0)+jac(i+1,0)+jac(i,1)+jac(i+1,1))/4
  end do
  call MPI_ALLREDUCE(dum,jacp,1,  &
       MPI_REAL8,MPI_SUM,           &
       tube_comm,ierr)

  call MPI_ALLREDUCE(dum1,dum2,1,  &
       MPI_REAL8,MPI_SUM,           &
       tube_comm,ierr)

  totvol = dx*ly*dz*jacp    
  !!!new tmm(1) definition (htc)
  n0=(real(tmm(1))/totvol)*real(ntube)
  n0e=mme*numprocs/totvol
  if(myid==0)then
     write(*,*)'totvol,jacp,dum2=',totvol,jacp,dum2
  end if

  !  calculate the volume of each radial subdomain
  do k = 1,nsubd
     dum = 0.
     do i = (k-1)*im/nsubd,k*im/nsubd-1
        r = xg(i)-0.5*lx+lr0
        j = int((r-rin)/dr)
        j = min(j,nr-1)
        wx0 = (rin+(j+1)*dr-r)/dr
        wx1 = 1.-wx0
        dum = dum+(jac(i,0)+jac(i+1,0)+jac(i,1)+jac(i+1,1))/4
     end do
     call MPI_ALLREDUCE(dum,jacp,1,  &
          MPI_REAL8,MPI_SUM,           &
          tube_comm,ierr)
     vol(k) = dx*ly*dz*jacp    
  end do
  call weight
  !     initialize particle quantities...
  if( cut.eq.0.) cut=1.
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  call ccfft('x',0,imx,0.0,tmpx,coefx,workx,0)
  call ccfft('y',0,jmx,0.0,tmpy,coefy,worky,0)
  call ccfft('z',0,kmx,0.0,tmpz,coefz,workz,0)
  call dsinf(1,x,1,0,1,0,imx*2,1,1.0,aux1,50000,aux2,20000)

  ncurr = 1
  call blendf
end subroutine initialize
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine loader_wrapper
  use gem_com
  use gem_equil
  implicit none

  integer :: n,i,j,k,ip,ns

  do ns = 1,nsm
     if(isuni.eq.0)call loadi(ns)
  end do
  if(idg.eq.1)write(*,*)'past loader'
end subroutine loader_wrapper
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine accumulate(n,ip)
  use gem_com
  use gem_equil
  implicit none

  integer :: n,i,j,k,ip
  grid1_start_tm=grid1_start_tm+MPI_WTIME()
  call grid1(ip,n)
  grid1_end_tm=grid1_end_tm+MPI_WTIME()

  if(idg.eq.1)write(*,*)'pass grid1'
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine accumulate
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine poisson(n,ip)
  use gem_com
  use gem_equil
  implicit none
  integer :: n,i,i1,j,k,ip,it,iter=1
  real :: myrmsphi,rmp(20),myavap(0:imx-1)
  real :: myjaca(0:imx-1),jaca(0:imx-1)

  do it = 1,iter
     call gkps_adiabatic_electron(n,ip)

     if(idg.eq.1)write(*,*)'pass gkps in poisson'
     myrmsphi=0.
     rmp(it)=0.
     do k=0,mykm-1
        do j=0,jm-1
           do i1=0,im-1
              myrmsphi=myrmsphi+phi(i1,j,k)*phi(i1,j,k)
           enddo
        enddo
     enddo
     call MPI_ALLREDUCE(myrmsphi,rmp(it),1, &
          MPI_REAL8,                               &
          MPI_SUM,TUBE_COMM,ierr)
     rmp(it)=sqrt(rmp(it)/(im*jm*km))
  end do
  if(idg.eq.1)write(*,*)'pass iter loop in poisson'
  rmsphi(n)=rmp(iter)        
  if(ip==1)ipred = 1
  if(ip==0)icorr = 1
  if(iter==1)goto 100
  if(rmp(iter)/rmp(iter-1)>1.1)then
     phi(:,:,:) = 0.
     if(ip==1)ipred = 0
     if(ip==0)icorr = 0
     rmsphi(n)=0
  end if
  if(n>100.and.rmp(iter)/rmpp>1.5)then
     phi(:,:,:) = 0.
     if(ip==1)ipred = -1
     if(ip==0)icorr = -1
     rmsphi(n)=0           
  end if
  rmpp = rmp(iter)
100 continue


  !remove zonal field
  if(izonal==0)then
     do i=0,im-1
        myavap(i) = 0.
        myjaca(i) = 0.
        do j=0,jm-1
           myavap(i)=myavap(i)+phi(i,j,0)*jac(i,0)
           myjaca(i)=myjaca(i)+jac(i,0)
        enddo
     enddo
     call MPI_ALLREDUCE(myavap,avap,imx, &
          MPI_REAL8,                               &
          MPI_SUM,TUBE_COMM,ierr)
     call MPI_ALLREDUCE(myjaca,jaca,imx, &
          MPI_REAL8,                               &
          MPI_SUM,TUBE_COMM,ierr)

     avap(0:imx-1)=avap(0:imx-1)/jaca(0:imx-1)

     do i = 0,imx-1
        do j = 0,jmx
           do k = 0,1
              phi(i,j,k) = phi(i,j,k)-avap(i)
           end do
        end do
     end do
  end if

  myrmsphi=0.
  do k=0,mykm-1
     do j=0,jm-1
        do i1=0,im-1
           myrmsphi=myrmsphi+phi(i1,j,k)*phi(i1,j,k)
        enddo
     enddo
  enddo
  call MPI_ALLREDUCE(myrmsphi,rmsphi(n),1, &
       MPI_REAL8,                               &
       MPI_SUM,TUBE_COMM,ierr)
  rmsphi(n)=sqrt(rmsphi(n)/(im*jm*km))

  if(idg.eq.1)write(*,*)'pass poisson'
  !$acc update device(phi)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine poisson
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine field(n,ip)
  use gem_com
  use gem_equil
  implicit none
  integer :: n,i,j,k,ip
  call grad(ip)
  if(idg.eq.1)write(*,*)'pass grad'
  call eqmo(ip)
  if(idg.eq.1)write(*,*)'pass eqmo'
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !$acc update device(phi,ex,ey,ez)        
end subroutine field
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine push_wrapper(n,ip)
  use gem_com
  use gem_equil
  implicit none
  integer :: n,i,j,k,ip,ns

  do ns = 1,nsm
     ppush_start_tm=ppush_start_tm+MPI_WTIME()
     if(ip.eq.1.and.ision==1)call ppush(n,ns)
     ppush_end_tm=ppush_end_tm+MPI_WTIME()
     cpush_start_tm=cpush_start_tm+MPI_WTIME()
     if(ip.eq.0.and.ision==1)call cpush(n,ns)
     cpush_end_tm=cpush_end_tm+MPI_WTIME()
  end do

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)          
end subroutine push_wrapper
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine diagnose(n)
  use gem_com
  use gem_equil
  implicit none
  integer :: n,i,j,k,ip

  call modes2(phi,pmodehis,n)
  if(idg.eq.1)write(*,*)'pass modes2'  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)          
  if(idg.eq.1)write(*,*)'pass yvec'    
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)        

end subroutine diagnose
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine reporter(n)
  use gem_com
  use gem_equil
  implicit none
  integer :: n,i,j,k,ip

  if(modulo(n,xnplt).eq.0) then
     call spec(n)
  endif
13 format(1x,i6,7(2x,i7))
  if(myid.eq.master)then
     open(16, file='indicator', status='unknown',position='append')
     write(16,13)n,ipred,icorr,jpred,jcorr,nopz,noen,nowe
     close(16)
  end if
  !        if(myid==0)write(*,13)n,ipred,icorr,jpred,jcorr,nopz,noen
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !     save particle arrays for restart if iput=1...
  !     do this before the code crashes due to graphics problems
  if((iput.eq.1).and.modulo(n+1,500).eq.0)call restart(2,n)

  !     periodically make output for plots
  call outd(n)
  if(idg.eq.1)write(*,*)'pass outd'

end subroutine reporter
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine dcmpy(u,v)
  use gem_com
  use gem_fft_wrapper
  implicit none
  INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep,ip,INFO
  INTEGER :: l1,m1,myk,myj,ix,ikx,id
  INTEGER :: recvcnt(0:ntube-1)
  real :: u(0:imx-1,0:jmx-1,0:1)
  complex :: v(0:imx-1,0:jcnt-1,0:1),myv(0:imx-1,0:jcnt-1,0:1)
  complex :: sbuf(0:imx*jmx*2-1),rbuf(0:imx*jcnt*2-1)
  complex :: temp3d(0:imx-1,0:jmx-1,0:1)
  real :: kx,ky,kx0,th,shat,sgny

  do k=0,mykm
     do j=0,jm-1
        do i=0,imx-1
           temp3d(i,j,k)=u(i,j,k)
        enddo
     enddo
     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = temp3d(i,j,k)
        end do
        call ccfft('y',-1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3d(i,j,k) = tmpy(j)
        end do
     end do
  enddo

  do m = 0,jcnt-1
     myv(:,m,:) = temp3d(:,jft(m),:)
  end do

  cnt = 2*jcnt*imx
  call mpi_allreduce(myv,v,cnt,MPI_DOUBLE_COMPLEX,mpi_sum, &
       grid_comm,ierr)

  !      return
end subroutine dcmpy
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
real function en3(s)
  real :: s

  if(0 .le. s .and. s .le. 1.0)then
     en3 = s*s/2
  else if(1.0 .le. s .and. s .le. 2.0) then
     en3 = -3./2.+3*s-s*s
  else if(2.0 .le. s .and. s .le. 3.0) then
     en3 = (3-s)*(3-s)/2.0
  else
     en3 = 0.
  end if

  return
end function en3
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine blendf
  use gem_com
  use gem_equil
  implicit none

  INTEGER :: m,i,j,k,m1,m2
  complex :: dum,dum1
  complex :: tmp(nb,nb),work(100)
  integer :: IPIV(10),INFO
  real :: r,qr,s1,s2,s3,dth1,wx0,wx1
  real :: aky(0:jmx-1),dely(0:imx-1)

  dth1 = pi2/nb
  do j = 0,im-1
     r = rin+xg(j)
     i = int((r-rin)/dr)
     i = min(i,nr-1)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     dely(j) = modulo(-pi2*lr0/q0*qr*sign(1.0,q0)+8000.*ly,ly)*tor
  enddo
  do j = 0,jm-1
     if(j.ge.(jm/2+1)) then
        aky(j) = -2.*pi*real(jm-j)/ly
     else
        aky(j) = 2.*pi*real(j)/ly
     end if
  enddo
  do j = 0,jm-1
     do i = 0,im-1
        pfac(i,j) = exp(IU*aky(j)*dely(i))
     end do
  enddo

  do m = 1,nb
     do i = 0,imx-1
        do j = 0,jmx-1
           do k = 0,kmx
              s1 = real(k)*pi2/kmx/dth1-m
              s2 = real(k)*pi2/kmx/dth1-m-nb
              s3 = real(k)*pi2/kmx/dth1-m+nb
              pol(m,i,j,k) = en3(s1)+en3(s2)*pfac(i,j)+en3(s3)/pfac(i,j)
           end do
        end do
     end do
  end do

  do i = 0,imx-1
     do j = 0,jmx-1
        do m1 = 1,nb
           do m2 = 1,nb
              dum = 0.
              do k = 0,km-1
                 dum = dum+conjg(pol(m1,i,j,k))*pol(m2,i,j,k)
              end do
              tmp(m1,m2) = dum
              pmtrx(i,j,m1,m2) = dum
           end do
        end do
        call ZGETRF(nb,nb,tmp,nb,IPIV,INFO)
        call ZGETRI(nb,tmp,nb,IPIV,work,100,INFO)
        do m1 = 1,nb
           do m2 = 1,nb
              pmtrxi(i,j,m1,m2) = tmp(m1,m2)
           end do
        end do
     end do
  end do

  return
  if(myid==0)then
     i = 5
     j = 1
     do k = 0,kmx
        write(*,10)(pol(m,i,j,k),m=1,nb)
     end do

     write(*,*)'pfac= ', pfac(i,j)
     do m = 1,nb
        write(*,10)pol(m,i,j,km)-pol(m,i,j,0)*pfac(i,j)
     end do
10   format(12(1x,e10.3))

     write(*,*)' '
     write(*,*)' '

     do m1=1,nb
        write(*,11)(pmtrx(i,j,m1,k),k=1,nb)
     end do
     write(*,*)' '
     write(*,*)' '
     do m1=1,nb
        write(*,11)(pmtrxi(i,j,m1,k),k=1,nb)
     end do

     do m1=1,nb
        do m2=1,nb
           dum = 0.
           do k = 1,nb
              dum = dum+pmtrx(i,j,m1,k)*pmtrxi(i,j,k,m2)
           end do
           tmp(m1,m2) = dum
        end do
     end do

     write(*,*)' '
     write(*,*)' '
     do m1=1,nb
        write(*,11)(tmp(m1,k),k=1,nb)
     end do
  end if

11 format(12(1x,e10.3))

  !      return
end subroutine blendf
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine filtbl(u)   
  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  complex :: lbfs(0:imx-1,0:jmx-1)
  complex :: rbfr(0:imx-1,0:jmx-1)
  complex :: u(0:imx-1,0:jmx-1,0:1)
  complex :: bm(1:nb),cm(1:nb),dum
  real :: qr
  real :: aky(0:jmx-1),dely(0:imx-1),dklz(0:imx-1,0:jmx-1)
  INTEGER :: i,j,k,l,m,n,ind,m1,m2
  INTEGER :: myk,mynum,id
  complex,dimension(:),allocatable :: holdu,rbuf,sbuf

  !      write(*,*)'enter filtor'
  mynum = imx*jmx/numprocs      
  allocate(holdu(0:imx*jmx*kmx/numprocs-1),  &
       rbuf(0:mynum-1), &
       sbuf(0:mynum-1))

  !      return
  !pack data to send

  do id = 0,last
     if(id.ne.myid)then
        do ind = 0,mynum-1 !id*mynum,(id+1)*mynum
           j = (id*mynum+ind)/im
           i = id*mynum+ind-j*im
           sbuf(ind) = u(i,j,0)
        end do
        !             send sbuf to id
        call MPI_SENDRECV(sbuf(0),mynum, &
             mpi_double_complex,id,10, &
             rbuf(0),mynum, &
             mpi_double_complex,id,10, &
             mpi_comm_world,stat,ierr)
        !unpack and put into holdu
        do ind = 0,mynum-1
           holdu(ind*km+int(id/ntube)) = rbuf(ind)
        end do
     end if
  end do

  !put own u(:,:,) in holdu
  do ind = 0,mynum-1
     j = (myid*mynum+ind)/im
     i = myid*mynum+ind-j*im
     holdu(ind*km+int(myid/ntube)) = u(i,j,0)
  end do

  do ind = 0,mynum-1
     j = (myid*mynum+ind)/im
     i = myid*mynum+ind-j*im
     do k = 0,km-1
        tmpz(k) = holdu(ind*km+k)
     end do

     !do blending function filtering
     do m1 = 1,nb
        dum = 0.
        do k = 0,km-1
           dum = dum+conjg(pol(m1,i,j,k))*tmpz(k)
        end do
        bm(m1) = dum
     end do

     do m1 = 1,nb
        dum = 0.
        do m2 = 1,nb
           dum = dum+pmtrxi(i,j,m1,m2)*bm(m2)
        end do
        cm(m1) = dum
     end do

     do k = 0,km-1
        dum = 0.
        do m1 = 1,nb
           dum = dum+cm(m1)*pol(m1,i,j,k)
        end do
        tmpz(k) = dum
     end do

     do k = 0,km-1
        holdu(ind*km+k) = tmpz(k)  
     end do
  end do

  do id = 0,last
     if(id.ne.myid)then
        do ind = 0,mynum-1
           sbuf(ind) = holdu(ind*km+int(id/ntube))
        end do
        !             send sbuf to id
        call MPI_SENDRECV(sbuf(0),mynum, &
             mpi_double_complex,id,20, &
             rbuf(0),mynum, &
             mpi_double_complex,id,20, &
             mpi_comm_world,stat,ierr)
        !unpack and put into u
        do ind = 0,mynum-1
           j = (id*mynum+ind)/im
           i = id*mynum+ind-j*im
           u(i,j,0) = rbuf(ind)
        end do
     end if
  end do

  !put own holdu in u
  do ind = 0,mynum-1
     j = (myid*mynum+ind)/im
     i = myid*mynum+ind-j*im
     u(i,j,0) = holdu(ind*km+int(myid/ntube))
  end do

  !       write(*,*)'before assign',myid
  !assign u(:,:,1)

  do i = 0,im-1
     do j = 0,jm-1
        lbfs(i,j) = u(i,j,0)
     end do
  end do

  call MPI_SENDRECV(lbfs(0:imx-1,0:jmx-1),imx*jmx, &
       MPI_double_complex,lngbr,10, &
       rbfr(0:imx-1,0:jmx-1),imx*jmx, &
       MPI_double_complex,rngbr,10, &
       TUBE_COMM,stat,ierr)
  !       write(*,*)'after send_recv',myid
  if(gclr.ne.glst)then                  
     do i = 0,im-1
        do j = 0,jm-1
           u(i,j,1) = rbfr(i,j)
        end do
     end do
  end if
  if(gclr==glst)then
     do i = 0,im-1
        do j = 0,jm-1
           u(i,j,1) = rbfr(i,j)*pfac(i,j)
        end do
     end do
  end if

  !       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !       return
end subroutine filtbl
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine gam(u,v)
  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  !     
  complex :: u(0:imx-1),v(0:imx-1)
  INTEGER :: n,i,j,k,k1,i1

  call ccfft('z',0,kmx,1.0,tmpz,coefz,workz,0)

  do i1 = 0,imx-1
     tmpz(0) = u(i1)

     if(GCLR.ne.master)then
        call MPI_SEND(tmpz(0),mykm,MPI_DOUBLE_COMPLEX,master, &
             gclr,tube_comm,ierr)
     end if

     if(gclr.eq.master)then
        do i = 1,GLST
           call MPI_RECV(tmpz(i*mykm),mykm,MPI_DOUBLE_COMPLEX,i, &
                i,tube_comm,stat,ierr)
        end do
     end if

     if(GCLR.eq.master) then
        call ccfft('z',1,kmx,1.0,tmpz,coefz,workz,0)
     end if

     call MPI_BCAST(tmpz,kmx,MPI_DOUBLE_COMPLEX,master, &
          tube_comm,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     do k = 0,kmx-1
        k1 = k
        if(k>kmx/2)k1=kmx-k
        if(k1>0)tmpz(k) = 0.
     end do
     call ccfft('z',-1,kmx,1.0,tmpz,coefz,workz,0)
     tmpz = tmpz/kmx

     v(i1) = tmpz(gclr)
  end do

  !      return
end subroutine gam
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ftcamp
  use gem_com
  use gem_fft_wrapper
  implicit none

  complex :: mycampf(0:mynf-1)
  real :: aomega(0:nfreq-1)
  integer :: i,j,k,thisf,nsize,ir
  real :: domega,om,dum1,dum2,dum3,dum4,dum5,x,gam,peak(0:6,1:5)
  complex :: sbuf(0:mynf-1),rbuf(0:nfreq-1)

  if(idg==1)write(*,*)'enter ftcamp'
  nsize=nm/ifskp
  domega = 2*frmax/nfreq
  do i = 0,nfreq-1
     aomega(i) = -frmax+i*domega
  end do
  !      ir = irlk
  do ir = 0,6
     !make the data stationary
     dum1 = 0.
     do i = 100,500
        dum1 = dum1+abs(camp(ir,i))
     end do
     dum1 = dum1/400
     dum2 = 0.
     do i = nsize-401,nsize-1
        dum2 = dum2+abs(camp(ir,i))
     end do
     dum2 = dum2/400
     gam = log(dum2/dum1)/(nsize-400)
     !!!a meaningless loop (htc)
     !do i = 0,nsize-1
        !         camp(:,i) = camp(:,i)/exp(gam*i)
        !         camp(:,i) = exp(IU*1.5e-3*dt*ifskp*i)*exp(gam*i)  !test FT effect
     !end do

     do j = 0,nfreq-1
        om = aomega(j)
        campf(ir,j) = 0.
        do i = 0,nsize-1
           campf(ir,j) = campf(ir,j)+camp(ir,i)*exp(-IU*om*dt*ifskp*i)
        end do
     end do

     !find 5 peaks
     dum1 = 0.
     do i=0,nfreq-1
        x = abs(campf(ir,i))
        if(x>dum1)then
           dum1 = x
        end if
     end do
     peak(ir,1)=dum1
     dum2 = 0.
     do i=0,nfreq-1
        x = abs(campf(ir,i))
        if(x>dum2.and. x .ne. dum1)then
           dum2 = x
        end if
     end do
     peak(ir,2)=dum2
     dum3 = 0.
     do i=0,nfreq-1
        x = abs(campf(ir,i))
        if(x>dum3 .and. x.ne.dum1 .and. x.ne.dum2)then
           dum3 = x
        end if
     end do
     peak(ir,3)=dum3
     dum4 = 0.
     do i=0,nfreq-1
        x = abs(campf(ir,i))
        if(x>dum4 .and. x.ne.dum1 .and. x.ne.dum2 .and. x.ne.dum3)then
           dum4 = x
        end if
     end do
     peak(ir,4)=dum4
     dum5 = 0.
     do i=0,nfreq-1
        x = abs(campf(ir,i))
        if(x>dum5 .and. x.ne.dum1 .and. x.ne.dum2 .and. x.ne.dum3 .and. x.ne.dum4)then
           dum5 = x
        end if
     end do
     peak(ir,5)=dum5
  end do

  do k = 0,glst
     if(tclr==0 .and. gclr==k)then
        open(15, file='freq', status='unknown',position='append')
        do ir=0,6
           do i = 0,nfreq-1
              x = abs(campf(ir,i))
              if(x==peak(ir,1)) write(15,10)i,aomega(i),abs(campf(ir,i))**2
              if(x==peak(ir,2)) write(15,10)i,aomega(i),abs(campf(ir,i))**2
              if(x==peak(ir,3)) write(15,10)i,aomega(i),abs(campf(ir,i))**2
              if(x==peak(ir,4)) write(15,10)i,aomega(i),abs(campf(ir,i))**2
              if(x==peak(ir,5)) write(15,10)i,aomega(i),abs(campf(ir,i))**2
           end do
           write(15,*)'frequency, gclr, ir=', gclr,ir
           do i = 0,nfreq-1
              write(15,10)i,aomega(i),abs(campf(ir,i))**2
           end do
           write(15,*)'end ir=',ir
        end do
        close(15)
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
10   format(1x,i6, 3(2x,e12.5))
  end do
  close(15)
  !      return
end subroutine ftcamp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

