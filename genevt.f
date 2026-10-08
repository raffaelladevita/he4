c***************************************************************************

	subroutine genevt(w,jch)
  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c
c 									    c
c  this subroutine extracts teta and phi angles of two-body decay 
c  products with respect the the rest frame of the decaying particle.
c 
c  EXTRACTED kinem values
c  REFER TO MESON for: rho-N, phi-N, omega-N, eta-N, f2-N channels
c  REFER TO DELTA for: delta-N channels
c  REFER TO MESON for: delta decay
c  REFER TO POSITIVE MESON for: rho decay, phi decay 
c  AND ARE ASSIGNED TO PARTIC #1
c    
c  101) delta++ decay
c  102) delta+  decay
c  103) delta0  decay
c  104) delta-  decay
c  105) rho+    decay
c  106) rho0    decay
c  107) rho-    decay
c  108) omega     decay -> special case
c  109) phi     decay
c  110) lambda  decay
c  112) f2      decay
c--------------------------------------------------------------------------c

	implicit none
	integer nchan,nchan1,nchan2
	parameter (nchan  = 100 )
	parameter (nchan1 = 101 )
	parameter (nchan2 = 200 )
	real pigr,twopigr,deltang
	real degrad
	real rran
	integer np,kgenev,ij
	real ecm,amass,pcm,wt
        common/genin/np,ecm,amass(18),kgenev
        common/genout/pcm(5,18),wt
	character*12 fragment(10)
	common/fragm/fragment
	real rm1,rm2
	real t_ref,target_p(4),proton_diff(4),t_mandel,scal_prod,extr_mode
	common /t_refer/ t_ref,extr_mode
        integer nchain,ichain(nchan)
        common /channels/nchain,ichain
	real vecgen(0:nchan2)
	common/geneven/vecgen

c             *** aggiunta marco ***

c    definisco in common la flag per l' estrazione della omega 

	integer ds_omega_fn
        common/omega/ds_omega_fn

c ************************************************************

	real xchan(nchan,50,0:18)
	real xpcle(nchan1:nchan2,0:20,0:40)
	common /sigmas/ xchan, xpcle
	real hmax,h
	integer ith,iph
	integer jord
	real dummy1(5),dummy2(5)
	real theta_m,phi_m,etot_m,p_m,u1,u2,u3
	real sigrif(18),dtheta,phi,teta,yran,xran,x,y,w,sigriftot
	real prob(0:18),sigmodel(18),w1,w2,dw,tetaj
	integer jch,j,k,i
	real ds_fn,dsdt_fn
	real Eg_cm_x,Eo_cm_x,po_cm_x,t_max_x,t_min_x,t_extr	



        data pigr/3.141592/
	data twopigr/6.283184/
	data deltang/0.157080/
        data degrad/0.0175/
	



	jord=0

c	vecgen = 0  =>  phase-space
c	vecgen = 1  =>  tabulated diff. cross section

        if(vecgen(jch).eq.0.) then
c ***	vecgen = 0 .....
	   ecm = w
c	   call genbod
	   if(np.ne.3) call genbod
	   if(np.eq.3) call gen3n1bod
	   return
	endif

c ***	vecgen = 1 .....
c       this is for sure a two-body channel...
c       decay angles are assigned to the partic #1
	rm1=amass(1)
	rm2=amass(2)
	if(fragment(1).eq.'proton'.or.
     +     fragment(1).eq.'neutron') then
	   rm2=amass(1)
	   rm1=amass(2)
c	   change order of p.cles....
	   jord=1
	endif
c	if(fragment(2).eq.'delta++'.or.
c     +     fragment(2).eq.'delta+'.or.
c     +     fragment(2).eq.'delta0'.or.
c     +     fragment(2).eq.'delta-') then
c	   rm2=amass(1)
c	   rm1=amass(2)
c	   jord=1
c	endif

        if(jch.le.nchan.and.vecgen(jch).eq.1.) then
c ***	channel decay .....
	    do j=1,18
	        sigrif(j) = 0.
	        sigmodel(j) = 0.
	        prob(j) = 0.
	    end do
            dtheta = 10.*degrad
            phi = rran()*twopigr
	    sigriftot = 0.
	    w1=1.0
	    dw=0.05
	    i=1
 9	    w2=w1+dw
	    if (w.ge.w1.and.w.lt.w2) goto 1	
	    w1=w2
	    i=i+1
	    goto 9
1	    continue
            do j=1,18
	        tetaj = j*dtheta
                sigmodel(j) =  xchan(jch,i,j)
	        sigrif(j) = sigmodel(j)*sin(tetaj)
	        if(sigmodel(j).lt.0.) sigmodel(j)=0.
	        sigriftot = sigriftot + sigrif(j)
            end do
	    prob(0) = 0.
	    do j=1,18
	        sigrif(j) = sigrif(j)/sigriftot 
	        prob(j) = prob(j-1)+sigrif(j)
	    end do
            yran = rran()
	    do 2 k=1,18
	    if (yran.gt.prob(k-1).and.yran.le.prob(k)) goto 3
2	    end do
3	    teta = 10*(k-1)+10/(prob(k)-prob(k-1))*(yran-prob(k-1))


c             *************  aggiunta marco  ***************

c      ** caso canale della omega ed estrazione secondo funzione continua **
c	    write (*,*) jch,ds_omega_fn,w
	    if ((jch.eq.19.or.jch.eq.35).and.ds_omega_fn.eq.1.and.w.ge.2.5.and.w.le.2.9) then

c      ** debugging di entrata nella routine **
c	       write(*,*), ds_omega_fn
c	       write(*,*) 'sono nel loop'
 999	       continue
	     if (extr_mode.eq.1) then ! extracting according theta and W distribution ds_fn(w,xran)
c     ** genero numeri casuali xran(0:180) , yran (0:massimo della f(w,teta)*sin (teta)) **
1000	      xran=rran()*180.
	      yran=rran()*10.
c     **      yran va da 0. a 10. perche' ho notato che il max di ds_fn=f(w,teta)*sin(teta) 
c           varia poco con w e vale circa 10 =>lo assumo come max assoluto
c	    write(*,*) 'x=',xran
c	    write(*,*) 'y=',yran
c	    write(*,*) 'f(w,teta)*sin(teta)=',ds_fn(w,xran)
c	    write(*,*) 'w=',w
	    


c    ** scarta xran se yran > f(w,xran) e se t > -1  **  
	      if (yran.le.ds_fn(w,xran)) goto 1003
	      goto 1000
	    else if(extr_mode.eq.2) then
c	       Mo_x=rm1 ! current omega mass
c	       Mp_x=rm2 ! proton mass
	     Eg_cm_x=(w**2-rm2**2)/2/w
             Eo_cm_x=(w**2+rm1**2-rm2**2)/2/w
	     po_cm_x=sqrt(Eo_cm_x**2-rm1**2)
c            defining min and max (-t) 
	     t_max_x=-(rm1**2-2*Eg_cm_x*Eo_cm_x-2*Eg_cm_x*po_cm_x) ! 180 deg
	     t_min_x=-(rm1**2-2*Eg_cm_x*Eo_cm_x+2*Eg_cm_x*po_cm_x) ! 0 deg
 1010	     t_extr=(rran()*(t_max_x-t_min_x))+t_min_x
c checking and  Rejecting low -t values
	     if(-t_extr.gt.t_ref) goto 1010
	     yran=rran()
	      if (yran.le.dsdt_fn(w,t_extr)) then ! dsdt_fn is the (w,-t) distribution
c	     write(*,*) t_extr,dsdt_fn(w,t_extr)
	       xran=acos((-t_extr-(rm1**2-2*Eg_cm_x*Eo_cm_x))/2/Eg_cm_x/po_cm_x)*180/3.1415
               goto 1003
	      else 
		 continue
		 goto 1010
	      endif	 
	    endif
   
c    ** adesso xran e' il mio angolo teta in gradi **
1003	    teta=xran
c	    write(*,*) 'teta =',teta
	    endif
c  ***************************************************************************

            theta_m = teta*degrad
	    phi_m = phi
 	    etot_m = (w**2+rm1**2-rm2**2)/(2*w)
	    p_m=sqrt(etot_m**2-rm1**2)
	    u1=sin(theta_m)*cos(phi_m)
	    u2=sin(theta_m)*sin(phi_m)
	    u3=cos(theta_m)
	    pcm(1,1)=p_m*u1
	    pcm(2,1)=p_m*u2
	    pcm(3,1)=p_m*u3
	    pcm(4,1)=etot_m
	    pcm(5,1)=p_m
	    pcm(1,2)=-pcm(1,1)
	    pcm(2,2)=-pcm(2,1)
	    pcm(3,2)=-pcm(3,1)
	    pcm(4,2)=w-etot_m
 	    pcm(5,2)=p_m
c****************
C+ Checking t>t_ref
*********** MOVED IN MC_CLAS
c ***   target proton 4-momentum assignment in cmh
c	target_p(1) = 0
c	target_p(2) = 0
c	target_p(3) = -(w**2-0.938**2)/2./w
c	target_p(4) = sqrt(target_p(3)**2+0.938**2)
c+ t-definition in HCM
c	do ij=1,4
c	 proton_diff(ij)=pcm(ij,2)-target_p(ij)
c	print *,ij,proton_diff(ij)
c	enddo 
c++ Calculating t as scalar product:
c	t_mandel=scal_prod(proton_diff,proton_diff)
c	t_mandel=proton_diff(4)*proton_diff(4)-t_mandel
c++ Rejecting low -t values
c	print *,t_mandel
c        if(t_mandel.gt.t_ref) goto 999
c	print *,'shdkajshakjfhkj',t_mandel
c-- End t-check routine
c*******************************************
	    goto 99
	endif

        if(jch.gt.nchan.and.vecgen(jch).eq.1.) then
c Omega is a special case
	 if(jch.eq.108) then
          call omega_decay(w,pcm)
	  goto 99
	 endif
c	resonance decay .....
	     hmax=xpcle(jch,0,0)
	     do ith=0,20
		do iph=0,40
		   if(xpcle(jch,ith,iph).gt.hmax) hmax=xpcle(jch,ith,iph)
		enddo
	     enddo
 33	     h=rran()*hmax
c	     teta=acos(pigr*rran())
c	     cross sect is binned in teta (not cos(teta)!)  intervals....
	     teta=pigr*rran()
	     phi=twopigr*rran()
             ith = int(teta/deltang)
	     iph = int(phi/deltang)
	     if(h.gt.xpcle(jch,ith,iph)) goto 33
 	     etot_m = (w**2+rm1**2-rm2**2)/(2*w)
	     p_m=sqrt(etot_m**2-rm1**2)		
	     u1=sin(teta)*cos(phi)
	     u2=sin(teta)*sin(phi)
	     u3=cos(teta)
	     pcm(1,1)=p_m*u1
	     pcm(2,1)=p_m*u2
	     pcm(3,1)=p_m*u3
	     pcm(4,1)=etot_m
	     pcm(5,1)=p_m
	     pcm(1,2)=-pcm(1,1)
	     pcm(2,2)=-pcm(2,1)
	     pcm(3,2)=-pcm(3,1)
	     pcm(4,2)=w-etot_m
 	     pcm(5,2)=p_m
	     goto 99
	endif
 99	continue
	if (jord.eq.1) then
c	re-established order of p.cles....
	   do j=1,5
	     dummy1(j)= pcm(j,1)
	     dummy2(j)= pcm(j,2)
	   enddo
	   do j=1,5
	     pcm(j,1)=dummy2(j)
	     pcm(j,2)=dummy1(j)
	   enddo
	 endif

	return
	end


	subroutine omega_decay(Momega,pcm_omg)
	implicit none
	real th,ph,rran,Momega,pcm_omg(5,18),th_dipi
	real R0_00, ReR0_10,R0_1m1, decay_distr, pigr
	real m3,mpi,mu,q,p,q2,p2,a2,yx,Qav,Tz,Tm,Tp
	real x,y,phi_all
	real p0a(3),ppa(3),pma(3),p0b(3),ppb(3),pmb(3),p0c(3),ppc(3),pmc(3),p0int(3),ppint(3),pmint(3)
	real thf,phf,th1,ph1,th2,ph2,th3,ph3,test(5),mpip,mpi0
	integer j
	pigr=3.14159


c+ Extracting perpendicular of (pi+ pi-) plane
c++ Defining matrix element values (from ABBHHM and SLAC data) and deacy_distr where:
c++ th = theta pi+.vector.pi- plane
c++ ph = theta pi+.vector.pi- plane
	R0_00   = 0.2
	ReR0_10 = 0.05 
	R0_1m1  = -0.05

 765	th=pigr*rran()
	ph=2.*pigr*rran()
c+ TO CORRECT TEH THETA EXTRACTION!!!!!
	yx=rran()*0.112 ! just a normalization
	decay_distr=3./4./pigr*(
     @                      0.5 * (1-R0_00)
     @                     +0.5 * (3*R0_00-1)*cos(th)**2
     @                     -1.414 * ReR0_10 * sin(2*th) * cos(ph)
     @                     -R0_1m1*           sin(th)**2  * cos(2*ph))


	if(yx.gt.decay_distr) goto 765
c	print *, decay_distr,th,ph
c--
c-
c++ Calculating the decay matrix element and the momentum according to omega quantum number (1-)
c++ In the omega rest frame, in the 3pi plane we have pion-dipion system
c++ vec(q) = dipion CMS
c++ vec(p) = remaining pion
c++ mu=vec(q)vec(p)/q/p (angle between p/q)
c++ therefore:
c++   pi0    vec(p)
c++   pi+   -vec(p)/2+vec(q)
c++   pi-   -vec(p)/2-vec(q)
c++ m3 = 3 pion mass system (M_omega) ad Q the available energy
	mpip=0.13957
	mpi0=0.13496
c++   pi0    vec(p)
c++   pi+   -vec(p)/2+vec(q)
c++   pi-   -vec(p)/2-vec(q)


	mpi = (2*mpip + mpi0)/3
c To be checked if we can live the current omega mass Momega otherwise decooment this line
c	Momega=0.7826	
	m3 = Momega/mpi
	Qav = (m3 - 3.)*mpi 
 976	mu = rran()*2.+(-1.)
	q2 = rran()*((m3-1.)**2/4.-1.)
	p2 = ((m3+1.)**2-4.*(q2+1))*((m3-1)**2-4.*(q2+1.))/4./m3**2
c++ A2 is the sqared matrix element assuming 
c++  A= (vec(pi0)extvec(pi+)) + (vec(pi+)extvec(pi-)) + (vec(pi-)extvec(pi0))
c++ compatible with JP=1- (omega quantum numbers) and minimal Lagrangian (l and L =1)
	A2 = 9.*p2*q2*(1-mu**2)
c++ 42 was found numerically to normalize the extraction
	yx=rran()*42.
	if (yx.gt.A2) goto 976
	p=0
	if(p2.gt.0) p=sqrt(p2)
	q=0
	if(q2.gt.0) q=sqrt(q2)
c++ Calculating kinetic energy for the 3 pions (Tz->pi0, Tm->pi+ Tp->pi-)
	Tz=(sqrt(p2+1.)-1.)*mpi0
	Tm=(sqrt(p2/4+q2-p*q*mu+1.)-1.)*mpi
	Tp=(sqrt(p2/4+q2+p*q*mu+1.)-1.)*mpi
c++ Defining normalized variables
	y=Tz/Qav
	x=(Tm-Tp)/sqrt(3.)/Qav
c	print *, x,y
c++ Now I can define all pi variables in the plane (label: a)
	th_dipi=acos(mu)
	p0a(1) = p*mpi0
	p0a(2) = 0.
	p0a(3) = 0.
	ppa(1) = (-p/2+q*cos(th_dipi))*mpip
	ppa(2) = (q*sin(th_dipi))*mpip
	ppa(3) = 0.
	pma(1) = (-p/2-q*cos(th_dipi))*mpip
	pma(2) = (-q*sin(th_dipi))*mpip
	pma(3) = 0.
c++ The 3-pions configuration can rotate by an overall phi angle in this plane (label: b)
	phi_all=2*pigr*rran()
	call ROTATION(0.,phi_all,p0a,p0b)
	call ROTATION(0.,phi_all,ppa,ppb)
	call ROTATION(0.,phi_all,pma,pmb)
c	th=60.*3.14159/180.
c++ rotating according the decay matrix (label: c)
	call ROTATION(th,0.,p0b,p0int)
	call ROTATION(th,0.,ppb,ppint)
	call ROTATION(th,0.,pmb,pmint)
c	ph=30.*3.14159/180.  
c	th=0.
	call ROTATION(0.,ph,p0int,p0c)
	call ROTATION(0.,ph,ppint,ppc)
	call ROTATION(0.,ph,pmint,pmc)

c	print *, th, ph
c	print *, ppb(1),ppb(2),ppb(3)     
c	print *, pmb(1),pmb(2),pmb(3)
c	call prod_vec(ppa,pma,thf,phf)
c	print *, 'th_in ',thf
c	print *, 'ph_in  ',phf
c	call prod_vec(ppb,pmb,thf,phf)
c	print *, 'th_f ',thf
c	print *, 'ph_f  ',phf

c+ filling output variables
	do j=1,3
	 pcm_omg(j,1)=ppc(j)
	 pcm_omg(j,2)=pmc(j)
	 pcm_omg(j,3)=p0c(j)
c+ To be shure that the enrgy is conserved, we got the pi0 as a funcion of pi+/pi-
	 pcm_omg(j,3)=-ppc(j)-pmc(j)
	enddo
	 pcm_omg(4,1)=sqrt(ppc(1)**2+ppc(2)**2+ppc(3)**2+mpip**2)
	 pcm_omg(5,1)=sqrt(ppc(1)**2+ppc(2)**2+ppc(3)**2)
	 pcm_omg(4,2)=sqrt(pmc(1)**2+pmc(2)**2+pmc(3)**2+mpip**2)
	 pcm_omg(5,2)=sqrt(pmc(1)**2+pmc(2)**2+pmc(3)**2)
	 pcm_omg(4,3)=sqrt(p0c(1)**2+p0c(2)**2+p0c(3)**2+mpi0**2)
	 pcm_omg(5,3)=sqrt(p0c(1)**2+p0c(2)**2+p0c(3)**2)
c+ To be shure that the enrgy is conserved, we got the pi0 as a funcion of pi+/pi-
	 pcm_omg(4,3)=Momega-pcm_omg(4,1)-pcm_omg(4,2)
	 pcm_omg(5,3)=sqrt(pcm_omg(1,3)**2+pcm_omg(2,3)**2+pcm_omg(3,3)**2)
	 do j=1,3
	  test(j)=pcm_omg(j,1)+pcm_omg(j,2)+pcm_omg(j,3)
c	 print *, test(j)
	 enddo 
	  test(4)=pcm_omg(4,1)+pcm_omg(4,2)+pcm_omg(4,3)-Momega
c	 print *, test(4)


c	 call angles(0,ppc(1),ppc(2),ppc(3),th1,ph1)
c	 call angles(0,pmc(1),pmc(2),pmc(3),th2,ph2)
c	 call angles(0,p0c(1),p0c(2),p0c(3),th3,ph3)
c	 print *,'pip th ph  x y z e p', th1,ph1,pcm_omg(1,1),pcm_omg(2,1),pcm_omg(3,1),pcm_omg(4,1),pcm_omg(5,1)
c	 print *,'pim th ph  x y z e p', th2,ph2,pcm_omg(1,2),pcm_omg(2,2),pcm_omg(3,2),pcm_omg(4,2),pcm_omg(5,2)
c	 print *,'pi0 x y z e', pcm_omg(1,3),pcm_omg(2,3),pcm_omg(3,3),pcm_omg(4,3),pcm_omg(5,3)
c  	 call prod_vec(ppc,pmc,th1,ph1)
c	 print *, th1,ph1
c	print *, decay_distr,th/3.14159*180.,ph/3.14159*180.
	return
	end




							! *** end genevt ***



c     ******************** aggiunta marco **************************
   
        function ds_fn(w1,teta1)
	implicit none
	real ds_fn,w1,teta1
	ds_fn=sin(teta1*3.14/180)*(2.42-0.44*w1)*((-4.99+15.9*w1)*
     &  exp(-0.5*((teta1-(3.45+0.3*w1))/(39.63-9.27*w1))**2)+
     &  (0.52+0.17*w1)*exp(-0.5*((teta1-(193.6-22.4*w1))/(-256.9+116.3*w1))**2))
	return
	end
c     ***************************************************************

c     ******************** aggiunta marco **************************
   
        function dsdt_fn(w1,t1)
	implicit none
c dsdt_fn is the (w,-t) distribution
	real dsdt_fn,w1,t1,f_exp
	real ds0,ds1,ds2
	real tmin0,tmin1,tmin2
	real tmax0,tmax1,tmax2
	real dsx0,dsx1
	real a0,a1,b0,b1
	real*8 par(6),x
	logical bruno
c+ Choosing between old Brunoldi extraction or new=8/2/02 (fitted) one 
	 bruno = .false.

	if(bruno) then

  	 tmin0=0
	 tmax0=0.5
 	 a0=1
	 b0=-6.8
	 ds0=f_exp(a0,b0,t1)
	 dsx0=f_exp(a0,b0,tmax0)

	 tmin1=tmax0
	 tmax1=1.5
	 b1=-3.2
	 a1=dsx0*exp(-b1*tmin1)
	 ds1=f_exp(a1,b1,t1)
	 dsx1=f_exp(a1,b1,tmax1)  

	 tmin2=tmax1
	 tmax2=6.
	 ds2=dsx1/(5-tmax1)*(t1+5.-2*tmin2)	
 
	 if(t1.lt.tmax0) dsdt_fn=ds0
	 if(t1.gt.tmin1.and.t1.lt.tmax1) dsdt_fn=ds1
	 if(t1.gt.tmin2.and.t1.lt.tmax2) dsdt_fn=ds2
	else
c+ Fit performed using eff derived by previous simulated distr
c+ the bin energies are:
c+ Eg = 3.197 ->  W = 2.6225
c+ Eg = 3.377 ->  W = 2.6861
c+ Eg = 3.557 ->  W = 2.7482
c+ Eg = 3.737 ->  W = 2.8090
c+ Eg = 3.917 ->  W = 2.8685
	   x=t1
         dsdt_fn=0
	 if(w1.gt.2.6.and.w1.lt.2.6861) then
          PAR(1) = 21.1732
          PAR(2) = -6.11662
          PAR(3) = 0.380611
          PAR(4) = -0.279989
          PAR(5) = 0.0679245
          PAR(6) = -0.00433914
	 elseif(w1.gt.2.6861.and.w1.lt.2.7482) then
	  PAR(1) = 21.8214
	  PAR(2) = -6.11442
	  PAR(3) = 0.249341
	  PAR(4) = -0.134445
	  PAR(5) = 0.0159455
	  PAR(6) = 0.00115018
	 elseif(w1.gt.2.7482.and.w1.lt.2.8090) then
	  PAR(1) = 19.8848
	  PAR(2) = -6.07722
	  PAR(3) = 0.221753
	  PAR(4) = -0.143483
	  PAR(5) = 0.0309009
	  PAR(6) = -0.00181644
	 elseif(w1.gt.2.8090.and.w1.lt.2.87) then
	  PAR(1) = 19.5443
	  PAR(2) = -5.88904
	  PAR(3) = 0.233437
	  PAR(4) = -0.147611
	  PAR(5) = 0.0298939
	  PAR(6) = -0.00166267
	 endif   
        dsdt_fn=1./par(1)*(
     %   par(1)*exp(x*par(2))
     %  +par(3)+par(4)*x+par(5)*x**2+par(6)*x**3
c     %  +par(7)*x**4+par(8)**5
     %  )
	endif
	return
	end
c     ***************************************************************

        function f_exp(a,b,x)
	implicit none
	real f_exp,a,b,x
	f_exp=a*exp(x*b)
	return
	end














