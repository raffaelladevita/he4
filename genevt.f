c***************************************************************************

	subroutine genevt(w,jch)
  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:					    c
c 									    c
c  this subroutine extracts theta and phi angles of two-body decay 
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
c  113) ppbar      decay
c  114) pp      decay
c  115) zeta+      decay
c--------------------------------------------------------------------------c

	implicit none
	integer nchan,nchan1,nchan2
	parameter (nchan  = 100 )
	parameter (nchan1 = 101 )
	parameter (nchan2 = 200 )
	real pigr,twopigr,deltang
	real degrad
	real rran,dsdt_valery,dsdt_max,dsdt_min,t_random_exp
	real rmass
	integer np,kgenev,ij
	real ecm,amass,pcm,wt
        common/genin/np,ecm,amass(18),kgenev
        common/genout/pcm(5,18),wt
	character*12 fragment(10)
	common/fragm/fragment
	real rm0,rm1,rm2,tth,pph
	real t_ref,target_p(4),proton_diff(4),t_mandel,scal_prod,extr_mode
	common /t_refer/ t_ref,extr_mode
        integer nchain,ichain(nchan)
        common /channels/nchain,ichain
	real vecgen(0:nchan2)
	common/geneven/vecgen
	real t_org,w_org
	real fit_esp,fit_lstar_xsec
	common/original/t_org,w_org
	logical ds_lambdastar_fn,ds_lambdastar_fn_new

c        *** aggiunta marco ***

c    definisco in common la flag per l estrazione della omega 

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
	integer jch,j,k,i,ii,kk,iw
	real ds_fn,dsdt_fn
	real Eg_x,Eg_cm_x,Eo_cm_x,po_cm_x,t_max_x,t_min_x,t_extr,yyy	
	real Eg_cm,E_cm,p_cm,costheta_cm,t
c	logical first_adam
c	data first_adam/.true./

	double precision adam_e,adam_t,cttheta,ttheta,pphi,wm2p,wm2k
	double precision cross_pwa
	double precision W_ad,Eg_ad,E_ad,p_ad,cost2pCM,thcm2p
	common/dipion/adam_e,adam_t,cttheta,ttheta,pphi,wm2p,wm2k


        data pigr/3.141592/
	data twopigr/6.283184/
	data deltang/0.157080/
        data degrad/0.0175/
	

	jord=0

c	print *,jch
c Adam Scz... 2 k PWA generator
	if(jch.eq.83) then
c	   if(first_adam) then
c	      print *,' Initializing Adam 2pion event generator'
c	      adam_e=3.5
c	      adam_t=0.45
c	      call init(adam_e,adam_t)
c	      first_adam=.false.
c	   endif
c+ 0) Extracting  2pi system mass M2pi, th and ph decay
c+ cross_pwa is normalized (/0.07) only for e=3.5 and -t=0.45
c 565	   wm2p=(1.4-0.4)*rran()+0.4
c         cttheta=(1+1)*rran()-1.
c	 ttheta=acos(cttheta)
c	 pphi=rran()*twopigr
c	 yyy=rran()
c	 print *, 'ok'
c         if (yyy.gt.cross_pwa(3.5,0.45,wm2p,cttheta,pphi)/0.07) goto 565
c	 print*,  yyy,cross_pwa(3.5,0.45,wm2p,cttheta,pphi)/0.07,wm2p,cttheta,pphi
c	 goto 565
c+ 1) Given 2pi system mass (wm2p) gamma energy (adam_e) and -t (adam_t) calculating W, p_cm, thcm2pi
c	 w=sqrt(2*amass(2)*adam_e+amass(2)**2)
	 Eg_ad=(w**2-amass(2)**2)/2/w
	 E_ad=(w**2+wm2k**2-amass(2)**2)/2/w
	 p_ad=sqrt(E_ad**2-wm2k**2)
	 cost2pCM=(-adam_t-wm2k**2+2*Eg_ad*E_ad)/(2*Eg_ad*p_ad)
c	 print*,  wm2p,w,eg_ad,e_ad,p_ad,acos(cost2pCM)/degrad
c	 print*,ttheta,pphi


	 theta_m=acos(cost2pCM)
	 phi_m = rran()*twopigr
	 etot_m = E_ad
	 p_m=p_ad
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
	   goto 99
	endif


c Adam Scz... 2 pion PWA generator
	if(jch.eq.87) then
c	   if(first_adam) then
c	      print *,' Initializing Adam 2pion event generator'
c	      adam_e=3.5
c	      adam_t=0.45
c	      call init(adam_e,adam_t)
c	      first_adam=.false.
c	   endif
c+ 0) Extracting  2pi system mass M2pi, th and ph decay
c+ cross_pwa is normalized (/0.07) only for e=3.5 and -t=0.45
c 565	   wm2p=(1.4-0.4)*rran()+0.4
c         cttheta=(1+1)*rran()-1.
c	 ttheta=acos(cttheta)
c	 pphi=rran()*twopigr
c	 yyy=rran()
c	 print *, 'ok'
c         if (yyy.gt.cross_pwa(3.5,0.45,wm2p,cttheta,pphi)/0.07) goto 565
c	 print*,  yyy,cross_pwa(3.5,0.45,wm2p,cttheta,pphi)/0.07,wm2p,cttheta,pphi
c	 goto 565
c+ 1) Given 2pi system mass (wm2p) gamma energy (adam_e) and -t (adam_t) calculating W, p_cm, thcm2pi
c	 w=sqrt(2*amass(2)*adam_e+amass(2)**2)

	 Eg_ad=(w**2-amass(2)**2)/2/w
	 E_ad=(w**2+wm2p**2-amass(2)**2)/2/w
	 p_ad=sqrt(E_ad**2-wm2p**2)
	 cost2pCM=(-adam_t-wm2p**2+2*Eg_ad*E_ad)/(2*Eg_ad*p_ad)
c	 print*,  wm2p,w,eg_ad,e_ad,p_ad,acos(cost2pCM)/degrad
c	 print*,ttheta,pphi


	 theta_m=acos(cost2pCM)
	 phi_m = pphi
	 etot_m = E_ad
	 p_m=p_ad
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
	   goto 99
	endif




c	if(jch.eq.87) then
c	   wm2p=0.5
c	  do iw=0,50
c	     wm2p=(1.4-0.4)*iw/50.+0.4
c	   do ii=1,180,9
c	    do kk=1,360,9  
c	     cttheta=cos(1.*ii*degrad)
c	     pphi=1.*kk*degrad
c             print*, wm2p,cttheta,pphi,cross_pwa(3.5,0.45,wm2p,cttheta,pphi)/0.07
c             print*, wm2p,cttheta,pphi,cross_pwa(3.5,0.55,wm2p,cttheta,pphi)/91000
c	     W_ad=sqrt(2*amass(2)*adam_e+amass(2)**2)
c             Eg_ad=(W_ad**2-amass(2)**2)/2/W_ad
c	     E_ad=(W_ad**2+wm2p**2-amass(2)**2)/2/W_ad
c	     p_ad=sqrt(E_ad**2-wm2p**2)
c	     cost2pCM=(-adam_t-wm2p**2+2*Eg_ad*E_ad)/(2*Eg_ad*p_ad)
c	     print*,  wm2p,W_ad,eg_ad,e_ad,p_ad,acos(cost2pCM)/degrad
c	     thcm2p=acos(cost2pCM)/degrad

c+ Norm factor=0.07 (safer) for Eg=3.5 t=0.45  
c+ Norm factor=120000 (safer) for Eg=3.5 t=0.55 lmax8-UNNORM   
c	    enddo 
c	   enddo
c	  enddo
c	endif 


c       This is the special case for pi0 pi0 (84) and pi0 eta (85)
c	if((jch.eq.84.or.jch.eq.85).and.vecgen(jch).eq.1) then
	if((jch.eq.84.or.jch.eq.85)) then
c+ 1) Given s (w*w) calculating gamma energy in CM, proton energy in CM, p_cm, thcm2
	 Eg_cm=(w**2-amass(2)**2)/2/w
	 E_cm=(w**2+amass(1)**2-amass(2)**2)/2/w
	 p_cm=sqrt(E_cm**2-amass(1)**2)
	 CALL T_MIN_MAX(t_min_x,t_max_x,0.,amass(2)**2,amass(1)**2,amass(2)**2,w**2)
	 t=t_random_exp(3.0,t_min_x,t_max_x)
c       previous instruction returns t positive, as positive are t_min_x and t_max_x. In the formula above, note the -t
	 costheta_cm=(-t-amass(1)**2+2*Eg_cm*E_cm)/(2*Eg_cm*p_cm)
c	 print*,  wm2p,w,eg_ad,e_ad,p_ad,acos(cost2pCM)/degrad
c	 print*,ttheta,pphi


	 theta_m=acos(costheta_cm)
	 phi_m = rran()*2*pigr
	 etot_m = E_cm
	 p_m=p_cm
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
	   goto 99
	endif



c    This is the special case for single eta (32)
c     if((jch.eq.32).and.vecgen(jch).eq.1) then
        if(jch.eq.32) then
c+ 1) Given s (w*w) calculate gamma energy in cm, proton energy in cm, momentum in cm
         Eg_cm=(w**2-amass(2)**2)/2/w
         E_cm=(w**2+amass(1)**2-amass(2)**2)/2/w
         p_cm=sqrt(E_cm**2-amass(1)**2)
c        print *,amass(2),' ',amass(1), ' ',Eg_cm,' ',E_cm,' ',p_cm
         CALL T_MIN_MAX(t_min_x,t_max_x,0.,amass(2)**2,amass(1)**2,amass(2)**2,w**2)
         t=t_random_exp(3.5,t_min_x,t_max_x)
c        previous instruction returns t positive, as positive are t_min_x and t_max_x. In the formula above, note the -t
         costheta_cm=(-t-amass(1)**2+2*Eg_cm*E_cm)/(2*Eg_cm*p_cm)


         theta_m=acos(costheta_cm)
         phi_m = rran()*2*pigr
         etot_m = E_cm
         p_m=p_cm
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
           goto 99
        endif





c-
c	vecgen = 0  =>  phase-space
c	vecgen = 1  =>  tabulated diff. cross section
c+ Special caselambda*
c	if(jch.eq.36.and.w.ge.2.1.and.w.le.2.8) vecgen(jch)=1
c-
c+ Special case Dipion decay
	if(jch.eq.117) vecgen(jch)=1
c-
c+ Special case Dikaon decay
	if(jch.eq.119) vecgen(jch)=1
c-
        if(vecgen(jch).eq.0.) then
c ***	vecgen = 0 .....
	   ecm = w
c	   call genbod
c	   print *, 'calling genbod',np,ecm,amass(1),amass(2)
	   if(np.ge.4.) call genbod

	   if(np.eq.3) call gen3n1bod
	   if(np.eq.2) call gen2n1bod
c	   if(np.eq.2) call genbod
c          print *, pcm(1,1), pcm(2,1),pcm(3,1), pcm(4,1), pcm(5,1)
c	   print *, pcm(1,2), pcm(2,2),pcm(3,2), pcm(4,2), pcm(5,2)
c	   call angles(0,pcm(1,2),pcm(2,2),pcm(3,2),tth,pph)    
c	   print *,'offfirst', tth,pph,pcm(5,2),pcm(4,2),sqrt(pcm(4,2)**2-pcm(5,2)**2)

c	   print *,
c	   print*, w
	   return
	endif

 334	 continue
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
c       ***	channel decay .....
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
 9	   w2=w1+dw
	   if (w.ge.w1.and.w.lt.w2) goto 1	
	   w1=w2
	   i=i+1
	   goto 9
 1	   continue
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
 2	   end do
 3	   teta = 10*(k-1)+10/(prob(k)-prob(k-1))*(yran-prob(k-1))


c Lambda* xsec from fitted function 14/3/06
c Enbedded (no external) flag to switch it  on/off
	   ds_lambdastar_fn=.true.
c choose between old function (fit of dsigna/dcostheta of 14/3/06) or new one (fit of dsigma/dt of 30/3/10)and
	   ds_lambdastar_fn_new=.true.

c	   print *, jch
	
	   if (jch.eq.36.and.ds_lambdastar_fn) then
	      if(.not.ds_lambdastar_fn_new.and.w.ge.2.1.and.w.le.2.8) then
2010	         xran=2*rran()-1
	         yran=rran()
	         if (yran.le.fit_esp(xran,w)) then
		       continue
		       goto 2003
		    else 
		       continue
		       goto 2010
		    endif	 
c    ** adesso xran e' il mio angolo teta in gradi **
2003	         teta=acos(xran)*180/3.1415
c	         print *, teta
c	         vecgen(jch)=0
	      elseif(ds_lambdastar_fn_new.and.w.ge.2.015) then		
2011	         xran=2*rran()-1
	         yran=rran()
c		print *," "
c		print *,w,xran
c	        print *,fit_lstar_xsec(xran,2.5515)
	         if (yran.le.fit_lstar_xsec(xran,w)) then
c	         if (yran.le.fit_lstar_xsec(xran,2.1526)) then
c	         if (yran.le.0.8) then
		       continue
		       goto 2004
		    else 
		       continue
		       goto 2011
		    endif	 
c    ** adesso xran e' il mio angolo teta in gradi **
2004	         teta=acos(xran)*180/3.1415
c	         print *, teta
c	         vecgen(jch)=0	
	      endif
	   endif
c Ending Lambda* xsec from fitted function	   

c             *************  aggiunta marco  ***************

c      ** caso canale della omega ed estrazione secondo funzione continua **
c	    write (*,*) jch,ds_omega_fn,w
	   if ((jch.eq.19.or.jch.eq.35).and.ds_omega_fn.eq.1.and.w.ge.2.5.and.w.le.2.9) then

c      ** debugging di entrata nella routine **
c	       write(*,*), ds_omega_fn
c	       write(*,*) 'sono nel loop'
 999	      continue
	      if (extr_mode.eq.1) then ! extracting according theta and W distribution ds_fn(w,xran)
c     ** genero numeri casuali xran(0:180) , yran (0:massimo della f(w,teta)*sin (teta)) **
 1000		 xran=rran()*180.
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
 1010		 t_extr=(rran()*(t_max_x-t_min_x))+t_min_x
c checking and  Rejecting low -t values
		 if(-t_extr.gt.t_ref) goto 1010
		 yran=rran()
		 if (yran.le.dsdt_fn(w,t_extr)) then ! dsdt_fn is the (w,-t) distribution
c	     write(*,*) t_extr,dsdt_fn(w,t_extr)
		    xran=acos((-t_extr-(rm1**2-2*Eg_cm_x*Eo_cm_x))/2/Eg_cm_x/po_cm_x)*180/3.1415
c+ Storing original Egamma and t value
		    w_org=w
		    t_org=t_extr
		    goto 1003
		 else 
		    continue
		    goto 1010
		 endif	 
	      endif
   
c    ** adesso xran e' il mio angolo teta in gradi **
 1003	      teta=xran
c	    write(*,*) 'teta =',teta
	   endif
c  ***************************************************************************


c  ***************************************************************************
c--- added  Valery's routine for zeta+kappa* angular distribution 
c--- according to Diakonov calculation  (R. De Vita, 07/19/03)
	   if(jch.eq.90) then
	      
	      rm0=rmass('proton      ')
	      Eg_x   =(w**2-rm0**2)/2/rm0
	      Eg_cm_x=(w**2-rm0**2)/2/w
	      Eo_cm_x=(w**2+rm1**2-rm2**2)/2/w
	      po_cm_x=sqrt(Eo_cm_x**2-rm1**2)	   
    
	      CALL T_MIN_MAX(t_min_x,t_max_x,0.,rm0**2,rm1**2,rm2**2,w**2)
	      dsdt_max=dsdt_valery(t_min_x+0.001,Eg_x,rm1,rm2)
	      dsdt_min=dsdt_valery(t_max_x-0.001,Eg_x,rm1,rm2)
	      
 1020	      t_extr=(rran()*(t_max_x-t_min_x))+t_min_x
	      yran = rran()
	      if(yran.le.dsdt_valery(t_extr,Eg_x,rm1,rm2)/dsdt_max) then
		 xran=acos((-t_extr-(rm1**2-2*Eg_cm_x*Eo_cm_x))/2/Eg_cm_x/po_cm_x)*180/3.1415
	      else 
		 continue
		 goto 1020
	      endif
	      teta=xran
c       print *,t_extr,dsdt_valery(t_extr,Eg_x,rm1,rm2)/dsdt_max,teta
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

	   if(jch.eq.117) then
c+ Dipion decay alreday have all variables stored
c	  print*, 'decay ', ttheta,pphi
	  teta=ttheta
	  phi=pphi
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

	   if(jch.eq.119) then
c+ Dikaon decay alreday have all variables stored
c	  print*, 'decay ', ttheta,pphi
	  teta=ttheta
	  phi=pphi
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
c	print *, w, rm1,rm2,p_m

	   goto 99
	 endif

c Omega is a special case
	   if(jch.eq.108) then
	      call omega_decay(w,pcm,w_org,t_org)
	      goto 99
	   endif
c       Zeta+ is a special case
	   if(jch.eq.115) then
c	      print *,'inside' 
	      ecm=w
	      call genbod
c	      call angles(0,pcm(1,1),pcm(2,1),pcm(3,1),teta,phi)   ! 1=neutron variable; 2=k+  
c	      print *,'inside', teta,phi,pcm(5,1),pcm(4,1)
	      return
	   endif

c	resonance decay .....
	   hmax=xpcle(jch,0,0)
	   do ith=0,20
	      do iph=0,40
		 if(xpcle(jch,ith,iph).gt.hmax) hmax=xpcle(jch,ith,iph)
	      enddo
	   enddo
 33	   h=rran()*hmax
c       teta=acos(pigr*rran())
c       cross sect is binned in teta (not cos(teta)!)  intervals....
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
c	print*, 'return'
	return
	end
	
	
	subroutine omega_decay(Momega,pcm_omg,w_o,t_o)
	implicit none
	real th,ph,rran,Momega,pcm_omg(5,18),th_dipi
	real R0_00, ReR0_10,R0_1m1, decay_distr, pigr
	real m3,mpi,mu,q,p,q2,p2,a2,yx,Qav,Tz,Tm,Tp
	real x,y,phi_all
	real p0a(3),ppa(3),pma(3),p0b(3),ppb(3),pmb(3),p0c(3),ppc(3),pmc(3),p0int(3),ppint(3),pmint(3)
	real thf,phf,th1,ph1,th2,ph2,th3,ph3,test(5),mpip,mpi0
	real par1(5),par2(5),par3(5)
	real w_o,t_o,w0,w1
	integer j,jth,jph
	real maxx
	pigr=3.14159


c+ Extracting perpendicular of (pi+ pi-) plane
c++ Defining matrix element values (from ABBHHM and SLAC data) and deacy_distr where:
c++ th = theta pi+.vector.pi- plane
c++ ph = theta pi+.vector.pi- plane

c++ Set obteined fitting g6a and g6b data
c+ Fit performed using eff derived by previous simulated distr
c+ the bin energies are:
c+ Eg = 3.197 ->  W = 2.6225
c+ Eg = 3.377 ->  W = 2.6861
c+ Eg = 3.557 ->  W = 2.7482
c+ Eg = 3.737 ->  W = 2.8090
c+ Eg = 3.917 ->  W = 2.8685
c	t_o=3.6
	
	if(w_o.gt.2.6.and.w_o.lt.2.8685.and.
     %     t_o.gt.1.) then
          PAR1(1) = 0.8093e-1
          PAR1(2) = 0.1101
          PAR1(3) = -0.2574e-1
          PAR1(4) = -0.722e-2
          PAR1(5) = 0.2109e-2

          PAR2(1) =-0.2327
          PAR2(2) = 0.3091
          PAR2(3) = -0.686e-1
          PAR2(4) =-0.2066e-1
          PAR2(5) = 0.4965e-2

          PAR3(1) =2.222
          PAR3(2) = -4.066
          PAR3(3) = 2.522
          PAR3(4) =-0.639
          PAR3(5) = 0.5638e-1
	R0_00   =
     %   par1(1)+par1(2)*t_o+par1(3)*t_o**2+par1(4)*t_o**3+par1(5)*t_o**4
	R0_1m1  =
     %   par2(1)+par2(2)*t_o+par2(3)*t_o**2+par2(4)*t_o**3+par2(5)*t_o**4
	ReR0_10   =
     %   par3(1)+par3(2)*t_o+par3(3)*t_o**2+par3(4)*t_o**3+par3(5)*t_o**4
c+ find maximum
	maxx=-1.
	do jth=0,30
	   do jph=0,90
	     th=jth*3.1416/30.
	     ph=jph*6.2832/90.
	decay_distr=3./4./pigr*(
     @                      0.5 * (1-R0_00)
     @                     +0.5 * (3*R0_00-1)*cos(th)**2
     @                     -1.414 * ReR0_10 * sin(2*th) * cos(ph)
     @                     -R0_1m1*           sin(th)**2  * cos(2*ph))
     @              *sin(th)
	 maxx=max(maxx,decay_distr)
	   enddo
	enddo
	maxx=maxx*1.5 ! safe factor (to not have fine th ph binning)
	else
	  R0_00   = 0.2
	  R0_1m1  = -0.05
	  ReR0_10 = 0.05 
	  maxx=0.112
	endif


 765	th=pigr*rran()
	ph=2.*pigr*rran()
	yx=rran()*maxx ! just a normalization 
	decay_distr=3./4./pigr*(
     @                      0.5 * (1-R0_00)
     @                     +0.5 * (3*R0_00-1)*cos(th)**2
     @                     -1.414 * ReR0_10 * sin(2*th) * cos(ph)
     @                     -R0_1m1*           sin(th)**2  * cos(2*ph))
     @              *sin(th)

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
	yx=rran()*50.
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
c+++ 1) rotating the vectors lying on the plane by theta around X(in the hel frame X=x hcm) axis
	call ROTATION(th,0.,p0b,p0int)
	call ROTATION(th,0.,ppb,ppint)
	call ROTATION(th,0.,pmb,pmint)
c	ph=30.*3.14159/180.  
c	th=0.
c+++ 2) rotating the vectors lying on the plane by phi around Z(in the hel frame Z=omega dir in hcm) axis
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

	function fit_lstar_xsec(costhcm,wmass)
c	Lambda* cross section fitted from g11 data (II iteration) 3/10
	implicit none
	real fit_lstar_xsec
	real costhcm,wmass
        real esp_x,esp_y
	real x,y
	real parxy(4)
	real parx(3)
	real pary(4)
	real mp,mk,ml	
	real wt,ekcm,eg,egt,t1,t2,tc,tmax
	real jac

	data parx/11.86,-5.97143,1.42562/
	data pary/2.09321,1.63914,0.272665,0.635543/
	data parxy/0.493391,-0.278177,-0.163617,0.088683/
	data mp,mk,ml/0.9383,0.4936,1.518/

c	limit value of W to channel threshold	
	wt=mk+ml
	if (wmass.lt.wt) then
	    eg=(wt**2-mp**2)/2/mp
	else
	    eg=(wmass**2-mp**2)/2/mp
	endif
	
c	evaluate maximum value of -t compatible with photon energy
	eg=(wmass**2-mp**2)/2/mp
	ekcm=(wmass**2+mk**2-ml**2)/2/wmass
	tc=-(2*mp*eg*(ekcm-sqrt(ekcm**2-mk**2)*costhcm)/wmass-mk**2)
	t1=-(2*mp*eg*(ekcm-sqrt(ekcm**2-mk**2))/wmass-mk**2)
	t2=-(2*mp*eg*(ekcm+sqrt(ekcm**2-mk**2))/wmass-mk**2)
	tmax=int(-(t2-t1)*10)/10.+0.1
	tmax=min(tmax,2.5)
	jac=2*mp*eg*sqrt(ekcm**2-mk**2)/wmass   ! (dt->dcostheta)

	x=-(tc-t1)/tmax
	if(x.gt.1) then
	   fit_lstar_xsec=0
           return
        endif
	y=eg


c	define t slope according to data parameterization
	esp_x=exp(parx(1)+parx(2)*y)+parx(3)
c	define energy dependence according to first -t data point
	esp_y=pary(1)*exp(-(y-pary(2))**2/2/pary(3)**2)+pary(4)
	
        fit_lstar_xsec= 
     %  exp(-x*tmax*esp_x)*esp_y
     %  +parxy(1)
     %  +parxy(2)*x*tmax+parxy(3)*y+parxy(4)*y*x*tmax
	if(fit_lstar_xsec.lt.0.01) fit_lstar_xsec=0.01	
c	print *," x=",x," y=",y,"   xsec=",fit_lstar_xsec
	fit_lstar_xsec=fit_lstar_xsec/2.8
	fit_lstar_xsec=fit_lstar_xsec*jac

	if(fit_lstar_xsec.ge.0.95) print *," x=",x," y=",y,"   xsec=",fit_lstar_xsec

	return
	end



	function fit_esp(ct,w1)
c Lambda* xsec fitted from g11-CLAS data (I iteration) 3/06
	implicit none
	real pp(7),parw1(6)
        real  fit_esp,ct,w1

	if (w1.gt.2.8) w1=2.8
	if (w1.lt.2.1) w1=2.1

c	w1=2.79
c Parw11
	pp(1)=  -0.17907     
	pp(2)=   0.82659E-01 
	pp(3)=   0.42010E-01 
	pp(4)=   0.11988E-02   
	pp(5)=  -0.86457E-02 
	pp(6)=  -0.36488E-02 
	pp(7)=   0.15912E-02 

	parw1(1)=pp(1)+pp(2)*w1+pp(3)*w1**2+pp(4)*w1**3
     %         +pp(5)*w1**4+pp(6)*w1**5+pp(7)*w1**6
c Parw12


	pp(1)=  -0.85124E-01
	pp(2)=   -3.1830    
	pp(3)=   0.77828    
	pp(4)=   0.42877       
	pp(5)=   0.10976    
	pp(6)=  -0.18774E-01
	pp(7)=  -0.25288E-01


 	parw1(2)=pp(1)+pp(2)*w1+pp(3)*w1**2+pp(4)*w1**3
     %          +pp(5)*w1**4+pp(6)*w1**5+pp(7)*w1**6


c Parw13

	pp(1)=     12.229    
	pp(2)=    -11.162    
	pp(3)=     1.2456    
	pp(4)=    0.39078    	
	pp(5)=    0.78396    
	pp(6)=   -0.45397    
	pp(7)=    0.63527E-01

 	parw1(3)=pp(1)+pp(2)*w1+pp(3)*w1**2+pp(4)*w1**3
     %         +pp(5)*w1**4+pp(6)*w1**5+pp(7)*w1**6


c Parw14
	pp(1)=    10.760    
	pp(2)=   -10.323    
	pp(3)=    1.4584    
	pp(4)=   0.34996    	
	pp(5)=   0.73816    
	pp(6)=  -0.46566     
	pp(7)=   0.70227E-01 
 	parw1(4)=pp(1)+pp(2)*w1+pp(3)*w1**2+pp(4)*w1**3
     %          +pp(5)*w1**4+pp(6)*w1**5+pp(7)*w1**6

c Parw15
	pp(1)=      13.600    
	pp(2)=     -11.692    
	pp(3)=     0.97677    
	pp(4)=     0.38411        
	pp(5)=     0.83053    
	pp(6)=    -0.43451     
	pp(7)=     0.55589E-01 
 	parw1(5)=pp(1)+pp(2)*w1+pp(3)*w1**2+pp(4)*w1**3
     %         +pp(5)*w1**4+pp(6)*w1**5+pp(7)*w1**6

c Parw16

	pp(1)=    9.9232    
	pp(2)=   -10.161    
	pp(3)=    1.5231    
	pp(4)=   0.35753    	
	pp(5)=   0.73859    
	pp(6)=  -0.46506     
	pp(7)=   0.69407E-01 
 	parw1(6)=pp(1)+pp(2)*w1+pp(3)*w1**2+pp(4)*w1**3
     %          +pp(5)*w1**4+pp(6)*w1**5+pp(7)*w1**6





        fit_esp= 
     %  parw1(1)*exp(ct*parw1(2))
     %  +parw1(3)+parw1(4)*ct+parw1(5)*ct**2+parw1(6)*ct**3

	if (fit_esp.lt.0.0036) fit_esp=0.0036

c Normalization factor (max value in the kinematic domain)
	fit_esp=fit_esp/0.19825
	return
	end







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








c     ***************************************************************
c     ***************************************************************


      REAL FUNCTION DSDT_VALERY(tt,Egamma,Mkpi,Mtheta)
	IMPLICIT NONE

	REAL tt,Egamma,Mkpi,Mtheta
	REAL t,gknz,ggvp,e2,pi,dsdtp,PS,A2,C,hc
	REAL mp,mkstar,mk0,mz,tmin,tmax,S
	REAL FF1, FF2, LAMBDA1, LAMBDA2

c      mkstar=0.892
c      mz=1.540

	mkstar=Mkpi
	mz=mtheta
	mp=0.93827
	mk0=0.4976
	S=mp**2+2.*Egamma*mp
	DSDT_VALERY=0.
	IF(SQRT(S).LE.mkstar+mz)        RETURN
	CALL T_MIN_MAX(tmin,tmax,0.,mp**2,mkstar**2,mz**2,S)
	if((tt.le.tmin) .OR. (tt.GT.tmax)) RETURN

	t=-tt
	Lambda1=1.0
	Lambda2=1.0
	FF1=(Lambda1**2-mk0**2)/(Lambda1**2-t)
	FF2=(Lambda2**2-mk0**2)/(Lambda2**2-t)
	gknz=4.
	ggvp=1.02
	pi=3.1415
	e2=4.*pi/137.035



	PS=(1./(64*pi))*1./(4.*Egamma**2*mp**2)
	C=e2*(gknz*ggvp)**2
	A2=-(t-mkstar**2)**2 * (t-(mz-mp)**2)/(t-mk0**2)**2

	DSDTP=C*PS*A2*FF1**2*FF2**2
	hc=0.389379*1000	! GeV^2 ubarn
	DSDTP=DSDTP*hc
	DSDT_VALERY=DSDTP*1000		! GeV^2 nbarns
c      print *,tt,ps,c,a2,dsdtp,tmin,tmax,S
	RETURN
      END

      SUBROUTINE T_MIN_MAX(T_MIN,T_MAX,m1_2,m2_2,m3_2,m4_2,s)
C
C---    1+2 --> 3+4
C
      IMPLICIT NONE
      REAL t_min, t_max
      REAL M1_2,M2_2,M3_2,M4_2,S

      DOUBLE PRECISION P2CM, P4CM, E2CM,E4CM, sd,tmin,tmax

      SD=S
      E2CM=(SD+M2_2-M1_2)/(2.*SQRT(SD))
      E4CM=(SD+M4_2-M3_2)/(2.*SQRT(SD))

      P2CM=SQRT(E2CM**2-M2_2)
      P4CM=SQRT(E4CM**2-M4_2)

      tmin=(m1_2 - m3_2 - m2_2 + m4_2)**2/(4*sd)-(P2CM-P4CM)**2
      tmax=(m1_2 - m3_2 - m2_2 + m4_2)**2/(4*sd)-(P2CM+P4CM)**2

      t_min=-tmin
      t_max=-tmax

      RETURN
      END

c     ******************************
c     ***************************** 
c     This function is used to extract a random number x distributed as exp(-alpha*x) between xmin and xmax
c     It uses the inversion formula
      real function t_random_exp(alpha,tmin,tmax)
        implicit none

	real rran,alpha,tmin,tmax,u,I0
c ,t_random_exp
	u=rran()
	I0=(1/alpha)*(exp(-alpha*tmin)-exp(-alpha*tmax))
	t_random_exp=(-1/alpha)*log(exp(-alpha*tmin)-alpha*u*I0)
	return
	end







