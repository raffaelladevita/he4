controlling******************************************************************

	subroutine mc_clas
  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c
c 									    c
c    this subroutine generates the physical events.			    c
c									    c 
c    the angular distribution of outgoing particles is read in the file     c
c    he3_xdiff.mtx                                                          c
c                                                                           c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	implicit none

	integer nchan
	parameter (nchan = 100 )
	real eglab,fwhm
	real qran
	integer kk
	integer jch2(10),chan
	integer numpar
	common /beam1/eglab,fwhm
        common /beam/ ebeam,ethresh,e1tag,e2tag,fraz
        integer nchain,ichain(nchan)
	integer npcode(10)
        common /channels/nchain,ichain
	integer nevent
	real rran
	common /control/ nevent
	common/bosfile1/p,rm,egamma,q
	common/bosfile2/numpar,jch0,npcode
	common /control1/ true_evt
	integer bosout,hbookout
	common/output/bosout,hbookout
	integer ncheck
	common/check/ncheck
	real charge
	real eg_true
	integer np_code
	character*12 partic(nchan,0:8)
	character*12 pf(10)
	character*12 parent(10),daughter(10),partnew(10,10)
	integer npd
        real ebeam,ethresh,e1tag,e2tag,fraz
	integer numpart(nchan)
	common/deca/numpart,partic
	common/fragm/parent
	real true_evt
	integer np,kgenev
	real ecm,amass,pcm,wt
	character*12 target
	real bradius,tlenght,toffset
	common /target/target,bradius,tlenght,toffset
        common/genin/np,ecm,amass(18),kgenev
        common/genout/pcm(5,18),wt

	real wr(251),egl(251),sigr(nchan,251)
	common/sigtot/sigr,wr,egl

	real sigtotal(251),max_sigtotal,divdif,sighadr
	common/maxsig/sigtotal,max_sigtotal

	integer zp_final_chan,n_ks,n_kl,k0_final_chan,zp_ang_dist
	common/zetap_flags/zp_final_chan,n_ks,n_kl,k0_final_chan,zp_ang_dist

	real beta_X2cmh(3),beta_cmh2lab(3)
	real beta1(3),beta2(3)
	real p0(4),p0i(4),p0d(4,10),wd(10),rmd(10,10),bt(3,10),p0spec(4)

	integer k
	real delta(4)

	real p(5,10),rm(10),q(10)

	real vmod
	real tc
	real t_ref,extr_mode
	common /t_refer/ t_ref,extr_mode

        real plab(5,10),pn_t(4),pn_m(4),t,t1,t2
	real pcmh(5,10)
	
	real p1(4),p2(4)
	integer loop,jch,j,il,jl,jm
	integer jch0
	integer now,later,new,nd(10)
	real w,rmass,sigma,egamma,e1,e2,s,scal_prod,
     &	     totmass,yran,yval
c	integer nloop_t

	integer maxp
        integer npart
c	integer jp
	real px,py,pz,pp,ep,mp,qp,bet,code,ch
	real w0,eg,evt,teta,phi,t_mdst,M_nn,M_nnbar,t_proton,tth,pph

	logical fillbos_flag
	integer np_nd,np_dc
	real    pcm_nd(5,10),pcm_dc(5)
	real p_cm,th_cm,m_nd
	real m_dc,p_dc,th_dc,ph_dc
	integer old_ev
	data old_ev/0/

	real pigr
	data pigr/3.14159265/
      parameter (maxp = 16)
      common /cwn1/ npart,
     +              px(maxp), py(maxp),  pz(maxp), ep(maxp),
     +              pp(maxp), bet(maxp), mp(maxp),
     +              qp(maxp), code(maxp)
      common /cwncm/np_nd,m_nd(maxp),p_cm(maxp),th_cm(maxp) 
      common /cwndc/np_dc,m_dc(maxp),p_dc(maxp),th_dc(maxp),ph_dc(maxp) 
      common /cwn2/ ch, w0, eg, evt
      common /cwn3/ teta(maxp), phi(maxp)
      common /ct/ t_mdst
      common /mnnbar/ M_nn,M_nnbar
c ***   target proton 4-momentum and minimum momentum transfer tmin
	real pcm_t(4),tmin

	logical first_adam
	double precision adam_e,adam_t,cttheta,ttheta,pphi,wm2p
	real yyy,twopigr
	data first_adam/.true./
	common/dipion/adam_e,adam_t,cttheta,ttheta,pphi,wm2p
	double precision cross_pwa
	data twopigr/6.283184/




c ***   beginning of event loop
c	{
	do loop=1,nevent

	   n_ks=0
	   n_kl=0

99	 continue
	   now=0
	   later=0
	   new=0

c ***    variable initialization 
	 do il=1,6
	   do jl=1,5
	     plab(jl,il) = 0.
	     pcmh(jl,il) = 0.
	   end do
	 end do

c ***    defining evnt type: good or background (random)
c ***    defining photon energy in the bremsstrahlung spectrum
 98	 continue
         qran = rran()
	 eglab=e1tag*(e2tag/e1tag)**qran
c	  eglab=e1tag+(e2tag-e1tag)*qran
	 egamma = eglab
	 true_evt=1.
	 qran = rran()
	 if(qran.lt.fraz*0.01) then
	    continue 
	 else
	    qran=rran()
	    eglab=ethresh*(ebeam/ethresh)**qran
	    true_evt=0.
	 endif

c ***  weighting photon interaction with hadronic cross section
	 sighadr = divdif(sigtotal,egl,251,eglab,1)
	 qran=rran()
	 if(qran.gt.sighadr/max_sigtotal) goto 98
	 
c ***    defining photon energy with tagging resolution
         sigma = 0.425*fwhm*egamma/100.
         if (sigma.eq.0.) goto 11
10       continue
         e1 = 2*rran()-1
         e2 = 2*rran()-1
         s = e1*e1+e2*e2
         if (s.gt.1) goto 10
         egamma = e1*sqrt(-2*alog(s)/s)*sigma + egamma
	 
11	 continue

c ***    w and beta calculation 
c ***    choosing reaction channel
	 if (nchain.eq.1) then
	   jch0 = ichain(1)
	 else
c	print *, eglab,jch0
	   call channel(eglab,jch0)
	 end if
	 jch=jch0

c	print *, 'ok1'

c+ Dipion special case
c Adam Scz... 2 pion PWA generator
	if(jch.eq.87) then
	   if(first_adam) then
	      print *,' Initializing Adam 2pion event generator'
	      adam_e=egamma
	      adam_t=0.55
	      call init(adam_e,adam_t)
	      first_adam=.false.
	   endif
c+ 0) Extracting  2pi system mass M2pi, th and ph decay
c+ cross_pwa is normalized (/0.07) only for e=3.5 and -t=0.45
c+ Norm factor=100000 (safer) for Eg=3.5 t=0.55 lmax8-UNNORM   
 565	 wm2p=(1.4-0.4)*rran()+0.4
         cttheta=(1+1)*rran()-1.
	 ttheta=acos(cttheta)
	 pphi=rran()*twopigr
	 yyy=rran()
c	 print *, 'ok',wm2p
         if (yyy.gt.cross_pwa(adam_e,adam_t,wm2p,cttheta,pphi)/0.07) goto 565
cc         if (yyy.gt.cross_pwa(adam_e,adam_t,wm2p,cttheta,pphi)/120000.) goto 565
	 amass(1)=wm2p
	endif
c

c-

c ***    taking into account fermi motion (if target is a nucleus)
c ***    and calculating the exact w value
	 call fermi(jch,w,beta_cmh2lab,p0i)
c p0i = 4-momento del nucleone target rispetto al lab
c beta_cmh2lab e' il beta del cmh rispetto al lab
c beta_X2cmh   e' il beta di  X   rispetto al cmh
c gamma + "target" -> X -> .....
c X per ora e' proprio il cmh, pertanto
c++ Defining the spectator variables

	  p0spec(1)=-p0i(1)
	  p0spec(2)=-p0i(2)
	  p0spec(3)=-p0i(3)
	  p0spec(4)=sqrt(.938**2+p0spec(1)**2+p0spec(2)**2+p0spec(3)**2)

c--

	 do j=1,3
	    beta_X2cmh(j)=0.
	 enddo
c        w is fixed...
	 w0=w
c ***   target proton 4-momentum assignment in cmh
	pcm_t(1) = 0
	pcm_t(2) = 0
	pcm_t(3) = -(w**2-0.938**2)/2./w
	pcm_t(4) = sqrt(pcm_t(3)**2+0.938**2)

	tmin = -(0.45 + (4.14 - eglab)*0.04/0.2)

c p0 = 4-momento del cmh (4-momento target + 4-momento fascio)
	 p0(1)=p0i(1)
	 p0(2)=p0i(2)
	 p0(3)=p0i(3)+eglab
	 p0(4)=p0i(4)+eglab
 

c ***    start reaction channel....

	 np = numpart(jch)

c        decay products ("daughter") are defined: they will become
c        "parents for the next decay level...
	 do j=1,np
	    daughter(j)=partic(jch,j)
	    amass(j) = rmass(partic(jch,j))
c           mass is defined here and must not be changed in the following
c           a new mass sampling via use of function rmass (mainly if
c           resonance is broad) gives a non 4-momentum conservation
c           mass value in the following is recalled through the pcle
c           4-momentum:  sqrt(e**2 - p**2)
	 enddo

 20	 continue

c        "new" parents are "old" dauthters after a particle decay....
	 do j=1,np
	    parent(j)=daughter(j)
	 enddo

	 totmass = 0.
	 do jm=1,np 
	   totmass = totmass+amass(jm)
	 end do

 21	 continue
	 if (w.lt.totmass+0.00001) goto 99
c	   if (jch.eq.0.and.new.eq.1) print *, 'entering genevt', eglab,w,jch
c	 print*, new,w,jch
	 call genevt(w,jch)

c	 print *, 'exit genevt',w,jch

c----------------------------------------------------------------
c----    TEMPORARILY SAVING OF CM VARIABLE OF PRIMARY PARTICLES
	 if(loop.gt.old_ev) then
	    np_nd=np
	    do j=1,np
	       do il=1,5
		  pcm_nd(il,j)=pcm(il,j)
	       end do
	    end do
	    old_ev=loop
	 endif

	 if(jch.eq.115.or.(jch.eq.0.and.new.eq.1)) then ! saving Z+ 2nd part decay (usually 2=k+)
	                                                ! or second deacy part for genboded decay
	       do il=1,5
		  pcm_dc(il)=pcm(il,2)
	       end do
	  call angles(0,pcm(1,2),pcm(2,2),pcm(3,2),tth,pph)    
c	  print *,'first', tth,pph,pcm(5,2),pcm(4,2),sqrt(pcm(4,2)**2-pcm(5,2)**2)
	 endif
	 

c----------------------------------------------------------------
c ***    warning: we must boost before, and then rotate....
c ***    from beta of cm with respect to lab we define
c ***    beta1 which describes the lorentz boost
c ***    along the cm motion direction (z' axis) 
c ***    warning: we cannot do a lorentz-3d bost, because
c ***    beta and p are not given in the same ref.syst.
	 do j=1,np
	  do il=1,4
	   p1(il)=pcm(il,j)
	  end do
c ***     boosting lungo beta rispetto al cmh
	  beta1(1)=0.
	  beta1(2)=0.
	  beta1(3)=-vmod(beta_X2cmh,3)
c ***     "minus" sign here because our beta is
c ***     beta of "X" with respect to cmh
c ***     we need instead beta of cmh with respect to "X"
c ***     boost & rotation... 
	  if(beta1(3).ne.0.) then
             call lorentz(beta1,p1)
	     call rotvec(p1,beta_X2cmh)
	  endif
	  do il=1,4
	   pcmh(il,j)=p1(il)
	  end do
	  pcmh(5,j)=sqrt(scal_prod(p1,p1))
	 end do

c ***    t definition for two pi production (phase space)
c ***    check for acceptable values of t
c ***    same procedure adopted for diffractive antiproton production
        if(later.eq.0) then
	   if (jch.eq.42.or.jch.eq.43) then
c	   if ((jch.ge.15.and.jch.le.18).or.jch.eq.42.or.jch.eq.43) then
c	t calculation with spectator particle
	      do il=1,3
		 pn_t(il) = pcm_t(il) - pcm(il,3)
	      end do
	      pn_t(4) = pcm_t(4) - pcm(4,3)
c	   tc = t_ref
	      tc = 9.
	      t1 = scal_prod(pn_t,pn_t)
	      t1 = pn_t(4)*pn_t(4)-t1
c	t calculation with particle from pair
	      do il=1,3
		 pn_t(il) = pcm_t(il) - pcm(il,1)
	      end do
	      pn_t(4) = pcm_t(4) - pcm(4,1)
	      t2 = scal_prod(pn_t,pn_t)
	      t2 = pn_t(4)*pn_t(4)-t2
	      yval = (exp(tc*(t1-tmin)/2.) - exp(tc*(t2-tmin)/2.))**2
	      yran = rran()
	      if (yran.gt.yval) goto 21
c	invariant mass calculation: M with spectator, M1 with particle from pair
	      do il=1,4
		 pn_m(il) = pcm(il,1) + pcm(il,2)
	      end do
	      M_nnbar = scal_prod(pn_m,pn_m)
	      M_nnbar = pn_m(4)*pn_m(4) - M_nnbar
	      M_nnbar = sqrt(M_nnbar)
	      do il=1,4
		 pn_m(il) = pcm(il,1) + pcm(il,3)
	      end do
	      M_nn = scal_prod(pn_m,pn_m)
	      M_nn = pn_m(4)*pn_m(4) - M_nn
	      M_nn = sqrt(M_nn)
	   endif
	   if (jch.eq.40.or.jch.eq.41) then
c	t calculation for phase space in nucleon-antinucleon production
	      do il=1,3
		 pn_t(il) = pcm_t(il) - pcm(il,1)
	      end do
	      pn_t(4) = pcm_t(4) - pcm(4,1)
	      t1 = scal_prod(pn_t,pn_t)
	      t1 = pn_t(4)*pn_t(4)-t1
c	invariant mass calculation for phase space in nucleon-antinucleon production
	      do il=1,4
		 pn_m(il) = pcm(il,1) + pcm(il,2)
	      end do
	      M_nnbar = scal_prod(pn_m,pn_m)
	      M_nnbar = pn_m(4)*pn_m(4) - M_nnbar
	      M_nnbar = sqrt(M_nnbar)
	      do il=1,4
		 pn_m(il) = pcm(il,1) + pcm(il,3)
	      end do
	      M_nn = scal_prod(pn_m,pn_m)
	      M_nn = pn_m(4)*pn_m(4) - M_nn
	      M_nn = sqrt(M_nn)
	   endif
	   if (jch.eq.25.or.jch.eq.26) then
	      do il=1,4
		 pn_t(il) = pcm(il,2)
	      end do
	      tc = 3.
	      t = scal_prod(pn_t,pn_t)
	      t = pn_t(4)*pn_t(4)-t
	      yval = exp(tc*t)
	      yran = rran()
	      if (yran.gt.yval) goto 21
	   endif
	endif


c ***    writing process on control file, on bos file and on ntple file
	 do kk=1,np
	   call tree(parent(kk),pf,npd)
	   if(npd.eq.1) then
	      parent(kk)=pf(1)
	      call tree(parent(kk),pf,npd)
	   endif

	   if(npd.eq.0) then
c	   this pcle is stable.....
	     now=now+1
	     do j=1,5
		p(j,now)=pcmh(j,kk)
	     enddo
	     rm(now)=rmass(parent(kk))
c            we can use "rmass" since this pcle is stable
	     q(now)=charge(parent(kk))
	     npcode(now)=np_code(parent(kk))
	   else

c            this pcle decays....
	     later=later+1
c	     decaying pcle stack is updated
c            w, beta, amass, np, p0 ... storing ....
	     wd(later)=sqrt(pcm(4,kk)**2-pcm(5,kk)**2)
	     jch2(later)=chan(parent(kk))
	     do j=1,3
		bt(j,later)=pcmh(j,kk)/pcmh(4,kk)
		p0d(j,later)=pcmh(j,kk)
	     enddo
	     p0d(4,later)=pcmh(4,kk)
	     nd(later)=npd
	     do j=1,npd
	        rmd(j,later)=rmass(pf(j))
		partnew(j,later)=pf(j)
	     enddo
c	     pcle is temporarely stored...
	   endif
	 enddo
c	   if (jch.eq.0.and.new.eq.1) print *,'inside',npd

	 new=new+1
c        decay level stack is updated
	 if(new.le.later) then
c        there are still pending pcles...
	   np=nd(new)
	   w=wd(new)
	   jch=jch2(new)
	   do j=1,3
	      beta_X2cmh(j)=bt(j,new)
	      p0(j)=p0d(j,new)
	   enddo
	   p0(4)=p0d(4,new)
	   do j=1,np
	      amass(j)=rmd(j,new)
	      daughter(j)=partnew(j,new)
	   enddo
c          a deeper decay level starts..
c	     print *, new, wd(later)
c	     print *, 
c	   if (jch.eq.0.and.new.eq.1) print *,'leaving', tth,pph,pcm(5,2),pcm(4,2),sqrt(pcm(4,2)**2-pcm(5,2)**2)
	   
	   goto 20
	 else
c        all pcles have been scanned, the event is written on file(s)
c +++    adesso dobbiamo passare dal cmh al lab

	   do j=1,now
	     do il=1,4
	      p2(il)=p(il,j)
	     end do
c ***        boosting lungo beta rispetto al cmh
c ***        e rotazione...
	     beta2(1)=0.
	     beta2(2)=0.
	     beta2(3)=-vmod(beta_cmh2lab,3)
	     if(beta2(3).ne.0.) then
		call lorentz(beta2,p2)
	        call rotvec(p2,beta_cmh2lab)
	     endif
	     do il=1,4
	       p(il,j)=p2(il)
	     end do
	     p(5,j)=sqrt(scal_prod(p2,p2))
	   end do
c ***      controlling if energy and impulse are conserved (deltae = 1 Mev)
c ***      if they are not conserved, another extraction is performed
	   do j=1,4
	      delta(j)=0.
	   enddo
	   do j=1,4
	      do k=1,now
	         delta(j)=delta(j)+p(j,k)
	      enddo
	   enddo
	   do j=1,4
	      delta(j)=delta(j)-p0i(j)
	   enddo
	   delta(3)=delta(3)-eglab
	   delta(4)=delta(4)-eglab
	   do k=1,4
c	      print *, abs(delta(k))
	     if(abs(delta(k)).gt.0.005) goto 99
	   enddo
c+++ Checking if the event contain a proton:
c+++ If yes, checking if t<tmax
c+++ Only valid for the first generation of particle
	   if (new.eq.0) then
	      t_proton=1.
	   do k=1,now
	      if (npcode(k).eq.2212) then
		t_proton=-2.*0.93827*(p(4,k)-0.93827) 
		if (t_proton.gt.t_ref) goto 99
	      endif
	    enddo
	    endif
c---
	   if(ncheck.gt.0) then
       if (loop.le.ncheck)
     +     call writeout(jch,now,p,rm,q,npcode,egamma,p0i)
	   if (loop.eq.ncheck+1) close(unit=11)	 
	   if (loop.eq.ncheck+1) ncheck=0
	 endif

c-----------------------------------------------------------------------------------------
c++++ TEMPORARILY adding second nucleon (spectator) in case of 
c                 quasifree production on deuteron target
c
	 if(target.eq.'d'.or.target.eq.'deuterium'.AND.jch0.le.50) then
c+++
c	    p(1,now+1)=0.
c	    p(2,now+1)=0.
c	    p(3,now+1)=eglab
c	    p(4,now+1)=eglab+rmass('deuteron    ')
c	    do j=1,now
c	       p(1,now+1)=p(1,now+1)-p(1,j)
c	       p(2,now+1)=p(2,now+1)-p(2,j)
c	       p(3,now+1)=p(3,now+1)-p(3,j)
c	       p(4,now+1)=p(4,now+1)-p(4,j)
c	    enddo
c---
	    do j=1,4
c	     print*,j, p(j,now+1)-p0spec(j)
	       p(j,now+1)=p0spec(j)
	    enddo

	    p(5,now+1)=sqrt(p(1,now+1)**2+p(2,now+1)**2+p(3,now+1)**2)
	    if(partic(jch0,0).eq.'proton      ') then
	       rm(now+1)=rmass('neutron     ')
	       q(now+1) =0
	       npcode(now+1) =2112
	    else
	       rm(now+1)=rmass('proton      ')
	       q(now+1) =1
	       npcode(now+1) =2212
	    endif      
	    now=now+1
	 endif
c-----------------------------------------------------------------------------------------
	 numpar = now
	 if (bosout.eq.1) call fillbos(loop,fillbos_flag)
	 npart = now
	 if (hbookout.eq.1.AND.fillbos_flag) then
	    ch=jch0
	    eg=egamma
	    evt=true_evt
	    do j=1,npart
	       px(j)=p(1,j)
	       py(j)=p(2,j)
	       pz(j)=p(3,j)
	       ep(j)=p(4,j)
	       pp(j)=sqrt(px(j)**2+py(j)**2+pz(j)**2)
	       teta(j)=acos(pz(j)/pp(j))*180/pigr
	       if(px(j).gt.0.and.py(j).eq.0) then
		  phi(j)=+pigr/2
	       elseif(px(j).lt.0.and.py(j).eq.0) then
		  phi(j)=3*pigr/2
	       elseif(px(j).ge.0.and.py(j).gt.0) then
		  phi(j)=atan(py(j)/px(j))
	       elseif(px(j).ge.0.and.py(j).lt.0) then
		  phi(j)=2*pigr+atan(py(j)/px(j))
	       elseif(px(j).le.0.and.py(j).lt.0) then
		  phi(j)=atan(py(j)/px(j))+pigr
	       elseif(px(j).le.0.and.py(j).gt.0) then
		  phi(j)=atan(py(j)/px(j))+pigr
	       endif
	       phi(j)=phi(j)*180/pigr
	       bet(j)=pp(j)/ep(j)
	       mp(j)=rm(j)
	       qp(j)=q(j)
	       code(j)=npcode(j)
	    enddo
	    do j=1,np_nd
	       p_cm (j)=sqrt(pcm_nd(1,j)**2+pcm_nd(2,j)**2+pcm_nd(3,j)**2)
	       th_cm(j)=acos(pcm_nd(3,j)/p_cm(j))*180/pigr	       
	       m_nd (j)=sqrt(pcm_nd(4,j)**2-p_cm(j)**2)
	    enddo
	    np_dc=1
	    do j=1,np_dc
	       p_dc (j)=pcm_dc(5)
	       call angles(0,pcm_dc(1),pcm_dc(2),pcm_dc(3),th_dc(j),ph_dc(j))  
	       m_dc (j)=sqrt(pcm_dc(4)**2-pcm_dc(5)**2)
c	       print*, 'last',th_dc(j),ph_dc(j), m_dc (j)
c	       print*, 
	    enddo
	       t_mdst=t_proton
	    call hfnt(1) 
	 endif
       endif
c ***  end loop
c      }
	if(int(loop/2000.)*2000.eq.loop) print *,loop
       enddo
       return
       end

  							! *** end mc_clas ***







