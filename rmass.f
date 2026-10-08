c***************************************************************************c

	real function rmass(p)
	
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c
c 									    c
c    this function assigns the mass to each particle or nucleus             c
c    (mass in gev).							    c 
c 									    c
c 									    c
c function value:							    c
c 									    c
c    rmass								    c
c 									    c
c---------------------------------------------------------------------------c

	implicit none

	real res_mass,fun_bg,kstar_mass,hyper_mass,dipion_mass
	character*12 p
	rmass = -1.

	if(p.eq.'photon')        then
	   rmass = 0.0
	elseif(p.eq.'neutrino') then
	   rmass = 0.0
	elseif(p.eq.'none')     then
	   rmass = 0.0
	elseif(p.eq.'electron') then
	   rmass = 0.511e-03
	elseif(p.eq.'positron') then
	   rmass = 0.511e-03
	elseif(p.eq.'mu+')      then
	   rmass = 0.10566
	elseif(p.eq.'mu-')      then
	   rmass = 0.10566
	elseif(p.eq.'proton')   then
	   rmass = 0.93827
	elseif(p.eq.'antiproton')   then
	   rmass = 0.93827
	elseif(p.eq.'neutron')  then
	   rmass = 0.93956
	elseif(p.eq.'antineutron')  then
	   rmass = 0.93956
	elseif(p.eq.'lambda')   then
	   rmass = 1.116
	elseif(p.eq.'sigma+')   then
	   rmass = 1.189
	elseif(p.eq.'sigma0')   then
	   rmass = 1.192
	elseif(p.eq.'sigma-')   then
	   rmass = 1.197
	elseif(p.eq.'csi0')     then
	   rmass = 1.3149
	elseif(p.eq.'csi-')     then
	   rmass = 1.3213
	elseif(p.eq.'delta++')  then
	   rmass = res_mass(1.235, 0.115, 1.077)
	elseif(p.eq.'delta+')   then
	   rmass = res_mass(1.235, 0.115, 1.077)
	elseif(p.eq.'delta0')   then
	   rmass = res_mass(1.235, 0.115, 1.077)
	elseif(p.eq.'delta-')   then
	   rmass = res_mass(1.235, 0.115, 1.077)
	elseif(p.eq.'sigma+3/2')then
	   rmass = 1.382
	elseif(p.eq.'sigma03/2')then
	   rmass = 1.387
	elseif(p.eq.'sigma-3/2')then
	   rmass = 1.382
	elseif(p.eq.'csi03/2')  then
	   rmass = 1.532
	elseif(p.eq.'csi-3/2')  then
	   rmass = 1.535
	elseif(p.eq.'omega-')   then
	   rmass = 1.67245
	elseif(p.eq.'omega')    then
	   rmass = res_mass(0.7826, 0.00843, 0.41)
        elseif(p.eq.'omega_mod')then
           rmass = fun_bg(1.)
	elseif(p.eq.'pi+')      then
	   rmass = 0.13957
	elseif(p.eq.'pi0')      then
	   rmass = 0.13496
	elseif(p.eq.'pi-')      then
	   rmass = 0.13957
	elseif(p.eq.'kappa+')   then
	   rmass = 0.49367
	elseif(p.eq.'kappa0')   then
	   rmass = 0.49772
	elseif(p.eq.'kappa-')   then
	   rmass = 0.49367
	elseif(p.eq.'kappas')   then
	   rmass = 0.49772
	elseif(p.eq.'kappal')   then
	   rmass = 0.49772
	elseif(p.eq.'akappa0')  then
	   rmass = 0.49772
	elseif(p.eq.'eta')      then
	   rmass = 0.5488
	elseif(p.eq.'etap')     then
	   rmass = 0.95747
	elseif(p.eq.'rho+')     then
	   rmass = res_mass(0.765, 0.149, 0.275)
	elseif(p.eq.'rho0')     then
	   rmass = res_mass(0.770, 0.149, 0.275)
	elseif(p.eq.'rho-')     then
	   rmass = res_mass(0.765, 0.149, 0.275)
	elseif(p.eq.'f0')       then
	   rmass = res_mass(0.980, 0.040, 0.860)
	elseif(p.eq.'phi')      then
	   rmass = 1.020
	elseif(p.eq.'f2')       then
	   rmass = res_mass(1.275, 0.185, 0.550)
	elseif(p.eq.'a2')       then
	   rmass = res_mass(1.319, 0.105, 1.005)
	elseif(p.eq.'rho3')       then
	   rmass = res_mass(1.690, 0.200, 1.100)
	elseif(p.eq.'deuteron') then
	   rmass = 1.8755
	elseif(p.eq.'diproton') then
	   rmass = 1.8755
	elseif(p.eq.'h3')       then
	   rmass = 2.80875
	elseif(p.eq.'he3')      then
	   rmass = 2.80822
	elseif(p.eq.'he4')      then
	   rmass = 3.72715
	elseif(p.eq.'dibaryon') then
	   rmass = res_mass(2.23 , 0.017 ,1.88)	   
	elseif(p.eq.'pp')       then
	   rmass = res_mass(2.23 , 0.115 ,1.88)
	elseif(p.eq.'ppbar')    then
	   rmass = res_mass(2.23 , 0.115 ,1.88)
	elseif(p.eq.'lambda*')  then
	   rmass = res_mass(1.5195, 0.016, 1.47)
	elseif(p.eq.'zeta+')    then
	   rmass = 1.530
	elseif(p.eq.'zeta++')    then
	   rmass = 1.580
	elseif(p.eq.'kappa*')   then
	   rmass = kstar_mass(0.892, 0.050,0.822, 0.55,0.35,0.620)
	elseif(p.eq.'kappa0*')   then
	   rmass = res_mass(0.892 , 0.050 ,0.742)
	elseif(p.eq.'hyper0')   then
	   rmass = hyper_mass(2.08,-0.18,0.2,1.600,2.900)
	elseif(p.eq.'hyper++')  then
	   rmass = hyper_mass(2.08,-0.3,0.3,1.600,2.900)
	elseif(p.eq.'s11')  then
	   rmass = res_mass(1.54, 0.3, 1.5)
	elseif(p.eq.'dipion')  then
	   rmass = dipion_mass(0.4, 1.4)
	endif
	return
	end

						       ! *** end rmass ***

	real function res_mass(meanval,gamma,thresh)

	implicit none 

	real meanval,gamma,thresh
	real q,rran,pigr
	data pigr/3.14159265/
 1	continue
	q=rran()
        res_mass = (meanval + gamma/2.*tan(pigr*(q-0.5)))
	if(res_mass.lt.thresh) goto 1

	return
	end
				                     ! *** end res_rmass ***

        real function fun_bg(a)
	implicit none
c       It makes sense for 0.4<M<1.1

	real q,rran,y,fun,a
	
 1	continue
	q=rran()*1.1
	if (q.lt.0.4) goto 1
	y=rran()*3500. ! 2068 is just normalization
        fun = -1316.+9514.*q-21080.*q**2+14950.*q**3
	if(y.gt.fun) goto 1
	fun_bg=q
	

        return
	end



	real function kstar_mass(meanval1,gamma1,meanval2,gamma2,sig_noise,thresh)

	implicit none 

	real meanval1,meanval2,gamma1,gamma2,sig_noise,thresh
	real q,y,rran,pigr
	data pigr/3.14159265/
	y=rran()
 1	continue
	if(y.lt.sig_noise) then
	   q=rran()
	   kstar_mass = (meanval1 + gamma1/2.*tan(pigr*(q-0.5)))
	else
	   q=rran()
	   kstar_mass = (meanval2 + gamma2/2.*tan(pigr*(q-0.5)))
	endif
	if(kstar_mass.lt.thresh) goto 1

	return
	end


	real function hyper_mass(meanval,gamma1,gamma2,thr_l,thr_h)

	implicit none 

	real meanval,gamma1,gamma2,thr_l,thr_h
	real x,y,rran,pigr
	data pigr/3.14159265/

 1	continue
	x=rran()*(thr_h-thr_l)+thr_l
	y=exp(-(x-meanval)**2/2/(gamma1+x*gamma2)**2)
	if(rran().gt.y) goto 1
	hyper_mass=x

	return
	end




	real function dipion_mass(xmin,xmax)

	implicit none
	real xmin,xmax,rran

	dipion_mass=rran()*(xmax-xmin)+xmin
	
	return
	end
