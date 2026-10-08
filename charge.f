c***************************************************************************c

	real function charge(p)
	
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c
c 									    c
c    this function assigns a charge to each particle or nucleus	  	    c
c    (charge in e.u.).							    c
c 									    c
c function value:							    c
c 									    c
c    charge								    c
c 									    c
c---------------------------------------------------------------------------c

	implicit none
	character*12 p

	charge = 0.
	if(p.eq.'photon')     charge =  0.0
	if(p.eq.'neutrino')   charge =  0.0
	if(p.eq.'none')       charge =  0.0
	if(p.eq.'electron')   charge = -1.0
	if(p.eq.'positron')   charge =  1.0
	if(p.eq.'mu+')        charge =  1.0
	if(p.eq.'mu-')        charge = -1.0
	if(p.eq.'proton')     charge =  1.0
	if(p.eq.'antiproton') charge =  -1.0
	if(p.eq.'neutron')    charge =  0.0
	if(p.eq.'antineutron') charge =  0.0
	if(p.eq.'lambda')     charge =  0.0
	if(p.eq.'sigma+')     charge =  1.0
	if(p.eq.'sigma0')     charge =  0.0
	if(p.eq.'sigma-')     charge = -1.0
	if(p.eq.'csi0')       charge =  0.0
	if(p.eq.'csi-')       charge = -1.0
	if(p.eq.'delta')      charge =  0.0
	if(p.eq.'delta++')    charge =  2.0
	if(p.eq.'delta+')     charge =  1.0
	if(p.eq.'delta0')     charge =  0.0
	if(p.eq.'delta-')     charge = -1.0
	if(p.eq.'sigma+3/2')  charge =  1.0
	if(p.eq.'sigma03/2')  charge =  0.0
	if(p.eq.'sigma-3/2')  charge = -1.0
	if(p.eq.'csi03/2')    charge =  0.0
	if(p.eq.'csi-3/2')    charge = -1.0
	if(p.eq.'omega-')     charge = -1.0
	if(p.eq.'omega')      charge =  0.0
	if(p.eq.'omega_mod')  charge =  0.0
	if(p.eq.'pi+')        charge =  1.0
	if(p.eq.'pi0')        charge =  0.0
	if(p.eq.'pi-')        charge = -1.0
	if(p.eq.'kappa+')     charge =  1.0
	if(p.eq.'kappa0')     charge =  0.0
	if(p.eq.'kappa-')     charge = -1.0
	if(p.eq.'kappas')     charge =  0.0
	if(p.eq.'kappal')     charge =  0.0
	if(p.eq.'akappa0')    charge =  0.0
	if(p.eq.'eta')        charge =  0.0
	if(p.eq.'etap')       charge =  0.0
	if(p.eq.'rho+')       charge =  1.0
	if(p.eq.'rho0')       charge =  0.0
	if(p.eq.'rho-')       charge = -1.0
	if(p.eq.'phi')        charge =  0.0
	if(p.eq.'f2')         charge =  0.0
	if(p.eq.'deuteron')   charge =  1.0
	if(p.eq.'h3')         charge =  1.0
	if(p.eq.'he3')        charge =  2.0
	if(p.eq.'he4')        charge =  2.0
	if(p.eq.'dibaryon')   charge =  2.0
	if(p.eq.'pp') 	      charge =  2.0
	if(p.eq.'ppbar')      charge =  0.0
	return
	end

							 ! *** end charge ***  
