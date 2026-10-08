  
	integer function np_code(p)
	
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c
c 									    c
c    this function assigns a particle code to each particle or nucleus      c
c    (useful for input to fastmc).				            c
c									    c 
c									    c 
c function value:							    c
c 									    c
c    np_code								    c
c 									    c
c---------------------------------------------------------------------------c

	implicit none
	character*12 p
	
	if(p.eq.'photon')     np_code=22
	if(p.eq.'nu_e')       np_code=0
	if(p.eq.'nu_mu')      np_code=0
	if(p.eq.'nu_tau')     np_code=0
	if(p.eq.'antinu_e')   np_code=0
	if(p.eq.'antinu_mu')  np_code=0
	if(p.eq.'antinu_tau') np_code=0
	if(p.eq.'electron')   np_code=11
	if(p.eq.'positron')   np_code=-11
	if(p.eq.'mu+')        np_code=-13
	if(p.eq.'mu-')        np_code=13
	if(p.eq.'proton')     np_code=2212
	if(p.eq.'antiproton') np_code=-2212
	if(p.eq.'neutron')    np_code=2112
	if(p.eq.'antineutron')    np_code=-2112
	if(p.eq.'lambda')     np_code=3122
	if(p.eq.'sigma+')     np_code=3222
	if(p.eq.'sigma0')     np_code=3212
	if(p.eq.'sigma-')     np_code=3112
	if(p.eq.'csi0')       np_code=3322
	if(p.eq.'csi-')       np_code=3312
	if(p.eq.'delta++')    np_code=2224
	if(p.eq.'delta+')     np_code=2214
	if(p.eq.'delta0')     np_code=2114
	if(p.eq.'delta-')     np_code=1114
	if(p.eq.'omega')      np_code=223
	if(p.eq.'pi+')        np_code=211
	if(p.eq.'pi0')        np_code=111
	if(p.eq.'pi-')        np_code=-211
	if(p.eq.'kappa+')     np_code=321
	if(p.eq.'kappa0')     np_code=311
	if(p.eq.'kappa-')     np_code=-321
	if(p.eq.'kappas')     np_code=310
	if(p.eq.'kappal')     np_code=130
	if(p.eq.'eta')        np_code=221
	if(p.eq.'etap')       np_code=331
	if(p.eq.'rho+')       np_code=213
	if(p.eq.'rho0')       np_code=113
	if(p.eq.'rho-')       np_code=-213
	if(p.eq.'phi')        np_code=333
	if(p.eq.'f2')         np_code=999
	if(p.eq.'deuteron')   np_code=0
	if(p.eq.'h3')         np_code=0
	if(p.eq.'he3')        np_code=0
	if(p.eq.'he4')        np_code=0
	if(p.eq.'dibaryon')   np_code=0

	return
	end



