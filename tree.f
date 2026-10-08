
  
	subroutine tree(p,pf,np)
	
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c
c 									    c
c    this function assigns a number of possible dacay modes                 c 
c    to each particle                                                       c
c    one can add (or delete) any particle in the following list:            c
c    the net result will be that this particle will be considered           c
c    unstable (stable)                                                      c
c 									    c
c    p  = particle name                                                     c
c    pf = names of decaying products                                        c
c    np = number of decaying products (0 if particle p is stable)           c
c---------------------------------------------------------------------------c

	implicit none
	character*12 p,pf(10)
	integer np,j,nbr
	real br(10)
	real q,rran

	np=0
	do j=1,10
	 pf(j)='none'
	 br(j)=0.
	enddo

c	if(p.eq.'pi+') then
c	   np=2
c	   pf(1)='mu+'
c	   pf(2)='neutrino'
c	   return
c	endif

 	if(p.eq.'pi0') then
	   np=2
	   pf(1)='photon'
	   pf(2)='photon'
	   return
	endif

c	if(p.eq.'pi-') then
c	   np=2
c	   pf(1)='mu-'
c	   pf(2)='neutrino'
c	   return
c	endif

	if(p.eq.'lambda')     then
	   nbr=2
	   br(1)=0.64
	   br(2)=0.36
	   q=rran()
	   if(q.lt.br(1)) then
	      np=2
	      pf(1)='pi-'
	      pf(2)='proton'
	   else
	      np=2
	      pf(1)='pi0'
	      pf(2)='neutron'
	   endif
	   return
	endif

	if(p.eq.'delta++')    then
	   np=2
	   pf(1)='pi+'
	   pf(2)='proton'
	   return
	endif

	if(p.eq.'delta+')     then
	   nbr=2
	   br(1)=0.33
	   br(2)=0.67
	   q=rran()
	   if(q.lt.br(1)) then
	      np=2
	      pf(1)='pi+'
	      pf(2)='neutron'
	   else
	      np=2
	      pf(1)='pi0'
	      pf(2)='proton'
	   endif	 
	   return
	endif

	if(p.eq.'delta0')     then
	   nbr=2
	   br(1)=0.33
	   br(2)=0.67
	   q=rran()
	   if(q.lt.br(1)) then
	      np=2
	      pf(1)='pi-'
	      pf(2)='proton'
	   else
	      np=2
	      pf(1)='pi0'
	      pf(2)='neutron'
	   endif
	   return
	endif

	if(p.eq.'delta-')     then
	   nbr=1
	   np=2
	   pf(1)='pi-'
	   pf(2)='neutron'
	   return
	endif

	if(p.eq.'omega')     then
	   nbr=1
	   np=3
	   pf(1)='pi+'
	   pf(2)='pi-'
	   pf(3)='pi0'
	   return
	endif

	if(p.eq.'eta')        then
	   nbr=4
	   br(1)=0.39
	   br(2)=0.32
	   br(3)=0.24
	   br(4)=0.05
	   q=rran()
	   if(q.lt.br(1)) then
	      np=2
	      pf(1)='photon'
	      pf(2)='photon'
	   else if(q.lt.br(1)+br(2))  then
	      np=3
	      pf(1)='pi0'
	      pf(2)='pi0'
	      pf(3)='pi0'
	   else if(q.lt.br(1)+br(2)+br(3))  then
	      np=3
	      pf(1)='pi+'
	      pf(2)='pi0'
	      pf(3)='pi-'
	   else if(q.lt.br(1)+br(2)+br(3)+br(4))  then
	      np=3
	      pf(1)='pi+'
	      pf(2)='pi-'
	      pf(3)='photon'
	   endif
	   return
	endif
   
	if(p.eq.'rho+')     then
	   nbr=1
	   np=2
	   pf(1)='pi+'
	   pf(2)='pi0'
	   return
	endif

	if(p.eq.'rho0')     then
	   nbr=1
	   np=2
	   pf(1)='pi+'
	   pf(2)='pi-'
	   return
	endif

	if(p.eq.'rho-')     then
	   nbr=1
	   np=2
	   pf(1)='pi0'
	   pf(2)='pi-'
	   return
	endif

	if(p.eq.'phi')      then
	   nbr=5
	   br(1)=0.50
	   br(2)=0.35
	   br(3)=0.05
	   br(4)=0.05
	   br(5)=0.05
	   q=rran()
	   if(q.lt.br(1)) then
	      np=2
	      pf(1)='kappa+'
	      pf(2)='kappa-'
	   else if(q.lt.br(1)+br(2))  then
	      np=2
	      pf(1)='kappas'
	      pf(2)='kappal'
	   else if(q.lt.br(1)+br(2)+br(3))  then
	      np=2
	      pf(1)='pi-'
	      pf(2)='rho+'
	   else if(q.lt.br(1)+br(2)+br(3)+br(4))  then
	      np=2
	      pf(1)='pi0'
	      pf(2)='rho0'
	   else if(q.lt.br(1)+br(2)+br(3)+br(4)+br(5))  then
	      pf(1)='pi+'
	      pf(2)='rho-'
	   endif
	   return
	endif

c	if(p.eq.'phi')      then
c	   nbr=2
c	   br(1)=0.6
c	   br(2)=0.4
c	   q=rran()
c	   if(q.lt.br(1)) then
c	      np=2
c	      pf(1)='kappa+'
c	      pf(2)='kappa-'
c	   else
c	      np=2
c	      pf(1)='kappas'
c	      pf(2)='kappal'
c	   endif
c	   return
c	endif

	if(p.eq.'f2')      then
	   nbr=2
	   br(1)=0.9
	   br(2)=0.1
	   q=rran()
	   if(q.lt.br(1)) then
	      np=2
	      pf(1)='pi+'
	      pf(2)='pi-'
	   else
	      np=4
	      pf(1)='pi+'
	      pf(2)='pi-'
	      pf(3)='pi0'
	      pf(4)='pi0'
	   endif
	   return
	endif

	if(p.eq.'omega_mod')  then
	   np=3
	   pf(1)='pi+'
	   pf(2)='pi-'
	   pf(3)='pi0'
	   return
	endif

	if(p.eq.'dibaryon')  then
	   nbr=1
	   np=2
	   pf(1)='proton'
	   pf(2)='proton'
	   return
	endif


	if(p.eq.'pp')  then
	   nbr=1
	   np=2
	   pf(1)='proton'
	   pf(2)='proton'
	   return
	endif

	if(p.eq.'ppbar')  then
	   nbr=1
	   np=2
	   pf(1)='proton'
	   pf(2)='antiproton'
	   return
	endif



	return
	end










