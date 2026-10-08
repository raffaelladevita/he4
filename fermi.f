c***************************************************************************c

	subroutine fermi(jch,w,beta,p0)
  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c
c 									    c
c    this subroutine calculates the system invariant mass accounting        c
c    for fermi motion of nucleon into nucleus. it also defines the	    c
c    four_momentum of the three nucleons and c.m. beta.			    c
c									    c
c    momentum distribution of nucleon used is ciofi calculation             c
c    (perspectives in nuclear physics at intermediate energies ictp         c
c     18 - 22 may 1987 pg. 367).					    c
c									    c
c       <sig_p> =  0.040 gev/c in he3					    c
c---------------------------------------------------------------------------c

	implicit none

	integer nchan
	parameter (nchan = 100 )
        character*12 partic(nchan,0:8)
        character*12 target
        real eglab,fwhm
	real rran
	real pi
	real rmass
	integer numpart(nchan)
        common /beam1/eglab,fwhm
        common /target/target
	common/deca/numpart,partic

 	real pn(4,4),beta(3),rm(4),p0(4)
	real sig_p,q,a,b,etot,s,w
	integer jch,j

	data pi/6.2831596/

	rm(1) = rmass(partic(jch,0))

c ***   target is a nucleon at rest
c ***   or target is 3He/deuteron/diproton
	if(target.eq.'proton'.or.
     +     target.eq.'neutron'.or.
     +     (partic(jch,0).ne.'proton'.and.partic(jch,0).ne.'neutron'))
     +  then
	   p0(1)=0.
	   p0(2)=0.
	   p0(3)=0.
	   p0(4)=rm(1)
	endif

c ***   target is a proton/neutron inside 4He
	if(target.eq.'he4') then
	   sig_p=0.030
	   if(partic(jch,0).eq.'proton') then
	      rm(2)=rmass('proton      ')
	      rm(3)=rmass('neutron     ')
	      rm(4)=rmass('neutron     ')
	   else
	      rm(2)=rmass('proton      ')
	      rm(3)=rmass('neutron     ')
	      rm(3)=rmass('proton      ')
	   endif  
	   do j=2,4
	      a = acos(1.-2.*rran())
	      b = 2.*pi*rran()
	      call norran(q)
	      pn(1,j) = sig_p*abs(q)*sin(a)*cos(b)
	      pn(2,j) = sig_p*abs(q)*sin(a)*sin(b)
	      pn(3,j) = sig_p*abs(q)*cos(a)
	      pn(4,j) = sig_p*abs(q)*cos(a)
	   end do
	   do j=2,4
	      pn(4,j)=sqrt(pn(1,j)**2+pn(2,j)**2+pn(3,j)**2+rm(j)**2)
           end do
c ***      second and third nucleon have the above four_momentum
c ***      first nucleon interacts with photon. its momentum is:
	   pn(1,1)=-(pn(1,2)+pn(1,3)+pn(1,4))
	   pn(2,1)=-(pn(2,2)+pn(2,3)+pn(2,4))
	   pn(3,1)=-(pn(3,2)+pn(3,3)+pn(3,4))
	   pn(4,1)=sqrt(pn(1,1)**2+pn(2,1)**2+pn(3,1)**2+rm(1)**2)
	   p0(1)=pn(1,1)
	   p0(2)=pn(2,1)
	   p0(3)=pn(3,1)
	   p0(4)=pn(4,1)
	endif


c ***   target is a proton/neutron inside 3He
	if(target.eq.'he3') then
	   sig_p=0.040
	   if(partic(jch,0).eq.'proton') then
	      rm(2)=rmass('proton      ')
	      rm(3)=rmass('neutron     ')
	   else
	      rm(2)=rmass('proton      ')
	      rm(3)=rmass('proton      ')
	   endif  
	   do j=2,3
	      a = acos(1.-2.*rran())
	      b = 2.*pi*rran()
	      call norran(q)
	      pn(1,j) = sig_p*abs(q)*sin(a)*cos(b)
	      pn(2,j) = sig_p*abs(q)*sin(a)*sin(b)
	      pn(3,j) = sig_p*abs(q)*cos(a)
	   end do
	   do j=2,3
	      pn(4,j)=sqrt(pn(1,j)**2+pn(2,j)**2+pn(3,j)**2+rm(j)**2)
           end do
c ***      second and third nucleon have the above four_momentum
c ***      first nucleon interacts with photon. its momentum is:
	   pn(1,1)=-(pn(1,2)+pn(1,3))
	   pn(2,1)=-(pn(2,2)+pn(2,3))
	   pn(3,1)=-(pn(3,2)+pn(3,3))
	   pn(4,1)=sqrt(pn(1,1)**2+pn(2,1)**2+pn(3,1)**2+rm(1)**2)
	   p0(1)=pn(1,1)
	   p0(2)=pn(2,1)
	   p0(3)=pn(3,1)
	   p0(4)=pn(4,1)
	endif

c ***   target is a proton/neutron inside 3H
	if(target.eq.'h3'.or.target.eq.'tritium') then
	   sig_p=0.040
c	   sig_p=0.
	   if(partic(jch,0).eq.'proton') then
	      rm(2)=rmass('neutron     ')
	      rm(3)=rmass('neutron     ')
	   else
	      rm(2)=rmass('proton      ')
	      rm(3)=rmass('neutron     ')
	   endif   
	   do j=2,3
	      a = acos(1.-2.*rran())
	      b = 2.*pi*rran()
	      call norran(q)
	      pn(1,j) = sig_p*abs(q)*sin(a)*cos(b)
	      pn(2,j) = sig_p*abs(q)*sin(a)*sin(b)
	      pn(3,j) = sig_p*abs(q)*cos(a)
	   end do
	   do j=2,3
	      pn(4,j)=sqrt(pn(1,j)**2+pn(2,j)**2+pn(3,j)**2+rm(j)**2)
           end do 
c ***      second and third nucleon have this four_momentum
c ***      first nucleon interacts with photon. its momentum is:
	   pn(1,1)=-(pn(1,2)+pn(1,3))
	   pn(2,1)=-(pn(2,2)+pn(2,3))
	   pn(3,1)=-(pn(3,2)+pn(3,3))
	   pn(4,1)=sqrt(pn(1,1)**2+pn(2,1)**2+pn(3,1)**2+rm(1)**2)
	   p0(1)=pn(1,1)
	   p0(2)=pn(2,1)
	   p0(3)=pn(3,1)
	   p0(4)=pn(4,1)
	endif

c ***   target is e nucleon inside a deuteron
	if(target.eq.'d'.or.target.eq.'deuterium') then
	   sig_p=0.040
c	   sig_p=0.
	   if(partic(jch,0).eq.'proton') then
	      rm(2)=rmass('neutron     ')
	   else
	      rm(2)=rmass('proton      ')
	   endif   
	   a = acos(1.-2.*rran())
	   b = 2.*pi*rran()
	   call norran(q)
	   pn(1,2) = sig_p*abs(q)*sin(a)*cos(b)
	   pn(2,2) = sig_p*abs(q)*sin(a)*sin(b)
	   pn(3,2) = sig_p*abs(q)*cos(a)
	   pn(4,2)=sqrt(pn(1,2)**2+pn(2,2)**2+pn(3,2)**2+rm(2)**2)
c ***      second nucleon has this four_momentum
c ***      first nucleon interacts with photon. its momentum is:
	   pn(1,1)=-pn(1,2)
	   pn(2,1)=-pn(2,2)
	   pn(3,1)=-pn(3,2)
	   pn(4,1)=sqrt(pn(1,1)**2+pn(2,1)**2+pn(3,1)**2+rm(1)**2)
	   p0(1)=pn(1,1)
	   p0(2)=pn(2,1)
	   p0(3)=pn(3,1)
	   p0(4)=pn(4,1)
	endif

c ***   calculating w
c ***   s is lorentz invariant (same in c.m. and lab. system)
c ***   scalar product has only third component, assuming photon direction
c ***   corresponding with z axis

	s = rm(1)**2+2.*p0(4)*eglab-2.*p0(3)*eglab
	w = sqrt(s)

c ***   calculating beta
c ***   beta = beta of hadronic cm (cmh) with respect to lab
	etot = eglab+p0(4)
	beta(1)= p0(1)/etot
	beta(2)= p0(2)/etot
	beta(3)=(p0(3)+eglab)/etot

	return
	end
  
							  ! *** end fermi ***







