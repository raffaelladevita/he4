c***************************************************************************c

	subroutine sig_data
  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c 
c 									    c
c    this subroutine defines the total cross sections for the different     c
c    interaction channels       					    c
c---------------------------------------------------------------------------c
c cross setions are given as a function of W.
c
c for reactions photon-nucleon
c  1.0 GeV < W < 3.5 GeV (dW=0.01)
c
c for reactions photon-deuteron
c  1.9 GeV < W < 5.1 GeV  
c
c for reactions photon-helium3
c  2.9 GeV < W < 6.4 GeV 
c
c all corresponding to   -->  0.063 GeV < Eg_lab < 6 GeV
c--------------------------------------------------------------------------c
c  PHOTON INTERACTION OUTPUT CHANNELS
c 1)  gamma  +  p  --->  pi+  +  n
c 2)  gamma  +  p  --->  pi0  +  p
c 3)  gamma  +  n  --->  pi-  +  p
c 4)  gamma  +  n  --->  pi0  +  n
c 5)  gamma  +  p  --->  delta++ +  pi-
c 6)  gamma  +  p  --->  delta+  +  pi0
c 7)  gamma  +  p  --->  delta0  +  pi+
c 8)  gamma  +  n  --->  delta+  +  pi-
c 9)  gamma  +  n  --->  delta0  +  pi0
c 10) gamma  +  n  --->  delta-  +  pi+
c 11) gamma  +  p  --->  rho0  +  p
c 12) gamma  +  p  --->  rho+  +  n
c 13) gamma  +  n  --->  rho-  +  p
c 14) gamma  +  n  --->  rho0  +  n
c 15) gamma  +  p  --->  pi+  +  pi-  +  p                           
c 16) gamma  +  p  --->  pi0  +  pi+  +  n                           
c 17) gamma  +  n  --->  pi+  +  pi-  +  n                           
c 18) gamma  +  n  --->  pi0  +  pi-  +  p                           
c 19) gamma  +  p  --->  omega  +  p                                 
c 20) gamma  +  n  --->  omega  +  n                                 
c 21) gamma  +  p  --->  pi+  +  pi-  +  pi0  +  p                   
c 22) gamma  +  p  --->  pi+  +  pi+  +  pi-  +  n                   
c 23) gamma  +  n  --->  pi+  +  pi-  +  pi0  +  n                   
c 24) gamma  +  n  --->  pi+  +  pi-  +  pi-  +  p                   
c 25) gamma  +  p  --->  phi  +  p                                   
c 26) gamma  +  n  --->  phi  +  n                                   
c 27) gamma  +  p  --->  pi+  +  pi+  +  pi-  +  pi-  +  p                   
c 28) gamma  +  p  --->  pi+  +  pi+  +  pi-  +  pi-  +  pi0  +  p             
c 29) gamma  +  p  --->  pi+  +  pi+  +  pi+  +  pi-  +  pi-  +  n 
c 30) gamma  +  p  --->  lambda  +  kappa+                
c 31) gamma  +  p  --->  kappa+  +  kappa-  +  p
c 32) gamma  +  p  --->  eta  +  p
c 33) gamma  +  p  --->  f2  +  p
c 34) gamma  +  n  --->  f2  +  n
c 
c 51) gamma  +  3He  --->  p + p + n
c 52) gamma  +  3He  --->  p + p + n + pi0
c 53) gamma  +  3He  --->  p + p + p + pi-  
c 54) gamma  +  3He  --->  p + p + n + pi+  
c 55) gamma  +  3He  --->  3He + pi0  
c 56) gamma  +  3He  --->  3H  + pi+ 
c 57) gamma  +  3He  --->  d  + p 
c 58) gamma  +  3He  --->  d  + p + pi0 
c 59) gamma  +  3He  --->  d  + n + pi+ 
c 60) gamma  +  d    --->  p  + n 
c 61) gamma  +  d    --->  p  + n + pi0
c 62) gamma  +  d    --->  p  + p + pi-
c 63) gamma  +  d    --->  n  + p + pi+
c 70) gamma  +  pp   --->  p  + p
c 71) gamma  +  pp   --->  p  + p + pi0
c 72) gamma  +  pp   --->  p  + n + pi+
c 80) gamma  +  d    --->  dibaryon + pi-
c-------------------------------------------------------------------------
c RESONANCE DECAY CHANNELS
c 101) delta++ decay 
c 102) delta+  decay 
c 103) delta0  decay 
c 104) delta-  decay 
c 105) rho+    decay 
c 106) rho0    decay 
c 107) rho-    decay 
c 109) phi     decay 
c 110) lambda  decay
c 112) f2      decay 
c--------------------------------------------------------------------------
	implicit none

	integer nchan
	parameter ( nchan = 100 )
        character*12 target,partic(nchan,0:8)
        common /target/target
	integer numpart(nchan)
	common/deca/numpart,partic
	integer jch,j,k
	real wr(251),egl(251)
	real sigr(nchan,251)
	common/sigtot/sigr,wr,egl
	real mn
	data mn/0.9385/
	real dw
	data dw/0.01/

c  51) gamma  +  3He  --->  p + p + n
	data (sigr(51,j),j=1,251) /251*10./
c  52) gamma  +  3He  --->  p + p + n + pi0
	data (sigr(52,j),j=1,251) /251*10./
c  53) gamma  +  3He  --->  p + p + p + pi-  
	data (sigr(53,j),j=1,251) /251*10./
c  54) gamma  +  3He  --->  p + p + n + pi+  
	data (sigr(54,j),j=1,251) /251*10./
c  55) gamma  +  3He  --->  3He + pi0  
	data (sigr(55,j),j=1,251) /251*3./
c  56) gamma  +  3He  --->  3H  + pi+ 
	data (sigr(56,j),j=1,251) /251*3./
c  57) gamma  +  3He  --->  d  + p 
	data (sigr(57,j),j=1,251) /251*3./
c  58) gamma  +  3He  --->  d  + p + pi0 
	data (sigr(58,j),j=1,251) /251*3./
c  59) gamma  +  3He  --->  d  + n + pi+ 
	data (sigr(59,j),j=1,251) /251*3./
c  60) gamma  +  d    --->  p  + n 
	data (sigr(60,j),j=1,251) /251*10./
c  61) gamma  +  d    --->  p  + n + pi0
	data (sigr(61,j),j=1,251) /251*10./
c  62) gamma  +  d    --->  p  + p + pi-
	data (sigr(62,j),j=1,251) /251*10./
c  63) gamma  +  d    --->  n  + p + pi+
	data (sigr(63,j),j=1,251) /251*10./
c  70) gamma  +  pp   --->  p  + p
	data (sigr(70,j),j=1,251) /251*10./
c  71) gamma  +  pp   --->  p  + p + pi0
	data (sigr(71,j),j=1,251) /251*10./
c  72) gamma  +  pp   --->  p  + n + pi+
	data (sigr(72,j),j=1,251) /251*10./
c  80) gamma  +  d    --->  dibaryon + pi-
	data (sigr(80,j),j=1,251) /251*1./

	wr(1)=1.00
	do j=2,251
	   wr(j)=wr(j-1)+dw
	enddo
	do j=1,251
	   egl(j)=(wr(j)**2-mn**2)/(2*mn)
	enddo

c       redefing some channels....
	do k=1,251
	   sigr(7,k)  = sigr(5,k)/3.
	   sigr(10,k) = sigr(8,k)*2.
	   sigr(13,k) = sigr(12,k)
	   sigr(14,k) = sigr(12,k)*2./3.
	   sigr(16,k) = sigr(15,k)/2.
	   sigr(17,k) = sigr(15,k)
	   sigr(18,k) = sigr(15,k)/2.
	   sigr(20,k) = sigr(19,k)
	   sigr(23,k) = sigr(21,k)
	   sigr(24,k) = sigr(22,k)*5./3.
	   sigr(26,k) = sigr(25,k)*5/3.
	   sigr(33,k) = sigr(11,k)/16.
	   sigr(34,k) = sigr(14,k)/16.
	enddo

c       redefining nucleus cross section if target = 'he3' or 'h3'
c       because number of protons and neutrons is different....
	if(target.eq.'he3') then
	 do jch=1,nchan
	  if(partic(jch,0).eq.'proton') then
	   do j=1,251
	    sigr(jch,j)=2.*sigr(jch,j)
	   enddo
	  endif
	 enddo
	endif
	if(target.eq.'h3'.or.target.eq.'tritium') then
	 do jch=1,nchan
	  if(partic(jch,0).eq.'neutron') then
	   do j=1,251
	    sigr(jch,j)=2.*sigr(jch,j)
	   enddo
	  endif
	 enddo
	endif

	return		
	end
 
							! *** end sig_data ***





