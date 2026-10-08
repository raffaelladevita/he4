c***************************************************************************c
 
	subroutine decay
  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c
c 									    c
c    this subroutine defines final states of interactions.                  c
c    each particle is characterized by some numbers through variable        c
c    									    c
c---------------------------------------------------------------------------c
ccc     PHOTON INTERACTION OUTPUT CHANNELS

c  1)  gamma  +  p  --->  pi+  +  n
c  2)  gamma  +  p  --->  pi0  +  p
c  3)  gamma  +  n  --->  pi-  +  p
c  4)  gamma  +  n  --->  pi0  +  n
c  5)  gamma  +  p  --->  delta++ +  pi-
c  6)  gamma  +  p  --->  delta+  +  pi0
c  7)  gamma  +  p  --->  delta0  +  pi+
c  8)  gamma  +  n  --->  delta+  +  pi-
c  9)  gamma  +  n  --->  delta0  +  pi0
c  10) gamma  +  n  --->  delta-  +  pi+
c  11) gamma  +  p  --->  rho0  +  p
c  12) gamma  +  p  --->  rho+  +  n
c  13) gamma  +  n  --->  rho-  +  p
c  14) gamma  +  n  --->  rho0  +  n
c  15) gamma  +  p  --->  pi+  +  pi-  +  p                           
c  16) gamma  +  p  --->  pi0  +  pi+  +  n                           
c  17) gamma  +  n  --->  pi+  +  pi-  +  n                           
c  18) gamma  +  n  --->  pi0  +  pi-  +  p                           
c  19) gamma  +  p  --->  omega  +  p                                 
c  20) gamma  +  n  --->  omega  +  n                                 
c  21) gamma  +  p  --->  pi+  +  pi-  +  pi0  +  p                   
c  22) gamma  +  p  --->  pi+  +  pi+  +  pi-  +  n                   
c  23) gamma  +  n  --->  pi+  +  pi-  +  pi0  +  n                   
c  24) gamma  +  n  --->  pi+  +  pi-  +  pi-  +  p                   
c  25) gamma  +  p  --->  phi  +  p                                   
c  26) gamma  +  n  --->  phi  +  n                                   
c  27) gamma  +  p  --->  pi+  +  pi+  +  pi-  +  pi-  +  p     
c  28) gamma  +  p  --->  pi+  +  pi+  +  pi-  +  pi-  +  pi0  +  p  
c  29) gamma  +  p  --->  pi+  +  pi+  +  pi+  +  pi-  +  pi-  +  n  
c  30) gamma  +  p  --->  lambda  +  kappa+                
c  31) gamma  +  p  --->  kappa+  +  kappa-  +  p
c  32) gamma  +  p  --->  eta  +  p
c  33) gamma  +  p  --->  f2  +  p
c  34) gamma  +  n  --->  f2  +  n
c
c  51) gamma  +  3He  --->  p + p + n
c  52) gamma  +  3He  --->  p + p + n + pi0
c  53) gamma  +  3He  --->  p + p + p + pi-  
c  54) gamma  +  3He  --->  p + p + n + pi+  
c  55) gamma  +  3He  --->  3He + pi0  
c  56) gamma  +  3He  --->  3H  + pi+  
c  57) gamma  +  3He  --->  d  + p 
c  58) gamma  +  3He  --->  d  + p + pi0 
c  59) gamma  +  3He  --->  d  + n + pi+ 
c  60) gamma  +  d    --->  p  + n 
c  61) gamma  +  d    --->  p  + n + pi0
c  62) gamma  +  d    --->  p  + p + pi-
c  63) gamma  +  d    --->  n  + p + pi+
c  70) gamma  +  pp   --->  p  + p
c  71) gamma  +  pp   --->  p  + p + pi0
c  72) gamma  +  pp   --->  p  + n + pi+
c  80) gamma  +  d    --->  dibaryon + pi-


c  101) delta++ decay 
c  102) delta+  decay 
c  103) delta0  decay 
c  104) delta-  decay 
c  105) rho+    decay 
c  106) rho0    decay 
c  107) rho-    decay 
c  109) phi     decay 
c  110) lambda  decay
c  112) f2      decay 

	implicit none

	integer nchan
	parameter (nchan = 100 )
	character*12 partic(nchan,0:8)
	integer numpart(nchan),jch,n
	common/deca/numpart,partic
	do jch=1,nchan
	   numpart(jch)=0
	   do n=0,8
	      partic(jch,n)='none'
	   enddo
	enddo

	partic(1,0)='proton'
	partic(1,1)='pi+'
	partic(1,2)='neutron'
	numpart(1)=2

	partic(2,0)='proton'
	partic(2,1)='pi0'
	partic(2,2)='proton'
	numpart(2)=2

	partic(3,0)='neutron'
	partic(3,1)='pi-'
	partic(3,2)='proton'
	numpart(3)=2

	partic(4,0)='neutron'
	partic(4,1)='pi0'
	partic(4,2)='neutron'
	numpart(4)=2

	partic(5,0)='proton'
	partic(5,1)='pi-'
	partic(5,2)='delta++'
	numpart(5)=2

	partic(6,0)='proton'
	partic(6,1)='pi0'
	partic(6,2)='delta+'
	numpart(6)=2

	partic(7,0)='proton'
	partic(7,1)='pi+'
	partic(7,2)='delta0'
	numpart(7)=2

	partic(8,0)='neutron'
	partic(8,1)='pi-'
	partic(8,2)='delta+'
	numpart(8)=2

	partic(9,0)='neutron'
	partic(9,1)='pi0'
	partic(9,2)='delta0'
	numpart(9)=2

	partic(10,0)='neutron'
	partic(10,1)='pi+'
	partic(10,2)='delta-'
	numpart(10)=2

	partic(11,0)='proton'
	partic(11,1)='rho0'
	partic(11,2)='proton'
	numpart(11)=2

	partic(12,0)='proton'
	partic(12,1)='rho+'
	partic(12,2)='neutron'
	numpart(12)=2

	partic(13,0)='neutron'
	partic(13,1)='rho-'
	partic(13,2)='proton'
	numpart(13)=2

	partic(14,0)='neutron'
	partic(14,1)='rho0'
	partic(14,2)='neutron'
	numpart(14)=2

	partic(15,0)='proton'
	partic(15,1)='pi+'
	partic(15,2)='pi-'
	partic(15,3)='proton'
	numpart(15)=3

	partic(16,0)='proton'
	partic(16,1)='pi0'
	partic(16,2)='pi+'
	partic(16,3)='neutron'
	numpart(16)=3

	partic(17,0)='neutron'
	partic(17,1)='pi+'
	partic(17,2)='pi-'
	partic(17,3)='neutron'
	numpart(17)=3

	partic(18,0)='neutron'
	partic(18,1)='pi0'
	partic(18,2)='pi-'
	partic(18,3)='proton'
	numpart(18)=3

	partic(19,0)='proton'
	partic(19,1)='omega'
	partic(19,2)='proton'
	numpart(19)=2

	partic(20,0)='neutron'
	partic(20,1)='omega'
	partic(20,2)='neutron'
	numpart(20)=2

	partic(21,0)='proton'
	partic(21,1)='pi+'
	partic(21,2)='pi-'
	partic(21,3)='pi0'
	partic(21,4)='proton'
	numpart(21)=4

	partic(22,0)='proton'
	partic(22,1)='pi+'
	partic(22,2)='pi+'
	partic(22,3)='pi-'
	partic(22,4)='neutron'
	numpart(22)=4

	partic(23,0)='neutron'
	partic(23,1)='pi+'
	partic(23,2)='pi-'
	partic(23,3)='pi0'
	partic(23,4)='neutron'
	numpart(23)=4

	partic(24,0)='neutron'
	partic(24,1)='pi+'
	partic(24,2)='pi-'
	partic(24,3)='pi-'
	partic(24,4)='proton'
	numpart(24)=4

	partic(25,0)='proton'
	partic(25,1)='phi'
	partic(25,2)='proton'
	numpart(25)=2

	partic(26,0)='neutron'
	partic(26,1)='phi'
	partic(26,2)='neutron'
	numpart(26)=2

	partic(27,0)='proton'
	partic(27,1)='pi+'
	partic(27,2)='pi+'
	partic(27,3)='pi-'
	partic(27,4)='pi-'
	partic(27,5)='proton'
	numpart(27)=5

	partic(28,0)='proton'
	partic(28,1)='pi+'
	partic(28,2)='pi+'
	partic(28,3)='pi-'
	partic(28,4)='pi-'
	partic(28,5)='pi0'
	partic(28,6)='proton'
	numpart(28)=6

	partic(29,0)='proton'
	partic(29,1)='pi+'
	partic(29,2)='pi+'
	partic(29,3)='pi+'
	partic(29,4)='pi-'
	partic(29,5)='pi-'
	partic(29,6)='neutron'
	numpart(29)=6

	partic(30,0)='proton'
	partic(30,1)='kappa+'
	partic(30,2)='lambda'
	numpart(30)=2

	partic(31,0)='proton'
	partic(31,1)='kappa+'
	partic(31,2)='kappa-'
	partic(31,3)='proton'
	numpart(31)=3

	partic(32,0)='proton'
	partic(32,1)='eta'
	partic(32,2)='proton'
	numpart(32)=2

	partic(33,0)='proton'
	partic(33,1)='f2'
	partic(33,2)='proton'
	numpart(33)=2

	partic(34,0)='neutron'
	partic(34,1)='f2'
	partic(34,2)='neutron'
	numpart(34)=2

	partic(35,0)='proton'
	partic(35,1)='omega_mod'
	partic(35,2)='proton'
	numpart(35)=2

        partic(40,0)='proton'
	partic(40,1)='proton'
	partic(40,2)='antiproton'
	partic(40,3)='proton'
	numpart(40)=3

	partic(41,0)='proton'
	partic(41,1)='neutron'
	partic(41,2)='antineutron'
	partic(41,3)='proton'
	numpart(41)=3

	partic(42,0)='proton'
	partic(42,1)='proton'
	partic(42,2)='antiproton'
	partic(42,3)='proton'
	numpart(42)=3

	partic(43,0)='proton'
	partic(43,1)='neutron'
	partic(43,2)='antineutron'
	partic(43,3)='proton'
	numpart(43)=3

	partic(51,0)='he3'
	partic(51,1)='proton'
	partic(51,2)='proton'
	partic(51,3)='neutron'
	numpart(51)=3

	partic(52,0)='he3'
	partic(52,1)='proton'
	partic(52,2)='proton'
	partic(52,3)='neutron'
	partic(52,4)='pi0'
	numpart(52)=4

	partic(53,0)='he3'
	partic(53,1)='proton'
	partic(53,2)='proton'
	partic(53,3)='proton'
	partic(53,4)='pi-'
	numpart(53)=4

	partic(54,0)='he3'
	partic(54,1)='proton'
	partic(54,2)='neutron'
	partic(54,3)='neutron'
	partic(54,4)='pi+'
	numpart(54)=4

	partic(55,0)='he3'
	partic(55,1)='pi0'
	partic(55,2)='he3'
	numpart(55)=2

	partic(56,0)='he3'
	partic(56,1)='pi+'
	partic(56,2)='h3'
	numpart(56)=2

	partic(57,0)='he3'
	partic(57,1)='deuteron'
	partic(57,2)='proton'
	numpart(57)=2

	partic(58,0)='he3'
	partic(58,1)='deuteron'
	partic(58,2)='proton'
	partic(58,3)='pi0'
	numpart(58)=3

	partic(59,0)='he3'
	partic(59,1)='deuteron'
	partic(59,2)='neutron'
	partic(59,3)='pi+'
	numpart(59)=3

	partic(60,0)='deuteron'
	partic(60,1)='proton'
	partic(60,2)='neutron'
	numpart(60)=2

	partic(61,0)='deuteron'
	partic(61,1)='proton'
	partic(61,2)='neutron'
	partic(61,3)='pi0'
	numpart(61)=3

	partic(62,0)='deuteron'
	partic(62,1)='proton'
	partic(62,2)='proton'
	partic(62,3)='pi-'
	numpart(62)=3

	partic(63,0)='deuteron'
	partic(63,1)='neutron'
	partic(63,2)='neutron'
	partic(63,3)='pi+'
	numpart(63)=3

	partic(70,0)='diproton'
	partic(70,1)='proton'
	partic(70,2)='proton'
	numpart(70)=2

	partic(71,0)='diproton'
	partic(71,1)='proton'
	partic(71,2)='proton'
	partic(71,3)='pi0'
	numpart(71)=3

	partic(72,0)='diproton'
	partic(72,1)='proton'
	partic(72,2)='neutron'
	partic(72,3)='pi+'
	numpart(72)=3

	partic(80,0)='deuteron'
	partic(80,1)='dibaryon'
	partic(80,2)='pi-'
	numpart(80)=2

	partic(44,0)='proton'
	partic(44,1)='ppbar'
	partic(44,2)='proton'
	numpart(44)=2

	partic(45,0)='proton'
	partic(45,1)='pp'
	partic(45,2)='antiproton'
	numpart(45)=2

	return 
	end

							  ! *** end decay ***
