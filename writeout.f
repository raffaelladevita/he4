
	subroutine writeout(jch,now,p,rm,q,npcode,eglab,p0i)
 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c
c 									    c
c    this subroutine writes on output control file the variables.	    c
c 									    c
c---------------------------------------------------------------------------c

	implicit none

	integer j,k,now,jch,npcode(10)
	integer loop
        real p(5,10),rm(10),q(10),eglab,p0i(4),delta(4)
	data loop/0/

	 loop=loop+1
	 write (11,118) loop

         write (11,*) ' '
         if (jch.eq.1)
     &   write (11,*) ' gamma  +  p  --->  pai+  +  n'
         if (jch.eq.2)
     &   write (11,*) ' gamma  +  p  --->  pai0  +  p'
         if (jch.eq.3)
     &   write (11,*) ' gamma  +  n  --->  pai-  +  p'
         if (jch.eq.4)
     &   write (11,*) ' gamma  +  n  --->  pai0  +  n'
         if (jch.eq.5) 
     &   write (11,*) ' gamma  +  p  --->  delta++ +  pai-'
         if (jch.eq.6)
     &   write (11,*) ' gamma  +  p  --->  delta+  +  pai0'
         if (jch.eq.7)
     &   write (11,*) ' gamma  +  p  --->  delta0  +  pai+'
         if (jch.eq.8)
     &   write (11,*) ' gamma  +  n  --->  delta+  +  pai-'
         if (jch.eq.9)
     &   write (11,*) ' gamma  +  n  --->  delta0  +  pai0'
         if (jch.eq.10)
     &   write (11,*) ' gamma  +  n  --->  delta-  +  pai+'
         if (jch.eq.11)
     &   write (11,*) ' gamma  +  p  --->  rho0  +  p'
         if (jch.eq.12)
     &   write (11,*) ' gamma  +  p  --->  rho+  +  n'
         if (jch.eq.13)
     &   write (11,*) ' gamma  +  n  --->  rho-  +  p'
         if (jch.eq.14)
     &   write (11,*) ' gamma  +  n  --->  rho0  +  n'
         if (jch.eq.15)
     &   write (11,*) ' gamma  +  p  --->  pi+  +  pi-  +  p'
         if (jch.eq.16)
     &   write (11,*) ' gamma  +  p  --->  pi0  +  pi+  +  n'
         if (jch.eq.17)
     &   write (11,*) ' gamma  +  n  --->  pi+  +  pi-  +  n'
         if (jch.eq.18)
     &   write (11,*) ' gamma  +  n  --->  pi0  +  pi-  +  p'
         if (jch.eq.19)
     &   write (11,*) ' gamma  +  p  --->  omega  +  p'
         if (jch.eq.20)
     &   write (11,*) ' gamma  +  n  --->  omega  +  n'
         if (jch.eq.21)
     &   write (11,*) ' gamma  +  p  --->  pi+  +  pi-  +  pi0  +  p'
         if (jch.eq.22)
     &   write (11,*) ' gamma  +  p  --->  pi+  +  pi+  +  pi-  +  n'
         if (jch.eq.23)
     &   write (11,*) ' gamma  +  n  --->  pi+  +  pi-  +  pi0  +  n'
         if (jch.eq.24)
     &   write (11,*) ' gamma  +  n  --->  pi+  +  pi-  +  pi-  +  p'
         if (jch.eq.25)
     &   write (11,*) ' gamma  +  p  --->  phi  +  p'
         if (jch.eq.26)
     &   write (11,*) ' gamma  +  n  --->  phi  +  n'
         if (jch.eq.27) write (11,*)
     &   ' gamma  +  p  --->  pi+  +  pi+  +  pi-  +  pi-  +  p'
         if (jch.eq.28) write (11,*) 
     &   ' gamma  +  p  --->  pi+  +  pi+  +  pi-  +  pi-  +  pi0  +  p'
         if (jch.eq.29) write (11,*) 
     &   ' gamma  +  p  --->  pi+  +  pi+  +  pi+  +  pi-  +  pi-  +  n'
         if (jch.eq.30) write (11,*) 
     &   ' gamma  +  p  --->  lambda  +  kappa+'
         if (jch.eq.31) write (11,*) 
     &   ' gamma  +  p  --->  p  +  kappa+  +  kappa-'
         if (jch.eq.32) write (11,*) 
     &   ' gamma  +  p  --->  p  +  eta'
         if (jch.eq.33) write (11,*) 
     &   ' gamma  +  p  --->  p  +  f2'
         if (jch.eq.34) write (11,*) 
     &   ' gamma  +  n  --->  n  +  f2'
         if (jch.eq.51) write (11,*) 
     &   ' gamma  +  3He  --->  p + p + n'
         if (jch.eq.52) write (11,*) 
     &   ' gamma  +  3He  --->  p + p + n + pi0'
         if (jch.eq.53) write (11,*) 
     &   ' gamma  +  3He  --->  p + p + p + pi-'
         if (jch.eq.54) write (11,*) 
     &   ' gamma  +  3He  --->  p + n + n + pi+'
         if (jch.eq.55) write (11,*) 
     &   ' gamma  +  3He  --->  3He + pi0'
         if (jch.eq.56) write (11,*) 
     &   ' gamma  +  3He  --->  3H + pi+'
         if (jch.eq.57) write (11,*) 
     &   ' gamma  +  3He  --->  d + p'
         if (jch.eq.58) write (11,*) 
     &   ' gamma  +  3He  --->  d + p + pi0'
         if (jch.eq.59) write (11,*) 
     &   ' gamma  +  3He  --->  d + n + pi+'
         if (jch.eq.60) write (11,*) 
     &   ' gamma  +  d  --->  p + n'
         if (jch.eq.61) write (11,*) 
     &   ' gamma  +  d  --->  p + n + pi0'
         if (jch.eq.62) write (11,*) 
     &   ' gamma  +  d  --->  p + p + pi-'
         if (jch.eq.63) write (11,*) 
     &   ' gamma  +  d  --->  p + n + pi+'
         if (jch.eq.70) write (11,*) 
     &   ' gamma  +  pp  --->  p + p'
         if (jch.eq.71) write (11,*) 
     &   ' gamma  +  pp  --->  p + p + pi0'
         if (jch.eq.72) write (11,*) 
     &   ' gamma  +  pp  --->  p + n + pi+'
         if (jch.eq.80) write (11,*) 
     &   ' gamma  +  d  --->  dibaryon + pi-'

       	 write (11,*) ' '
         write (11,110) eglab
         write (11,*) ' '
         write (11,111) p0i
         write (11,*) ' '
         write (11,112)
         write (11,*) ' '
	 do k=1,now
           write (11,113,err=100)
     +     q(k),npcode(k),p(1,k),p(2,k),p(3,k),p(4,k),rm(k)
	 end do
	write(11,*) ' '
	write(11,114)
	write(11,115)
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
	write(11,116) delta
	write(11,*) ' '
	write(11,117)
	write(11,*) ' '

	return

 100	write(6,*) 'ERROR IN WRITING'
 110	format(2x,'tagged photon energy :',f10.5,' GeV')
 111	format(2x,'fermi motion 4-momentum :',4f10.5)
 112	format(4x,' q ',4x,'part.',6x,'px',8x,'py',8x,'pz',9x,' e',8x,' m')
 113	format(4x,f3.0,4x,i5,5(2x,f9.5))
 114	format(2x,'4-momentum conservation check :')
 115	format(4x,'pxf-pxi',10x,'pyf-pyi',10x,'pzf-pzi',10x,'Ef-Ei')
 116	format(1x,4(e12.5,5x))
 117	format(70('*'))
 118	format(2x,'event number',i5,/)

	end 

 
c						       c *** end writeout ***







