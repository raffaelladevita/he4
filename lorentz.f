c***************************************************************************c 

	subroutine lorentz(beta,p)
  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c									    c 
c functional description:						    c
c 									    c
c    this subroutine performs lorentz boost.				    c
c									    c 
c---------------------------------------------------------------------------c

	implicit none
	real beta(3),p(4)
	real v(4),b2,scal_prod,gam,betaxp
	integer j

	b2 = scal_prod(beta,beta)
	gam = 1./sqrt(1.-b2)
	betaxp = scal_prod(beta,p)
	v(4) = gam*(p(4)-betaxp)
	do j=1,3
	  v(j) = p(j)+beta(j)*((gam-1.)/b2*betaxp-gam*p(4))
	end do
	do j=1,4
	  p(j) = v(j)
	end do
	return
	end
  
							! *** end lorentz ***
