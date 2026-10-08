c***************************************************************************c

	real function scal_prod(a,b)
	
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c
c 									    c
c    this function performs the product between two four_momentum.	    c
c 									    c
c									    c 
c function value:							    c
c 									    c
c    scal_prod								    c
c 									    c
c---------------------------------------------------------------------------c

	implicit none
	real a(4),b(4)
	integer j

	scal_prod = 0.
	do j=1,3
	  scal_prod = scal_prod+a(j)*b(j)
	end do
	end

						      ! *** end scal_prod ***
