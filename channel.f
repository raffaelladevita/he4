c***************************************************************************c

	subroutine channel(eglab,jch)
  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c 
c 									    c
c    this subroutine selects the interaction channel using total cross	    c
c    section to produce weighted probabilities.				    c
c 									    c
c---------------------------------------------------------------------------c

	implicit none

	integer nchan
	parameter (nchan = 100 )
        character*12 partic(nchan,0:8)
        integer nchain,ichain(nchan)
	integer numpart(nchan)
        common /channels/nchain,ichain
	common/deca/numpart,partic

	real rran
	real eglab
	real q,sigrel(nchan),sigsum,prob(0:nchan)
	real sigtotrif(251),divdif
	integer jch,n,k
	integer j
	real wr(251),egl(251)
	real sigr(nchan,251)

	common/sigtot/sigr,wr,egl

c ***   variable initialization

	sigsum = 0.
	prob(0) = 0.
c	print *, nchain

	do j=1,nchain
	  sigrel(j) = 0.
	  prob(j) = 0.
	end do

	do n=1,nchain
	  do k=1,251
	     sigtotrif(k) = sigr(ichain(n),k)
	  end do
	  sigrel(n) = divdif(sigtotrif,egl,251,eglab,1)
	  sigsum = sigsum+sigrel(n)
	end do
	do n=1,nchain-1
	  sigrel(n) = sigrel(n)/sigsum
	  prob(n) = prob(n-1)+sigrel(n)
	end do

	prob(nchain)=1.
	q = rran()

	do j=1,nchain
	  if (q.gt.prob(j-1).and.q.le.prob(j)) then
	    jch = ichain(j)
	    return
	  end if
	end do
	return		
	end
 
							! *** end channel ***





