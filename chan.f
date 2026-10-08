c***************************************************************************c

	integer function chan(p)
	
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c
c 									    c
c    this function assigns to each particle or the channel decay            c 
c 									    c
c 									    c
c function value:							    c
c 									    c
c    chan 								    c
c    chan(p) = 0  means that the resonance decay is phase-space type
c---------------------------------------------------------------------------c

	implicit none
	character*12 p
	chan = 0
c	if(p.eq.'delta++')    chan = 101
c	if(p.eq.'delta+')     chan = 102
c	if(p.eq.'delta0')     chan = 103
c	if(p.eq.'delta-')     chan = 104
c	if(p.eq.'rho+')       chan = 105
	if(p.eq.'rho0')       chan = 106
c	if(p.eq.'rho-')       chan = 107
	if(p.eq.'omega')      chan = 108
c	if(p.eq.'phi')        chan = 109
c	if(p.eq.'lambda')     chan = 110
c	if(p.eq.'eta')        chan = 111
c	if(p.eq.'f2')         chan = 112
	if(p.eq.'ppbar')      chan = 113
	if(p.eq.'pp')         chan = 114
	if(p.eq.'zeta+')      chan = 115
	if(p.eq.'s11')        chan = 116
	if(p.eq.'dipion')     chan = 117
c	if(p.eq.'lambda*')    chan = 118

	return
	end

							  ! *** end chan ***
