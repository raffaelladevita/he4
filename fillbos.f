ccc programma per la scrittura delle BOS-bank
	subroutine fillbos(nstory)
	implicit none
	integer j,ihead,nbank,nr,mevt,imcev,imcvx,imctk,itagr,istatus1
	real p(5,10),rm(10),q(10)
	real true_evt,egamma
	integer sector
	real rran
	integer numpar,jch0,npcode(10)
	integer Eid_from_egamma
	common/control1/true_evt
	common/bosfile1/p,rm,egamma,q
	common/bosfile2/numpar,jch0,npcode
	real c(3,10)
	integer nstory,k
	real I_torus, bit_p,bit_pip,bit_pim,fac_cut,th,ph,vert_z,vert_x,vert_y,vert_r,vert_phi
	logical bit_10,bit_11,bit_12,bit_13,bit_110

	include "bcs.inc"

c+++ Init fid cuts
	I_torus = 3376.
	fac_cut = 1.5
	bit_p = 0.
	bit_pip = 0.
	bit_pim = 0.
	bit_10 = .false.
	bit_11 = .false.
	bit_12 = .false.
	bit_13 = .false.
	bit_110 = .false.
c+++ Extracting vertex position
        vert_r   = rran()*  0.5 ! Vertex radius
        vert_phi = rran()*  6.2832 ! Vertex phi
        vert_x = vert_r*cos(vert_phi)
        vert_y = vert_r*sin(vert_phi)
	vert_z = (rran()-.5)* 17.7 - 0.0 ! Target lenght + offset

	do j=1,numpar
	   if(p(5,j).ne.0.) then
	      do k=1,3
		 c(k,j)=p(k,j)/p(5,j)
	      enddo
	   else
	      do k=1,3
		 c(k,j)=0.
	      enddo
	   endif
c+++ Check of fiducial cuts
	  if((c(3,j).gt.1.or.c(3,j).lt.-1).or. 
     %       (c(1,j)+c(2,j)+c(3,j).eq.0)) then
            th = -1000
           print *, 'EVNT theta angle undefined for particle n. ',j
           continue
          else
	   th = (180./3.1416)*acos(c(3,j))
	  endif	

	  if((c(1,j).gt.1.or.c(1,j).lt.-1.or.
     %       c(2,j).gt.1.or.c(2,j).lt.-1).or.
     %       (c(1,j)+c(2,j)+c(3,j).eq.0)) then
           ph = -1000
           print *, 'EVNT x or y momentum component not well defined for particle n. ',j
           continue
          else		
	   ph = (180./3.1416)*atan(c(2,j)/c(1,j))
	   if(c(1,j).lt.0.and.c(2,j).gt.0) ph = ph + 180.
	   if(c(1,j).lt.0.and.c(2,j).lt.0) ph = ph + 180.
	   if(c(1,j).gt.0.and.c(2,j).lt.0) ph = ph + 360.
	  endif

	if (npcode(j).eq.2212)   ! proton
     &   call pseudo_spa(1,p(5,j),th,(mod(ph+30.,60.)-30.),I_torus,
     &                    fac_cut,sector(ph),bit_p)
	if (npcode(j).eq.211)  ! pip
     &   call pseudo_spa(1,p(5,j),th,(mod(ph+30.,60.)-30.),I_torus,
     &                    fac_cut,sector(ph),bit_pip)
	if (npcode(j).eq.-211)  ! pip
     &   call pseudo_spa(-1,p(5,j),th,(mod(ph+30.,60.)-30.),I_torus,
     &                    fac_cut,sector(ph),bit_pim)
c--- Ended fid cuts check
	enddo




C+++++ 	Writing BOS header
        iHEAD = NBANK("HEAD",0,8,1)           ! BOS HEADER
        iw(iHEAD+1)=2
        iw(iHEAD+2)=1
        iw(iHEAD+3)=nstory
        iw(iHEAD+4)=0
        iw(iHEAD+5)=-2
        iw(iHEAD+6)=0
        iw(iHEAD+7)=7
        iw(iHEAD+8)=0
c-----
C+++++ 	Writing	other MonteCarlo BOS banks	
	NR   = MEVT

        iMCEV = NBANK('MCEV',0,2,1)
        iw(iMCEV+1) = rran()*100000
        iw(iMCEV+2) = rran()*100000

	iMCVX = NBANK('MCVX',0,5,1)
	rw(iMCVX+1) = vert_x
	rw(iMCVX+2) = vert_y
	rw(iMCVX+3) = vert_z
	rw(iMCVX+4) = 0.
	iw(iMCVX+5) = 0
c-----
c+++++ Writing MCTK bank	
        iMCTK = NBANK('MCTK',0,11,numpar+1)
c-----        
c+++++ Looping and writing final state particles variables 
c++++ Filling initial  photon information        
           rw(imctk+1)=0
           rw(imctk+2)=0       
           rw(imctk+3)=0 ! to be changed in 0
           rw(imctk+4)=0
           rw(imctk+5)=egamma   ! In mass place there is the photon enrgy
           rw(imctk+6)=true_evt ! true or random flag
	   iw(imctk+7)=jch0      ! channel number
           iw(imctk+8)=0
           iw(imctk+9)=1
           iw(imctk+10)=0
           iw(imctk+11)=0
	   imctk = imctk + 11       
	DO j = 1,numpar
           rw(imctk+1)=c(1,j)
           rw(imctk+2)=c(2,j)           
           rw(imctk+3)=c(3,j)
           rw(imctk+4)=p(5,j)
           rw(imctk+5)=rm(j)
           rw(imctk+6)=q(j)
	   iw(imctk+7)=npcode(j)
           iw(imctk+8)=0
           iw(imctk+9)=1
           iw(imctk+10)=0
           iw(imctk+11)=0
	   imctk = imctk + 11
        ENDDO
c-----       
c+++++ Writing TAGR bank	
        iTAGR = NBANK('TAGR',0,6,1)
        rw(iTAGR+1)= egamma	! energia fotone
        rw(iTAGR+2)= 0.         ! Time of the photon reconstruc. in the Tagger
        rw(iTAGR+3)= 0.		! Time of the photon after RF correction
        iw(iTAGR+4)= 7          ! Status ( 7 or 15 are Good)
        iw(iTAGR+5)= 1          ! T counter Id
        iw(iTAGR+6)= Eid_from_egamma(egamma)          ! E counter Id 
     
c-----
        CALL FWBOS(IW,12,'E',iSTATUS1)
        CALL BDROP(IW,'E')
        CALL BGARB(IW)
       
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
        return
        end        


	function Eid_from_egamma(eg)
	implicit none
	integer Eid_from_egamma,j
	real eg
	real  Eaverage,Eboundary
	common /TAGR/Eaverage(767),Eboundary(767)
	do j = 1,766
c	 write (*,*)Eboundary(j)
	 
	  Eid_from_egamma = j
	  if (eg.gt.Eboundary(j+1).and.eg.lt.Eboundary(j)) goto 1
        enddo
c	stop
1	return
	end







