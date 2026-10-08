c***************************************************************************c  

	subroutine h_read
 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c functional description:						    c
c 									    c
c    this subroutine reads the input preferences, the tabulated partial     c
c    waves decomposition for resonant part and for background part, and     c
c    the breit-wigner parameters.                                           c
c									    c
c 									    c
c formal parameters:							    c
c  									    c
c     target:								    c
c	  'proton'              proton					    c
c	  'neutron'             neutron					    c
c	  'he3'                 helium3			                    c
c	  'd'  or 'deuterium'   deuterium	         		    c
c	  'h3' or 'tritium'     tritium                                     c
c	  'he4'                 helium4			                    c
c									    c
c     eglab:								    c
c     	  laboratory photon energy (gev)				    c
c									    c
c     fwhm:								    c
c	  tagging channel resolution (%)				    c
c									    c
c     nchain	  							    c
c         number of present channels					    c
c									    c
c     ichain(nchain)							    c
c	  present channel numbers					    c
c									    c
c									    c
c     iphas								    c
c	selection index							    c
c	1 = weighted distribution					    c
c	0 = phase space							    c
c     									    c
c     nevent								    c
c	  number of events						    c
c  									    c
c     nrandom								    c
c	  random number generators					    c
c  									    c
c---------------------------------------------------------------------------c

	implicit none

	integer nchan,nchan1,nchan2,ndim
	parameter (nchan  = 100 )
	parameter (nchan1 = 101 )
	parameter (nchan2 = 200 )
	parameter (ndim   = 201 )
        real ebeam,ethresh,e1tag,e2tag,fraz
	real ft1,ft2
	real q,rran
	character*12 target
	real eglab,fwhm
	integer ncheck
	integer iniz
	integer dummy_int
	real dummy_real
	real  Eaverage,Eboundary
	real bradius,tlenght,toffset
	common /TAGR/Eaverage(767),Eboundary(767)
	common/check/ncheck
        common /target/target,bradius,tlenght,toffset
	common /beam1/eglab,fwhm
	common /beam/ ebeam,ethresh,e1tag,e2tag,fraz
	integer bosout,hbookout
        common/output/bosout,hbookout
	integer flag_acc,zp_final_chan,n_ks,n_kl,k0_final_chan,zp_ang_dist
	real    I_torus
	common/acc_flags/flag_acc,I_torus
	common/zetap_flags/zp_final_chan,n_ks,n_kl,k0_final_chan,zp_ang_dist

c                ***aggiunta marco***

c  def una flag per la omega canale 19  (se = 0 =>il programma funziona come prc  ima,se = 1 =>il prog estrae con funzione continua)

        integer ds_omega_fn
	common/omega/ds_omega_fn

c   ***********************************************************************

	real wr(251),egl(251)
	real sigr(nchan,251)
	common/sigtot/sigr,wr,egl
	integer ith,iph,i
	real th,ph
	real w
        integer nchain,ichain(nchan)
	integer ichain_p(nchan),ichain_n(nchan)
     
        common /channels/nchain,ichain
	real vecgen(0:nchan2)
	common/geneven/vecgen
	integer nevent
	common /control/ nevent

	integer iphas

	real xchan(nchan,50,0:18)
	real xpcle(nchan1:nchan2,0:20,0:40)
	common /sigmas/ xchan, xpcle

	real t_ref,extr_mode
	common /t_refer/ t_ref,extr_mode

	integer n,jch,inw,inth,j,k
c	integer l
	integer nch
	character*40 fileinp,filebos,fileout,filepaw
c	character*60 fileinp,filebos,fileout,filepaw

c    proton channels:
 	data ichain_p/ 1, 2, 5, 6, 7,11,12,15,16,19,
     +                21,22,25,27,28,29,30,31,32,33,
     +                35,36,37,38,40,41,44,45,46,47,
     +                48,86,87,88,89,90,91,92,93,94,
     +                95,96,97,98,99,55*0/
c    neutron channels:
 	data ichain_n/ 3, 4, 8, 9,10,13,14,17,18,20,
     +                23,24,26,34,39,49,84*0/

	common/filenames/fileinp,filebos,fileout,filepaw
	character*120 he4_parms
	character*120 filename
 
 	parameter (he4_parms  = 'GENOVA_PARMS')


ccccccccccccccccccccccccccccccccccccccccccc
ccc da capire...
c	CALL revinm(he3_parms,'~/mystuff/jlab/gamma/he3_xdiff.mtx',filename)
c        open (unit=55,file=filename,status='old',form='formatted')
c     &          access='sequential',form='unformatted')
cccccccccccccccccccccccccccccccccccccccccccc



c       vecgen(0) = 0. ALWAYS (by defintion):
c       it means that the resonance decay is phase-space type
c	do not redefine it to 1. !!!
	data vecgen/ndim*0./

c       initializing disang matrices....
	do jch=1,nchan
	 do inw=1,50
	  do inth=0,18
	   xchan(jch,inw,inth)=0.
	  enddo
	 enddo
	enddo
	do jch=nchan+1,2*nchan
	 do ith=0,20
	  do iph=0,40
	   xpcle(jch,ith,iph)=0.
	  enddo
	 enddo
	enddo

c ***   opening otput control file....
c	open(unit=11,file=fileout,type='unknown')
	open(unit=11,file=fileout)

c ***   reading input preferences
c	open (unit=1,file=fileinp,type='unknown')
	open (unit=1,file=fileinp)
	read (1,*) target
	write (11,2323)
	if (target.eq.'proton') write (11,*)'target = proton' 
	if (target.eq.'neutron') write (11,*)'target = neutron' 
	if (target.eq.'he3') write (11,*)'target = he3' 
	if (target.eq.'h3') write (11,*)'target = tritium' 
	if (target.eq.'tritium') write (11,*)'target = tritium' 
	if (target.eq.'d') write (11,*)'target = deuterium' 
	if (target.eq.'deuterium') write (11,*)'target = deuterium' 
	if (target.eq.'he4') write (11,*)'target = he4' 

	read (1,*) bradius
	write (11,2323)
	write (11,*) 'beam radius = ',bradius

	read (1,*) tlenght
	write (11,2323)
	write (11,*) 'target lenght = ',tlenght

	read (1,*) toffset
	write (11,2323)
	write (11,*) 'target offset = ',toffset

	read (1,*) ebeam
	write (11,2323)
	write (11,*) 'maximum beam energy = ',ebeam

	read (1,*) ethresh
	write (11,2323)
	write (11,*) 'hardware CLAS threshold (GeV) = ',ethresh

	read (1,*) ft1,ft2
	write (11,2323)
	write (11,*) 'min. & max. tagging percent energy = ',ft1,ft2,' % '

	read (1,*) fraz
	write (11,2323)
	write (11,*) '% of good events = ',fraz,' % '

	read (1,*) fwhm
	write (11,2323)
	write (11,*) 'tagging resolution = ',fwhm,' % '

	write (11,2323)
	read (1,*,err=777) (ichain(n), n=1,100)
 777	continue
c	jch = -1 => all proton  channels
c	jch = -2 => all neutron channels
  	if(ichain(1).eq.-1) then
 	   nchain = 22
 	   do j=1,nchain
 	     ichain(j)=ichain_p(j)
 	   enddo
  	elseif(ichain(1).eq.-2) then
 	   nchain = 14
 	   do j=1,nchain
 	     ichain(j)=ichain_n(j)
 	   enddo
 	elseif(ichain(1).gt.0) then
 	   nchain=n-1
 	endif
	write (11,*) '# channels in = ',nchain
	write(11,2325)(ichain(n), n=1,nchain)
 2325	format(10i4)
	write (11,2323)
	read (1,*) t_ref
	write (11,2323)
	write (11,*) 't_limit =     ',t_ref
	read (1,*) iphas
	write (11,2323)
	if (iphas.eq.1) write (11,*) 'weighted distribution'
	if (iphas.eq.0) write (11,*) 'phase space distribution'
	write (11,2323)
	read (1,*) nevent
	write (11,*) 'number of total events = ',nevent
	write (11,2323)
	write (11,2323)
	read (1,*) ncheck
	if(ncheck.le.0) then
	   close(unit=11)
	else
	   write (11,*) 'start controlled events',ncheck
	   write(11,2323)
	   write(11,2324)
	endif

c                  ****      aggiunta marco  ****

c     legge la flag ds_omega_fn e scrive nel file output di check quale tipo di
c     estrazione ho usato
 
	read(1,*) ds_omega_fn
	if (ds_omega_fn.eq.0) write(11,*) 'omega standard extraction'
        if (ds_omega_fn.eq.1) write(11,*) 'omega ds function extraction'

	read(1,*) extr_mode
	if (extr_mode.eq.1) write(11,*) 'omega extraction according f(w,th)'
        if (extr_mode.eq.2) write(11,*) 'omega extraction according f(w,-t)'

c  ***********************************************************************

	read (1,*) iniz
	read (1,*) bosout
	if (bosout.eq.1) then
	   write(11,*) 'BOS file written'
	else
	   write(11,*) 'BOS file not written'
	endif
	read (1,*) hbookout
	if (hbookout.eq.1) then
	   write(11,*) 'HBOOK file written'
	else
	   write(11,*) 'HBOOK file not written'
	endif


c  ***********************************************************************
c--- ADDITION FOR ZETA+ ANALYSIS (R. De Vita 07/19/03)
	read (1,*) flag_acc      ! generate within fiducial cuts (1) or not (0)
	read (1,*) I_torus       ! torus current value (used for fiducial cuts)
	if(flag_acc.eq.1) then
	   write(11,*) 'Events will be generated inside the fiducial cuts'
	else
	   write(11,*) 'Events will be generated in the entire phase space'
	endif
	read (1,*) zp_final_chan ! fix final channel for zeta+ and related channels
	if(zp_final_chan.eq.1) then
	   write(11,*) 'only decay into k+k0 in final state'
	elseif(zp_final_chan.eq.2) then
	   write(11,*) 'only decay into k0k0 in final state'
	elseif(zp_final_chan.eq.3) then
	   write(11,*) 'only decay into k+k- in final state'
	else
	   write(11,*) 'all possible decay channels will be included'
	endif
	read (1,*) zp_ang_dist   ! fix angular zp angular distribution
	read (1,*) k0_final_chan ! fix final channel k0
	if(k0_final_chan.eq.1) then
	   write(11,*) 'only decay into pi+ pi- final state'
	else
	   write(11,*) 'all possible decays'
	endif
c  ***********************************************************************	

	close (unit=1)


c ***   reading file containing differential cross section

c       in file xdiff.mtx are stored only differential cross
c       sections of those channels needed to subroutine genevt
c       genbod doesn't neeed differential cross sections values...

	CALL revinm(he4_parms,'he4_chan_disang.mtx',filename)
        open(unit=55,file=filename,status='old',form='formatted')
      do nch=1,nchan
       read(55,*,end=333) jch
c      set vecgen(jch) = 1. -> not genbod but genevt ... 
       if(iphas.eq.1) vecgen(jch)=1.
       do inw=1,50
	read(55,*) w
	read(55,*) (xchan(jch,inw,inth),inth=0,18)
       end do
      end do
 333  continue
      close (unit=55)

	CALL revinm(he4_parms,'he4_pcle_disang.mtx',filename)
        open(unit=55,file=filename,status='old',form='formatted')
      do nch=nchan+1,nchan2
       read(55,*,end=334) jch
c      set vecgen(jch) = 1. -> not genbod but genevt ... 
       if(iphas.eq.1) vecgen(jch)=1.
       do ith=0,20
	  do iph=0,40
	     read(55,*) th,ph,xpcle(jch,ith,iph)
          end do
       end do
      enddo
 334  continue
      close (unit=55)
	
	CALL revinm(he4_parms,'he4_xtot.mtx',filename)
        open(unit=55,file=filename,status='old',form='formatted')
	do j=1,100
	   read(55,*,end=335) jch
	   read(55,*) (sigr(jch,k),k=1,251)
	enddo
 335	continue
	close(unit=55)
	do k=1,251
	   sigr(90,k)=sigr(90,k)/1000.
	enddo
c+ Reading TAGGER E-> T E correspondenceble


	CALL revinm(he4_parms,'tagE-boundaries.dat',filename)
        open(unit=55,file=filename,status='old',form='formatted')

7         format (i5,i6,i6,i6,i6,f9.1,f8.6,f8.6,f8.6,f8.6)

          print*, '** dat2flux: Reading energy E boundaries from local file:tagE-boundaries.dat'
	  read(55,7) dummy_int,dummy_int,dummy_int,dummy_int,dummy_int,
     +		dummy_real,Eaverage(1),dummy_real,Eboundary(2),Eboundary(1)
    	  do i=2,766
	    read(55,7) dummy_int,dummy_int,dummy_int,dummy_int,dummy_int,
     +		dummy_real,Eaverage(i),dummy_real,Eboundary(i+1),dummy_real
    	  enddo
	  do i=1,766
	   Eaverage(i)=Eaverage(i)* ebeam
	   Eboundary(i)=Eboundary(i)*ebeam
    	  enddo	   
	close(unit=55)

c	randomizing rran sequence.....
        CALL RLUXGO(4,iniz,0,0)    
c	do j=1,iniz
c	  q=rran()
c	enddo

	e1tag=ebeam*ft1/100.
	e2tag=ebeam*ft2/100.


100	format (a5,2f7.4,i3,i3,f8.4,a6,f9.4,f9.4,f9.4,f9.4)
101	format (a5,9f7.3)
102	format (a5,7f7.3)
2323    format (1x)
2324	format(70('*'))
	
	return
	end

						   ! ***  end h_readinput ***
















