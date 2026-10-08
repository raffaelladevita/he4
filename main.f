c***************************************************************************c

	program he3
 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c 									    c
c abstract:								    c
c 									    c
c	 this program simulates the photon interaction with helium3 or      c
c        nucleon in the case of production of real pions in the		    c
c        final state.						            c
c	 it's possible to select a single channel (with one or two pions    c
c	 in final state), or to consider all the processes weighted with    c
c	 their cross sections.						    c
c	 it's also possible to neglect some channels.			    c
c	 this program calculates the differential and total cross section   c
c	 for the sigle and double pion photoproduction.			    c
c	 it is based on the work of r.l. walker and w.j.metcalf (phys.rev.  c
c	 182 (1969) 1729 and nucl.phys. b76 (1974) 253).		    c
c	 the analysis is made in terms of a simple model in which the       c
c	 photoproduction amplitude consist of three separate contributions: c
c	 the born approximation, expressed in terms of cgln coefficients,   c
c        breit-wigner resonances and additional background contributions    c
c        in the low partial waves.			  		    c
c	 for double pion photoproduction the born and resonance part of     c
c	 total cross section is calculated.				    c
c	 the channels studied are delta-pai, rho-nucleon and phase space.   c 
c									    c
c	 the output files are:						    c
c									    c
c	 - paw file to see impulse distributions of initial state and       c
c	   final state particles.					    c
c									    c
c	 - bos file as input to fastmc code.				    c
c									    c
c 									    c
c author(s):								    c
c 									    c
c	 p. corvisiero  l. mazzaschi					    c
c 									    c
c creation date:							    c
c 									    c
c	 26 - 10 - 1992							    c
c 									    c
c common blocks:							    c
c  									    c
c     /control/								    c
c     /control_gen/							    c
c     /pawc/								    c
c     /quest/								    c
c  									    c
c modification history:							    c
c 									    c
c	 date     | name  | description					    c
c ----------------+-------+-------------------------------------------------c
c  4 -  1 - 1993  |  l.m. | change of bos file filling and paw file        c
c		  |	  | filling. inclusion of new channels for he3.     c
c ----------------+-------+-------------------------------------------------c
c [change_entry]							    c
c---------------------------------------------------------------------------c

c	implicit none

c	integer nw,ibid,m0,m8,m1,j2,mevt,run,istatus1
	integer lrecl
	integer nfpaw,nfbos,pawid
	integer bosout,hbookout
	common/output/bosout,hbookout

      integer*4 maxpages
      parameter ( maxpages= 100000)
     
      integer*4 iquest
      integer*4 ipawc
      common /quest/iquest(100)
      common /pawc/ipawc(maxpages*128)

	include "bcs.inc"
	
	data nfpaw /13/
	data nfbos /12/
	data pawid /1/
	integer nzlink
        parameter (nzlink = 32000)

	integer j,nch
      character filename*(40),char(40)*1
      character fileinp*(40),filebos*(40),fileout*(40),filepaw*(40)
      character inp(40)*1,bosf(40)*1,out(40)*1,paw(40)*1
      equivalence(filename,char)
      equivalence(fileinp,inp)
      equivalence(filebos,bosf)
      equivalence(fileout,out)
      equivalence(filepaw,paw)

      common/filenames/fileinp,filebos,fileout,filepaw

	integer maxp,np
	real px,py,pz,p,e,m,q,beta,code,ch
	integer ierr,icycle
	real w,eg,evt,teta,phi,t_mdst,M_nn,M_nnbar

	integer np_nd,np_dc
	real p_cm,th_cm,m_nd
	real m_dc,p_dc,th_dc,ph_dc

      parameter (maxp = 16)
      common /cwn1/ np,
     +              px(maxp), py(maxp),  pz(maxp), e(maxp),
     +              p(maxp),  beta(maxp), m(maxp),
     +              q(maxp),  code(maxp) 
      common /cwncm/np_nd,m_nd(maxp),p_cm(maxp),th_cm(maxp)
      common /cwndc/np_dc,m_dc(maxp),p_dc(maxp),th_dc(maxp),ph_dc(maxp) 
      common /cwn2/ ch, w, eg, evt
      common /cwn3/ teta(maxp), phi(maxp)
      common /ct/ t_mdst
      common /mnnbar/ M_nn,M_nnbar

	integer istat,id,lun,nh
	data lrecl/1024/ 

      call getarg(1,filename)
      do j=1,40
        if(char(j).eq.' ') goto 2
      enddo
 2    continue
      char(j)='.'
      nch=j
      do j=1,nch
         inp(j)=char(j)
         bosf(j)=char(j)
         out(j)=char(j)
         paw(j)=char(j)
      enddo
      inp(nch+1)='i'
      inp(nch+2)='n'
      inp(nch+3)='p'
      bosf(nch+1)='b'
      bosf(nch+2)='o'
      bosf(nch+3)='s'
      out(nch+1)='o'
      out(nch+2)='u'
      out(nch+3)='t'
      paw(nch+1)='h'
      paw(nch+2)='b'
      paw(nch+3)='o'
      paw(nch+4)='o'
      paw(nch+5)='k'

        iquest(10) = nzlink
	nh=maxpages*128
        call hlimit(nh)


c ***   reading and general initialization
	call h_read
	call sig_data
	call sig_max
	call decay
c ***   opening control outp file

c ***   opening filebos
	if (bosout.eq.1) then
c          open (unit=11,file=fileout,type='unknown')
          open (unit=11,file=fileout)
c	  inizializzazione e filling BOS banks ....... 
          NW     = 11            ! number of colums in event bank
          IBID   = 12            ! BOS output device number
          M0     = 0             ! number of header bank
          M8     = 8             ! number of colums in header bank
          M1     = 1             ! number of rows in header bank
          J2     = 1             ! number of rows in event bank
          MEVT   = 1             ! event number        
          RUN    = 1       
          call BOS(iw,Nbcs)
          call BKFMT("HEAD","I")
          call BKFMT('MCEV','I')
          call BKFMT('MCVX','(4F,I)')     ! MC vertex parameters 
          call BKFMT('MCTK','(6F,5I)')
          call BKFMT('TAGR','(3F,I)')
          call BLIST(iw,"E=","HEAD")
          call BLIST(iw,'E+','MCEV')
          call BLIST(iw,'E+','MCVX')
          call BLIST(iw,'E+','MCTK')
          call BLIST(iw,'E+','TAGR')	  
          CALL FPARM(
     >    'OPEN UNIT=12 FILE='//filebos(1:index(filebos,'.')+3)//' WRITE RECL=32760 '//
     >    'ACTION=WRITE STATUS=NEW FORM=BINARY')
	endif

c ***   booking ntuple
	if (hbookout.eq.1) then
	   id=1
	   lun=1
           call hropen(lun,'cebaf',filepaw,'n',lrecl,istat)
           call hbnt(id,'clas',' ')
           call hbset('bsize',1024,ierr)
           call hbname(id,'evnt', ch   , 'ch:r')
           call hbname(id,'evnt', w    , 'w:r')
           call hbname(id,'evnt', eg   , 'eg:r')
           call hbname(id,'evnt', evt  , 'evt:r')
           call hbname(id,'kine', np   , 'np[0,16]:i')
           call hbname(id,'kine', px   , 'px(np):r'  )
           call hbname(id,'kine', py   , 'py(np):r'  )
           call hbname(id,'kine', pz   , 'pz(np):r'  )
           call hbname(id,'kine', e    , 'e(np):r'  )
           call hbname(id,'kine', p    , 'p(np):r'  )
           call hbname(id,'kine', beta , 'beta(np):r'  )
           call hbname(id,'kine', m    , 'm(np):r' )
           call hbname(id,'kine', teta , 'teta(np):r' )
           call hbname(id,'kine', phi  , 'phi(np):r' )
           call hbname(id,'kine', q    , 'q(np):r'  )
           call hbname(id,'kine', code , 'code(np):r'  )
           call hbname(id,'kine',t_mdst, 't_mdst:r' )
           call hbname(id,'kine',M_nn, 'M_nn:r' )
           call hbname(id,'kine',M_nnbar,'M_nnbar:r' )
           call hbname(id,'kine', np_nd, 'np_nd[0,16]:i'  )
	   call hbname(id,'kine', m_nd , 'm_nd(np_nd):R'  )
           call hbname(id,'kine', p_cm , 'p_cm(np_nd):r'  )
           call hbname(id,'kine', th_cm, 'th_cm(np_nd):r'  )

           call hbname(id,'kine', np_dc, 'np_dc[0,16]:i'  )
	   call hbname(id,'kine', m_dc , 'm_dc(np_dc):R'  )
           call hbname(id,'kine', p_dc , 'p_dc(np_dc):r'  )
           call hbname(id,'kine', th_dc, 'th_dc(np_dc):r'  )
           call hbname(id,'kine', ph_dc, 'ph_dc(np_dc):r'  )

	endif

c ***   process simulation
	call mc_clas

c ***   closing filebos
	if(bosout.eq.1) then
	   call fwbos(iw,ibid,'0',istatus1)	
	   call fclos()
	   close(3)
	endif

c ***   closing hbook file
	if(hbookout.eq.1) then
	   call hrout(id,icycle,' ')
	   call hrend('cebaf')
	endif

	stop
	end

							  ! ***  end main *** 
c***************************************************************************c


