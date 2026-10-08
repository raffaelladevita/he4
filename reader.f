 	SUBROUTINE reader(unit,fl,eg,pr,pip,pim,miss2,vser)
	implicit none
	logical fl
        integer jj, unit
        real eg,th_p,ph_p,th_pip,ph_pip,th_pim,ph_pim
	real p_p, p_pip,p_pim,miss2
        real pr(4),pip(4),pim(4),dummy
	real tsne1,tsne2
	real vser(10) ! service vector

	fl=.false.
	goto 99
c	PRINT *, 'OKKK'
c          read(unit,*,*) dummy,eg
          read(unit,*,ERR=5656,END=5656) miss2, eg

c	PRINT *, 'OKKK'
c           print *, 'egamma=', miss2, eg
           DO jj=1,4
            read(unit,*,ERR=5656,END=5656) pr(jj),pip(jj),pim(jj)
c           print *,'p pip pim = ',jj, pr(jj),pip(jj),pim(jj)
           enddo 
	   p_p=sqrt(pr(1)**2+pr(2)**2+pr(3)**2)
	   p_pip=sqrt(pip(1)**2+pip(2)**2+pip(3)**2)
	   p_pim=sqrt(pim(1)**2+pim(2)**2+pim(3)**2)
	   call angles(0,pr(1),pr(2),pr(3),th_p,ph_p)
	   call angles(0,pip(1),pip(2),pip(3),th_pip,ph_pip)
	   call angles(0,pim(1),pim(2),pim(3),th_pim,ph_pim)
c           print *,'P p pip pim =', p_p, p_pip, p_pim
c           print *,'Theta p pip pim =', th_p, th_pip, th_pim
c           print *,'Phi p pip pim =', ph_p, ph_pip, ph_pim

c Reader foooooooor multilayer analysis May 21
c [gamma', 'px', 'py', 'pz', 'pipx', 'pipy', 'pipz', 'pimx', 'pimy', 'pimz', 'pe', 'pipe', 'pime', 'M_{p\pi^+}', 'M_{\pi+\pi^-}', 'tsne1', 'tsne2']
c [gamma', 'px', 'py', 'pz', 'pipx', 'pipy', 'pipz', 'pimx', 'pimy', 'pimz', 'pe', 'pipe', 'pime', 'pim_mass', 'pca1','pca2', 'tsne1', 'tsne2']
	   
 99	   continue
          read(unit,*,ERR=5656,END=5656) eg
     %       ,pr(1),pr(2),pr(3)
     %       ,pip(1),pip(2),pip(3)
     %       ,pim(1),pim(2),pim(3)
     %       ,pr(4),pip(4),pim(4)
     %       ,vser(3) !pi missing mass
     %       ,vser(4),vser(5) !pca1, pca2
     %       ,vser(1),vser(2) !tsne1, tsne2
     %       ,vser(6),vser(7),vser(8) !tsne1, tsne2label1,label2,label3
	  print*, vser(6),vser(7),vser(8)
c	  print *, 'egamma=', eg
c 	  print *, pr(1),pr(2),pr(3),tsne1
	   p_p=sqrt(pr(1)**2+pr(2)**2+pr(3)**2)
	   p_pip=sqrt(pip(1)**2+pip(2)**2+pip(3)**2)
	   p_pim=sqrt(pim(1)**2+pim(2)**2+pim(3)**2)
	   call angles(0,pr(1),pr(2),pr(3),th_p,ph_p)
	   call angles(0,pip(1),pip(2),pip(3),th_pip,ph_pip)
	   call angles(0,pim(1),pim(2),pim(3),th_pim,ph_pim)
c           print *,'P p pip pim =', p_p, p_pip, p_pim
c           print *,'Theta p pip pim =', th_p, th_pip, th_pim
c           print *,'Phi p pip pim =', ph_p, ph_pip, ph_pim
   


	   goto 100
5656    fl=.true.

 100	return


	end      	
      
 
