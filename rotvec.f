c
	subroutine rotvec(v4,pprime)
c
c	trasforma le componenti di un 4-vettore dal sdr s' al sdr s
c	dove s' = sdr in volo, fotone virtuale, ecc.ecc
c	ed   s  = cm hadronico
c
c	v4   = coordinate vettore nel sdr s' in volo (input)
c	v4   = coordinate vettore nel sdr s  del lab (output)
c	w4   = coordinate vettore nel sdr s  del lab
c	pprime = vettore che rappresenta l'asse z
c		 del sdr s' rispetto ad s
c
	dimension ir(3)
	real a(3,3),ainv(3,3)
	real v4(4),w4(4)
	real cosdir_xprime(3),cosdir_yprime(3),cosdir_zprime(3)
	real pprime(3)
	real assexprime(3),asseyprime(3),assezprime(3)
	real assez(3)

c	assez sono i coseni direttori dell' asse z del sdr lab
c
	data assez /0., 0., 1./
c
c	assez' = direzione fotone virtuale rispetto al lab
	do i=1,3
	 assezprime(i)=pprime(i)
	enddo
c
c       il piano z'z coincide con il piano adronico....
c	assey' = perpendicolare al piano z'z
c	call cross(assezprime,assez,asseyprime)
	call cross(assez,assezprime,asseyprime)
c
c	assex' = perpendicolare ad assey' e assez'
	call cross(asseyprime,assezprime,assexprime)
c
c	normalizzo ad 1 per avere i coseni direttori
	do i=1,3
	 if(assexprime(i).ne.0) then
	  cosdir_xprime(i)=assexprime(i)/vmod(assexprime,3)
	 else
	  cosdir_xprime(i)=0
	 endif
	 if(asseyprime(i).ne.0) then
	  cosdir_yprime(i)=asseyprime(i)/vmod(asseyprime,3)
	 else
	  cosdir_yprime(i)=0
	 endif
	 if(assezprime(i).ne.0) then
	  cosdir_zprime(i)=assezprime(i)/vmod(assezprime,3)
	 else
	  cosdir_zprime(i)=0
	 endif
	enddo
c
c	riempio la matrice...
	do j=1,3
	    a(1,j)=cosdir_xprime(j)
	    a(2,j)=cosdir_yprime(j)
	    a(3,j)=cosdir_zprime(j)
	enddo
c
c	i calcoli che seguono servono per invertire una matrice 3x3
c	e ruotare un vettore dal sdr s' (v4) al sdr s (w4)
c
c !!!!	le matrici in fortran sono riempite nell'ordine:
c	a11,a21,a31,.....an1, a12,a22,a32,.....an2,.....
c
c	una matrice di nxm elementi (n=righe, m=colonne)
c	e' riempita per colonne, cioe' cosi':
c	1->n		prima colonna
c	n+1->2n		seconda colonna
c	2n+1->3n	terza colonna
c	...............
c	n(m-1)->mn	ultima colonna
c
c	la matrice di rotazione di un vettore da un sdr s ad un sdr s'
c	e' la seguente:
c
c			|				|
c			|	u	v	w	|
c			|	 x'	 x'	 x'	|
c			|				|
c	a	=	|	u	v	w	|
c			|	 y'	 y'	 y'	|
c			|				|
c			|	u	v	w	|
c			|	 z'	 z'	 z'	|
c
c	dove u    e' il primo   coseno direttore dell'asse x' nel sdr s
c	      x'
c
c	dove v    e' il secondo coseno direttore dell'asse x' nel sdr s
c	      x'
c
c	dove w    e' il terzo   coseno direttore dell'asse x' nel sdr s
c	      x'
c
c	dove u    e' il primo   coseno direttore dell'asse y' nel sdr s
c	      y'
c
c	dove v    e' il secondo coseno direttore dell'asse y' nel sdr s
c	      y'
c
c	ecc. ecc. ......
c
c	in questoo caso: s=lab,  s'=sdr fotone virtuale
c	oppure         : s=lab,  s'=sdr del c.m.
c	z' = direzione del fotone virtuale, oppure del vettore beta del cm
c	y' = perpendicolare al piano individuato dagli assi z e z' con
c	verso dato da : (z')x(z)   (x = prodotto vettore)
c	x' = perpendicolare al piano y'-z' con verso dato da : (y')x(z')
c
c	la matrice a, applicata al vettore di coordinate x,y,z nel sdr s
c	lo esprime secondo le coordinate x',y',z' nel sdr s'
c	a noi serve la trasformazione inversa -> devo invertire la matrice a
c	infatti: a = (s->s')  ==>  a^-1 = (s'->s)
c
c
	call rinv(3,a,3,ir,ifail)

c	se la matrice inversa non e' definita univocamente
c	come nei casi di particolari valori degli angoli => 99
	if (ifail .ne. 0) goto 99

	do i=1,3
	 do j=1,3
	  ainv(i,j)=a(i,j)
	 enddo
	enddo
c
c	ottenuta la matrice della trasformazione inversa s' -> s
c	la applico al vettore per trasformare le sue componenti
c	note nel sistema s' (v4) in quelle nel sistema s (w4)
c
	call vmatr(v4,ainv,w4,3,3)
	w4(4)=v4(4)

	do i=1,4
	   v4(i)=w4(i)
	enddo

	return
c

99	continue
c
c	casi particolari:
c
c
c	l'asse z' coincide con  +x  =>  a(3,1)= +1
c
	if(a(3,1).eq.1)  then
	 w4(1)= v4(3)
	 w4(2)= v4(2)
	 w4(3)=-v4(1)
	 w4(4)= v4(4)
	endif
c
c
c	l'asse z' coincide con  -x  =>  a(3,1)= -1
c
	if(a(3,1).eq.-1) then
	 w4(1)=-v4(3)
	 w4(2)= v4(2)
	 w4(3)= v4(1)
	 w4(4)= v4(4)
	endif
c
c
c	l'asse z' coincide con  +y  =>  a(3,2)= +1
c
	if(a(3,2).eq.1)  then
	 w4(1)= v4(1)
	 w4(2)= v4(3)
	 w4(3)=-v4(2)
	 w4(4)= v4(4)
	endif
c
c
c	l'asse z' coincide con  -y  =>  a(3,2)= -1
c
	if(a(3,2).eq.-1) then
	 w4(1)= v4(1)
	 w4(2)=-v4(3)
	 w4(3)= v4(2)
	 w4(4)= v4(4)
	endif
c
c
c	l'asse z' coincide con  +z  =>  a(3,3)= +1
c
	if(a(3,3).eq.1)  then
	 w4(1)= v4(1)
	 w4(2)= v4(2)
	 w4(3)= v4(3)
	 w4(4)= v4(4)
	endif
c
c
c	l'asse z' coincide con  -z  =>  a(3,3)= -1
c
	if(a(3,3).eq.-1) then
	 w4(1)=-v4(1)
	 w4(2)= v4(2)
	 w4(3)=-v4(3)
	 w4(4)= v4(4)
	endif

	do i=1,4
	   v4(i)=w4(i)
	enddo

	return
	end
