

C*******************************************************************************
	SUBROUTINE ROTATION
     @  (thx,phx,v,vout)
C*******************************************************************************
	implicit none
	REAL V(3),VPRIME(3),DROT(3,3)
	REAL V1,V2,V3,VPRIME1,VPRIME2,VPRIME3
	REAL THETA_E,PHI_E1
	real vout(3),thx,phx
	integer l,k

C**************************************************************
C	Evaluation of the components of the rotation matrix
C**************************************************************

	DROT(1,1) = COS(THX)*COS(PHX)
	DROT(1,2) = -COS(THX)*SIN(PHX)
	DROT(1,3) =   SIN(THX)

	DROT(2,1) = SIN(PHX)
	DROT(2,2) = COS(PHX)
	DROT(2,3) = 0.

	DROT(3,1) = -SIN(THX)*COS(PHX)
	DROT(3,2) = -SIN(THX)*SIN(PHX)
	DROT(3,3) = COS(THX)




	DO L = 1,3
	 DO K = 1,3
	IF(ABS(DROT(L,K)-1.).LT.1.E-6) DROT(L,K) = 1.
	IF(ABS(DROT(L,K)+1.).LT.1.E-6) DROT(L,K) = -1.
	IF(ABS(DROT(L,K)).LT.1.E-4) DROT(L,K) = 0.
	 ENDDO
	ENDDO


C**************************************************************
C	Evaluation of the components of the rotated vector
C**************************************************************

	DO L = 1,3
	VPRIME(L) = 0.
	 DO K = 1,3
	 VPRIME(L) = VPRIME(L) + V(K)*DROT(L,K)
	 ENDDO
	ENDDO

C--------------------------------------------------------
C	Fills the output values with the rotated vector
C--------------------------------------------------------

	do l=1,3
	   vout(l)=VPRIME(l)
	enddo

	RETURN

	END

	SUBROUTINE prod_vec(va,vb,thout,phout)
	implicit none 
	real va(3),vb(3),vout(3)
	real mod_vout,thout,phout
	
	mod_vout=sqrt((va(2)*vb(3)-va(3)*vb(2))**2
     $               +(va(3)*vb(1)-va(1)*vb(3))**2
     %               +(va(1)*vb(2)-va(2)*vb(1))**2)

	vout(1)=(va(2)*vb(3)-va(3)*vb(2))/mod_vout
	vout(2)=(va(3)*vb(1)-va(1)*vb(3))/mod_vout
	vout(3)=(va(1)*vb(2)-va(2)*vb(1))/mod_vout

	call angles(0,vout(1),vout(2),vout(3),thout,phout)
	RETURN
	END

c*******************************************************************************
	subroutine angles(iflag,v1,v2,v3,theta,phi)
c**********	
c+ INPUT	
c++	iflag = 0 calculating th and phi in degrees
c++	      = 1 calculating th and phi in rad
c++	v1,v2,v3 = input vector components
c+ OUTPUT
c++     theta, phi = calcalated angles      
c*******************************************************************************

	implicit none
	integer iflag
	real v1,v2,v3,v_m,theta,phi,php1

	v_m = sqrt(v1**2+v2**2+v3**2)

c-----------------------------------------
c	angles in the new frame (degrees)
c-----------------------------------------

	if(iflag.eq.0) then

	  if(v_m.gt.0.) then

	    theta = acos((v3/v_m))*180./3.14159265

	    if (v1.eq.0..and.v2.eq.0.) phi = 0

	    if (v1.gt.0..and.v2.ge.0.) then
 	    php1 = v2/v1
	    phi = atan(php1)*180./3.14159265
	    endif

	    if (v1.eq.0..and.v2.gt.0.) phi = 90.

	    if (v1.lt.0..and.v2.ge.0.) then
	    php1 = v2/v1
	    phi = atan(php1)*180./3.14159265 + 180.
	    endif
	
	    if (v1.lt.0..and.v2.lt.0.) then
	    php1 = v2/v1
	    phi = atan(php1)*180./3.14159265 + 180.
	    endif

	    if (v1.eq.0..and.v2.lt.0.) phi = 270.

	    if (v1.gt.0..and.v2.lt.0.) then
	    php1 = v2/v1
	    phi = atan(php1)*180./3.14159265 + 360.
	    endif

	  else

	    theta = 0.
	    phi = 0.

	  endif


	endif

c-----------------------------------------
c	angles in the new frame (radians)
c-----------------------------------------

	if(iflag.eq.1) then
	

	  if(v_m.gt.0.) then

	    theta = acos(v3/v_m)

	    if (v1.eq.0..and.v2.eq.0.) phi = 0.

	    if (v1.gt.0..and.v2.ge.0.) then
	    php1 = v2/v1
	    phi = atan(php1)
	    endif

	    if (v1.eq.0..and.v2.gt.0.) phi = 1.5708

	    if (v1.lt.0..and.v2.ge.0.) then
	    php1 = v2/v1
 	    phi = atan(php1) + 3.14159265
	    endif
	
	    if (v1.lt.0..and.v2.lt.0.) then
	    php1 = v2/v1
	    phi = atan(php1) + 3.14159265
	    endif

	    if (v1.eq.0..and.v2.lt.0.) phi = 4.7124

	    if (v1.gt.0..and.v2.lt.0.) then
	    php1 = v2/v1
	    phi = atan(php1) + 6.2831852
	    endif

	  else

	    theta = 0.
	    phi = 0.

	  endif


	endif
c-------------------------------------------------------------------------------


	return

	end



c------------------------------------------------------------------------------
