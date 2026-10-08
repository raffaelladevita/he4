
      subroutine gen2n1bod
      implicit none
      integer np,kgenev
      real    ecm,amass,pcm,wt
      common /genin /np,ecm,amass(18),kgenev
      common /genout/pcm(5,18),wt
      real    r,m1,m2,m3,theta,phi,ac,rran,q,theta2,phi2
	real pigr,twopigr,deltang
        data pigr/3.141592/
	data twopigr/6.283184/
	data deltang/0.157080/

      wt=1.0
      m1=amass(1)
      m2=amass(2)
c      if (m1.gt.m2) then
c         m2=amass(1)
c         m1=amass(2)
c      endif
      ac=rran()*2-1.
      theta =acos(ac) 
      phi = rran()*twopigr
      q=((ecm**2-(m1+m2)**2)*(ecm**2-(m2-m1)**2))**0.5/2/ecm
c      print *, theta,phi,ecm,q

      pcm(1,1)=q*sin(theta)*cos(phi)
      pcm(2,1)=q*sin(theta)*sin(phi)
      pcm(3,1)=q*cos(theta)
      pcm(5,1)=sqrt(pcm(1,1)**2+pcm(2,1)**2+pcm(3,1)**2)
      pcm(4,1)=sqrt(pcm(5,1)**2+amass(1)**2)
      theta2=pigr-theta
      phi2=pigr+phi
      pcm(1,2)=q*sin(theta2)*cos(phi2)
      pcm(2,2)=q*sin(theta2)*sin(phi2)
      pcm(3,2)=q*cos(theta2)
      pcm(5,2)=sqrt(pcm(1,2)**2+pcm(2,2)**2+pcm(3,2)**2)
      pcm(4,2)=sqrt(pcm(5,2)**2+amass(2)**2)
c      print *, ecm,pcm(4,1)+pcm(4,2)
c      print *, m1,m2
      return
      end

c======================================================================
c     gen3n1bod 
c     input/output through commons (see genbod)
c     np - must be 3 - the num of particles -- input
c     tecm = w -- input 
c     amass(18)(1,2,3) masses -- input
c     pcm(1,i),pcm(2,i),pcm(3,i)-px,py,pz momenta for i-th part.--output
c     orther parameters in /genin / and /genout/ are not used
c======================================================================

      subroutine gen3n1bod
      implicit none
      integer np,kgenev
      real    tecm,amass,pcm,wt
      common /genin /np,tecm,amass(18),kgenev
      common /genout/pcm(5,18),wt

      integer nc
      real    r,m1,m2,m3
      real    pii,w,s12,s23,s12min,s12max,s23min,s23max
      real    cteta,fi,psi,teta
      real    q1(0:3),q2(0:3),q3(0:3)

      pii = acos(-1.)

      if(np.ne.3) then
        print *,' sub gen3n1bod: np.ne.3=',np
        stop
      endif
      w=tecm      
      m1=amass(1)
      m2=amass(2)
      m3=amass(3)
      if(w.lt.m1+m2+m3) then
        print *,' sub gen3n1bod: w is too small=',w
        stop
      endif

 1001 continue
      s12min=(m1+m2)**2
      s12max=(w-m3)**2
      s23min=(m2+m3)**2
      s23max=(w-m1)**2
      call ranlux(r,1)
      s12=s12min + r*(s12max-s12min)
      call ranlux(r,1)
      s23=s23min + r*(s23max-s23min)
      call ranlux(r,1)
      cteta=-1. + r*(2.)
      teta=acos(cteta)
      call ranlux(r,1)
      fi = 0.  + r*(2.*pii-0.)
      call ranlux(r,1)
      psi= 0. +  r*(2.*pii-0.)

      call kkk3genb3n1(q1,q2,q3,m1,m2,m3,w,s12,s23,teta,fi,psi, nc)
      if(nc.ne.0) goto 1001

      pcm(1,1)=q1(1)
      pcm(2,1)=q1(2)
      pcm(3,1)=q1(3)
      pcm(5,1)=sqrt(pcm(1,1)**2+pcm(2,1)**2+pcm(3,1)**2)
      pcm(4,1)=sqrt(pcm(5,1)**2+m1**2)

      pcm(1,2)=q2(1)
      pcm(2,2)=q2(2)
      pcm(3,2)=q2(3)
      pcm(5,2)=sqrt(pcm(1,2)**2+pcm(2,2)**2+pcm(3,2)**2)
      pcm(4,2)=sqrt(pcm(5,2)**2+m2**2)


      pcm(1,3)=q3(1)
      pcm(2,3)=q3(2)
      pcm(3,3)=q3(3)
      pcm(5,3)=sqrt(pcm(1,3)**2+pcm(2,3)**2+pcm(3,3)**2)
      pcm(4,3)=sqrt(pcm(5,3)**2+m3**2)
  
      return
      end


c------------------------------kkk3genb--------------------------------
c  for:
c           w -> 1(pi-) + 2(pi+) + 3(delta) (in cms)
c  calculates:
c    4-momenta of the particles
c  starting from: w,s12,s23,teta1,fi1,fi23
c  ww      -                                               -- input
c  p1(0:3),p2(0:3),p3(0:3) - 4-momenta                     -- output
c  m1,m2,m3 - masses
c  nc - exit code  0 - ok                                  -- output
c                  1 - out of allowed kin. region
c  look for the comments in kinema3
c----------------------------------------------------------------------

      subroutine kkk3genb3n1(p1,p2,p3,m1,m2,m3,
     &                       ww,s12,s23,teta1,fi1,fi23, nc)
      implicit none
      real*4 ww,s12,s23
      real*4 p1  (0:3),p2  (0:3),p3  (0:3),p1_m,  p2_m,  p3_m
      real*4 p1_1(0:3),p2_1(0:3),p3_1(0:3),p1_1_m,p2_1_m,p3_1_m
      real*4 th1,ph1, th2,ph2, th3,ph3
      real*4 m1,m2,m3,teta1,fi1,fi23,m1_2,m2_2,m3_2
      integer*4 nc
      real*4 zero,quasizero,pii, s,s31, g, a,b
      real*4 g_byckling4,x,y,z,u,v,w
      

c ----- g-function of byckling, byckling, p.89 -----
      g_byckling4(x,y,z,u,v,w) =
     &   (x**2)*y+x*(y**2) + (z**2)*u+z*(u**2) + (v**2)*w+v*(w**2)
     &   + x*z*w + x*u*v + y*z*v + y*u*w
     &   - x*y*(z+u+v+w) - z*u*(x+y+v+w) - v*w*(x+y+z+u)


      zero=0.0d+0
      pii = acos(-1.d+0)
      nc=0
      s = ww**2
      s31 = s - s12 - s23 + m1**2 + m2**2 + m3**2
      m1_2=m1**2
      m2_2=m2**2
      m3_2=m3**2

c ----- -----
      g = g_byckling4(s12,s23,s,m2_2,m1_2,m3_2)
      if(g.gt.zero) then
        nc=1
        return
      endif

c----------------------------------------
c     moduli of the momenta and energies
c----------------------------------------
      p1(0) = (s+m1_2-s23)/2./ww
      p2(0) = (s+m2_2-s31)/2./ww
      p3(0) = (s+m3_2-s12)/2./ww
      p1_m = sqrt(p1(0)**2-m1_2)
      p2_m = sqrt(p2(0)**2-m2_2)
      p3_m = sqrt(p3(0)**2-m3_2)

c------------------------------------------------------------------
c     angles and three momenta of the particles in the frame 
c     with the z-axes along the p1, x-axes perp. to the z situated
c     on the a-b-1-plane
c------------------------------------------------------------------
      th1 = 0.d+0
      ph1 = 0.d+0
      a = ((m1_2+m2_2 + 2.*p1(0)*p2(0) - s12)/2./p1_m/p2_m)
      if(abs(a).gt.1.)then
        nc=1
        return
      endif
      th2 = acos(a)
      ph2 = fi23
      a = ((m1_2+m3_2 + 2.*p1(0)*p3(0) - s31)/2./p1_m/p3_m)
      if(abs(a).gt.1.)then
        nc=1
        return
      endif
      th3 = acos(a)
      ph3 = fi23 + pii
      if(ph3 .gt. 2.*pii) ph3 = ph3 - 2.*pii

c------------------------------------
c     three momenta of the particles 
c------------------------------------

      p1(1) = 0.0
      p1(2) = 0.0
      p1(3) = p1_m
      p2(1) = p2_m*sin(th2)*cos(ph2)
      p2(2) = p2_m*sin(th2)*sin(ph2)
      p2(3) = p2_m*cos(th2)
      p3(1) = p3_m*sin(th3)*cos(ph3)
      p3(2) = p3_m*sin(th3)*sin(ph3)
      p3(3) = p3_m*cos(th3)

c----------------------------------------------------------------------
c     rotation to the hadronic cms 
c     (z axis along q vector but x axis is still on the hadronic plane)
c----------------------------------------------------------------------

      p1_1(1) =  cos(teta1)*p1(1) + sin(teta1)*p1(3)
      p1_1(2) =  p1(2)
      p1_1(3) = (-sin(teta1))*p1(1) + cos(teta1)*p1(3)
      p1_1(0) = p1(0)
      p1_1_m  = p1_m
      p2_1(1) =  cos(teta1)*p2(1) + sin(teta1)*p2(3)
      p2_1(2) =  p2(2)
      p2_1(3) = (-sin(teta1))*p2(1) + cos(teta1)*p2(3)
      p2_1(0) = p2(0)
      p2_1_m  = p2_m
      p3_1(1) = - p1_1(1) - p2_1(1)
      p3_1(2) = - p1_1(2) - p2_1(2)
      p3_1(3) = - p1_1(3) - p2_1(3)
      p3_1(0) = p3(0)
      p3_1_m  = p3_m

c----------------------------------------------------------
c     rotation to the lab plane (z axis is still along the 
c     virtual photon but x axis on the ee' plane)
c----------------------------------------------------------

      p1(1) = cos(fi1)*p1_1(1) - sin(fi1)*p1_1(2)
      p1(2) = sin(fi1)*p1_1(1) + cos(fi1)*p1_1(2)
      p1(3) = p1_1(3)
      p1(0) = p1_1(0)
      p1_m  = p1_1_m
      p2(1) = cos(fi1)*p2_1(1) - sin(fi1)*p2_1(2)
      p2(2) = sin(fi1)*p2_1(1) + cos(fi1)*p2_1(2)
      p2(3) = p2_1(3)
      p2(0) = p2_1(0)
      p2_m  = p2_1_m 
      p3(1) = - p1(1) - p2(1)
      p3(2) = - p1(2) - p2(2)
      p3(3) = - p1(3) - p2(3)
      p3(0) = p3_1(0)
      p3_m  = p3_1_m 

c--------------------------------------------------
c     check for kinematics 
c--------------------------------------------------

c ----- 4-momentum conservation:  -----
      quasizero=(ww)*1.d-5
      if(
     &   abs(p1(1)+p2(1)+p3(1)).gt.quasizero .or.
     &   abs(p1(2)+p2(2)+p3(2)).gt.quasizero .or.
     &   abs(p1(3)+p2(3)+p3(3)).gt.quasizero .or.
     &   abs(p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2-m1**2).gt.quasizero.or.
     &   abs(p2(0)**2-p2(1)**2-p2(2)**2-p2(3)**2-m2**2).gt.quasizero.or.
     &   abs(p3(0)**2-p3(1)**2-p3(2)**2-p3(3)**2-m3**2).gt.quasizero
     &  ) then
        print *,' sub. kkk3: kinematics is bad'
        print *,' p1(1)+p2(1)+p3(1)=',p1(1)+p2(1)+p3(1)
        print *,' p1(2)+p2(2)+p3(2)=',p1(2)+p2(2)+p3(2)
        print *,' p1(3)+p2(3)+p3(3)=',p1(3)+p2(3)+p3(3)
        print *,' m1,3v_p1**2=',m1,
     &        abs(p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2-m1**2)
        print *,' m2,3v_p2**2=',m2,
     &        abs(p2(0)**2-p2(1)**2-p2(2)**2-p2(3)**2-m2**2)
        print *,' m3,3v_p3**2=',m3,
     &        abs(p3(0)**2-p3(1)**2-p3(2)**2-p3(3)**2-m3**2)
        stop
      endif
      if(abs(p1(0)+p2(0)+p3(0)-ww).gt.quasizero ) then
        print *,' sub. kkk3: kinematics is bad'
        print *,' ww=',ww
        print *,' p1(0)+p2(0)+p3(0)=',p1(0)+p2(0)+p3(0)
        stop
      endif

c ----- check s12,s23 -----
      quasizero = (ww)*1.d-5
      a = (p1(0)+p2(0))**2
     &  - (p1(1)+p2(1))**2 - (p1(2)+p2(2))**2 - (p1(3)+p2(3))**2
      b = (p2(0)+p3(0))**2
     &  - (p2(1)+p3(1))**2 - (p2(2)+p3(2))**2 - (p2(3)+p3(3))**2
      if(abs(a-s12).gt.quasizero .or. abs(b-s23).gt.quasizero) then
        print *,' ...sub. kkk3 kinematics is bad'
        print *,' s12_calc,s12=',a,s12
        print *,' s23_calc,s23=',b,s23
        stop
      endif

      return
      end
      

