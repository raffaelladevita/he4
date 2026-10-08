      function adam_cros(eg,t,xm12,tgj,phi)

c     adam_cros = d\sigma/dM_kk dcos\theta d\phi) 
c     input: 
c     eg = lab fram photon energy
c     t = NEGATIVE of momentum transfer (that is t is a positive numebr e.g. 0.4)
c     xm12 = KK invariant mass
c     tgj = cos(theta) [-1:1] polar anganle in helicity frame 
c     phi [0:2*pi] azymuthal 
c     helicity frame: rest frame of KK with z axis along the negative of the recoil proton
c     y axis perendiculat to production plane. 


      implicit double precision (a-h,o-z)
      double complex Rfactor
      double complex ppn,ppf,p0n,p0f
      double complex apn,apf,amn,amf
      double complex Y11
      double complex xi
c	print*, xm12

	adam_cros=0
      pi = 4.d0*datan(1.d0) 
      xi = (0.d0,1.d0)

      xmn = 0.9383     
      xmk = 0.49367
      s = xmn**2 + 2.*xmn*eg
      xmf = 1.019455
      xmf = 1.08
      gf = 0.00426

      ee1 = (s - xmn**2 - t)/(2.*xm12)
      ee2 = (s - xmn**2 - xm12**2)/(2.*xm12)
      qqq = (xm12**2 + t)/(2.*xm12)
	if ((ee2**2 - xmn**2).lt.0) return
      ppp = dsqrt(ee2**2 - xmn**2)
      zz2 = (ee1**2 - ee2**2 - qqq**2)/(2*qqq*ppp)

      cc1 = (xmn**2-ee1*ee2+t/2.)/(ppp*dsqrt(ee1**2-xmn**2))
	if (((1.+cc1)/2.).lt.0) return
      tc = dsqrt((1.+cc1)/2.)
	if (((1.-cc1)/2.).lt.0) return
  
      ts = dsqrt((1.-cc1)/2.)
 
      scale = 1./15.
      Rfactor = -2.*scale*dsqrt(xm12**2/4.-xmk**2)
     $ /(xmf**2 - xm12**2 - xi*xmf*gf)

c      Rfactor = 1.

      ppn = 2.*Rfactor*xm12*ts
      ppf = 2.*Rfactor*xm12*tc

      p0n = -Rfactor*ts*dsqrt(2.*t)
      p0f = -Rfactor*tc*dsqrt(2.*t)

      Y11 = -dsqrt(1.-tgj**2)*(dcos(phi)+xi*dsin(phi))/dsqrt(2.d0)
      Y10 = tgj

      apn = ppn*Y11 + p0n*Y10
      apf = ppf*Y11 + p0f*Y10 

      amn =  ppn*(-dconjg(Y11)) - p0n*Y10
      amf = -ppf*(-dconjg(Y11)) + p0f*Y10

      adam_cros = dble( apn*dconjg(apn) + apf*dconjg(apf)
     $     + amn*dconjg(amn) + amf*dconjg(amf) )
	if ((xm12**2/4.-xmk**2).lt.0) return
      adam_cros = dsqrt(xm12**2/4.-xmk**2)*adam_cros

c	print *, ee1, ee2, qqq , ppp, zz2
c	print *, 'cros=',cros
      return
      end 






      function adam_cros_sal(eg,t,xm12,tgj,phi)


c     adam_cros = d\sigma/dM_kk dcos\theta d\phi) 
c     input: 
c     eg = lab fram photon energy
c     t = NEGATIVE of momentum transfer (that is t is a positive numebr e.g. 0.4)
c     xm12 = KK invariant mass
c     tgj = cos(theta) [-1:1] polar anganle in helicity frame 
c     phi [0:2*pi] azymuthal 
c     helicity frame: rest frame of KK with z axis along the negative of the recoil proton
c     y axis perendiculat to production plane. 



      implicit double precision (a-h,o-z)
      double complex Rfactor
      double complex ppn,ppf,p0n,p0f
      double complex apn,apf,amn,amf
      double complex sn 
      double complex Y11
      double complex xi


	adam_cros_sal=0
      pi = 4.d0*datan(1.d0) 
      xi = (0.d0,1.d0)

      xmn = 0.9383     
      xmk = 0.49367
      s = xmn**2 + 2.*xmn*eg
      xmf = 1.019455
      gf = 0.00426
c      gf = 0.01
      ee1 = (s - xmn**2 - t)/(2.*xm12)
      ee2 = (s - xmn**2 - xm12**2)/(2.*xm12)
      qqq = (xm12**2 + t)/(2.*xm12)
	if ((ee2**2 - xmn**2).lt.0) return
      ppp = dsqrt(ee2**2 - xmn**2)
      zz2 = (ee1**2 - ee2**2 - qqq**2)/(2*qqq*ppp)

      cc1 = (xmn**2-ee1*ee2+t/2.)/(ppp*dsqrt(ee1**2-xmn**2))
	if (dabs(cc1).gt.1.) return 
	if (((1.+cc1)/2.).lt.0) return
      tc = dsqrt((1.+cc1)/2.)
	if (((1.-cc1)/2.).lt.0) return
      ts = dsqrt((1.-cc1)/2.)
 
      scale = 1./15.
      Rfactor = -2.*scale*dsqrt(xm12**2/4.-xmk**2)
     $ /(xmf**2 - xm12**2 - xi*xmf*gf)

c      Rfactor = 1.

      ppn = 4.*Rfactor*xm12*ts
      ppf = 4.*Rfactor*xm12*tc

      p0n = 0.25*Rfactor*ts*dsqrt(2.*t)
      p0f = 0.25*Rfactor*tc*dsqrt(2.*t)

      sn = 2.0*-Rfactor*dsqrt(2.*t)*ts*xm12
c      sn = 0.

      Y11 = -dsqrt(1.-tgj**2)*(dcos(phi)+xi*dsin(phi))/dsqrt(2.d0)
      Y10 = tgj
      Y00 = 1. 

      apn = ppn*Y11 + p0n*Y10
      apf = ppf*Y11 + p0f*Y10 + sn*Y00

      amn =  ppn*(-dconjg(Y11)) - p0n*Y10
      amf = -ppf*(-dconjg(Y11)) + p0f*Y10 + sn*Y00

  
      adam_cros_sal = dble( apn*dconjg(apn) + apf*dconjg(apf)
     $     + amn*dconjg(amn) + amf*dconjg(amf) )

	if ((xm12**2/4.-xmk**2).lt.0) return
      adam_cros_sal = dsqrt(xm12**2/4.-xmk**2)*adam_cros_sal



      return
      end 
