	subroutine bdisnew4he3(eic1,ep1,theta1,aa1,zz1,ann1,esep1,pmax1,
     &	innt1,innp1,f01,bigB1,ag1,bg1,alpha1,sigdeep,w1a,w2a)
        
        Implicit none
c	implicit real*8 (a-z)
	real*8 wp,ww1,ww2, w11, w22, nu,p,xp
	real*8 eic1,ep1,theta1,aa1,zz1,ann1,esep1,pmax1
	real*8 w11tmp,w22tmp
        real*8  innt1,innp1,ag1,bg1,cg1,dg1,sigdeep
c	modified version of BWF's routine bdis
c	ag,bg,cg,dg are coeff of a 2 gaussian fit to the f(y)
c	which was converted to n(k)
	real*8 rk,rho
	real*8 rq2		!to pass to w1w2

        integer inp
        real*8 f01,bigB1,alpha1,thec
        real*8 eic,ep,aa,zz,ann,esep,pmax
        real*8 innt,innp,ag,bg,f0,BigB,alpha
        real*8 agd,bgd,cgd,dgd
        real*8 theta,pi,s2,anu,q2,amp,x,amp2
        real*8 amn,ama_1,q3
        real*8 w1ap,w2ap,w1an,w2an
        real*8 w1ap1,w2ap1,w1an1,w2an1
        real*8 w1,w2,w1a,w2a,sigm,tt2
        real*8 arg1,arg21,arg22,arg2
        real*8 weight1,weight2
        real*8 wp_max,ama,du,akf,dp,u
        real*8 ei,anup,wpmin2,radical,fackf
        real*8 ip,iu
	real*8 r ! dummy variable for f1f2in09


c	common/interp/npts,nterms,rk(200),rho(200)
	
c	eic = real(eic1)
c	ep = real(ep1)
c	theta = real(theta1)
c	aa = real(aa1)
c	zz = real(zz1)
c	ann = real(ann1)
c	esep = real(esep1)
c	pmax = real(pmax1)
c	innt = real(innt1)
c	innp = real(innp1)
c	ag = real(ag1)
c	bg = real(bg1)
c	f0 = real (f01)
c	BigB = real (bigB1)
c	alpha=real(alpha1)


	eic = eic1
	ep = ep1
	theta = theta1
	aa = aa1
	zz = zz1
	ann = ann1
	esep = esep1
	pmax = pmax1
	innt = innt1
	innp = innp1
	ag = ag1
	bg = bg1
	f0 = f01
	BigB = bigB1
	alpha=alpha1
	
	f0=f0*1000
	bigB=bigB*1000
	alpha=alpha/1000
	ag=ag*1000
	bg=bg*1000

	
        if (aa .eq. 3)then       !he3
	  agd =    0.39167E+01
	  bgd =    0.78468E+02
	  cgd =    0.17924E+00
	  dgd =    0.66511E+01
        elseif (aa .eq. 4)then	!he4
	  agd =   0.35160E+01
	  bgd =   0.60049E+02
	  cgd =   0.38390E+00
	  dgd =   0.13421E+02
c        elseif ((aa .eq. 12).or.(aa.eq.9).or. 
c     >          (aa.eq.10).or.(aa.eq.11)) then      !carbon
	elseif ( (aa.gt.5.5) .and. (aa.lt.12.5)) then
	  agd =   0.28757E+01
	  bgd =   0.41922E+02
	  cgd =   0.33801E+00
	  dgd =   0.13824E+02
        elseif (aa .eq. 27)then      !alum
	  agd =   0.25660E+01
	  bgd =   0.36962E+02
	  cgd =   0.41675E+00
	  dgd =   0.13772E+02
        elseif ((aa.ge.39.0) .and. (aa.le.80.0))then   !iron
	  agd =   0.28316E+01
	  bgd =    0.44624E+02
	  cgd =    0.37850E+00
	  dgd =    0.12003E+02
        elseif(aa .gt. 80.0)then	!gold, Jerry, Gold!
	  agd =   0.24947E+01
	  bgd =   0.30614E+02
	  cgd =   0.42398E+00
	  dgd =   0.12272E+02
c	else
c	   write(6,*) 'cannot determine nucleus - crap!'
c	   stop
        endif

c       temp values for check.

c	ag =   0.75287E+01
c	bg =   0.18856E+03
c	cg =   0.22657E+00
c	dg =   0.14649E+02
c	eic incident energy in gev
c	ep final energy
c	thec scattering angle in radians
c	aa atomic number
c	zz number of protons
c	ann number of neutrons
c	esep spearation energy
c	pmax max p in integration of p > kf
c	innt integration steps over theta
c	integration steps over p
c       
c       
c	write(6,*)'E,EP,THETA, A, ZZ,ANN,ESP,PMAX,INNT,INNp,ag,bg,bigB, f0, alpha'
c	write(6,*)eic,ep,theta,aa,zz,ann,esep,pmax,innt,innp,ag,bg,bigB,
c	1 f0, alpha

	pi = 3.14159265
	thec=theta
	theta=theta*pi/180
	s2 = sin(theta/2.)**2
	anu = eic - ep
	q2 = 4.*eic*ep*s2
	
c	write (6,*) 'in bdis, got q2 of ', q2
	rq2=q2		! to pass to w1w2
c	write(6,*)eic,ep,thec,esep,q2, rq2,


	amp = 0.938273
	x = q2/2./amp/(eic-ep)
	amp2 = amp*amp
	amn = 0.939565 
	if(aa.eq.2) then
	  ama_1=0.93827
	elseif(aa .gt. 2 .and. aa .lt. 4)then !he3
	  ama_1 = 1.87609
	elseif (aa.eq.4) then
	  ama_1 = 2.8094
	elseif (aa.eq.6) then
	   ama_1=4.6733
	elseif (aa.eq.7) then
	   ama_1=5.6056
	elseif (aa.eq.9) then
	  ama_1 = 7.5402
	elseif (aa.eq.10) then
	  ama_1 = 8.4060
	elseif (aa.eq.11) then
	  ama_1 = 9.3341
	elseif (aa.eq.12) then
	  ama_1 = 10.25553
	elseif (aa.eq.27) then
	  ama_1 = 24.205
	elseif (aa.eq.40) then
	  ama_1 = 37.327
	elseif (aa.eq.48.and.zz.eq.20) then
	   ama_1 = 43.7378
	elseif (aa.eq.48 .and.zz.eq.22) then
	   ama_1 = 43.6581
	elseif(aa.eq.54) then
	   ama_1 = 49.3147
	elseif (aa.eq.56) then
	  ama_1 = 51.1743
	elseif(aa.eq.58) then
	   ama_1 = 53.0367
	elseif(aa.eq.64) then
	   ama_1=58.6188
	elseif (aa.eq.64) then
	  ama_1 = 57.6886
	elseif (aa.eq.108) then
	  ama_1 = 99.5473
	elseif (aa.eq.119) then
	   ama_1 = 109.6479
	elseif (aa.eq.197) then
	  ama_1 = 182.5394
	elseif (aa.eq.232) then
	   ama_1 = 215.2104
	endif
	ama = ama_1 + amp - esep
c	write (6,*) 'ama_1 is ', ama_1
c
	q3 = sqrt(q2 +anu*anu)
	w1ap = 0.0
	w2ap = 0.0
	w1an = 0.0
	w2an = 0.0
	du = 2./innt
	akf = (1. - exp(-aa/8.))*.22 + 0.04
c
	dp = pmax/innp
c
c	calculate proton and neutron free structure functions as a check if
c	desired
c	inp = 1
c	w0 = sqrt(amp2 + 2.*amp*amu - q2)
c	call w1w2(q2,w0,w1,w2,inp)
c	w2p = w2
c	inp = 2
c	call w1w2(q2,w0,w1,w2,inp)
c	w2n = w2
	wp_max = 0.0
c	do smearing correction
	do ip = 1,innp
	  p = (ip - .5)*dp
	  w1ap1 = 0.0
	  w2ap1 = 0.0
	  w1an1 = 0.0
	  w2an1 = 0.0
	  
	  do 10 iu = 1,innt

	    u = (iu - .5)*du
	    u = u - 1.
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c	i removed this dependence on kf 12/4/87 just to see the
c	effect
c       if(p.le.akf)then
c       
	    ei = ama - sqrt(p*p+ama_1*ama_1)
c       
c       else
c       ei = amd - sqrt(p*p + amp**2)
c       endif
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    anup = (ei*anu - p*u*q3)/amp
c	note that w1 and w2 are zero if wp < M + mpi
c	
	    wpmin2 = (amp + .135)**2
cdg	    if (iu.eq.1) then
c       write (6,*)'This is ', this
c	      write (6,*) 'wpmin2', wpmin2, ' thingie ', (ei*ei - p*p +
c	1	2.*ei*anu - 2.*p*u*q3 - q2), ei, p, anu, u, q3, q2
c       write (6,*) 'Ei is,',ei
cdg	    endif
	    if((ei*ei - p*p + 2.*ei*anu - 2.*p*u*q3 - q2).le
	1     .wpmin2) goto 10
	    radical =(ei*ei - p*p + 2.*ei*anu - 2.*p*u*q3 - q2)
	    wp = sqrt(radical)
				!write (6,*) 'wp (nu) is ',wp
				!write (6,*) 'real nu is ', anu
	    inp = 1
				!write(6,*)'In if stmt, about to go to w1w2'
	    w11=0.
	    w22=0.
	    nu = (wp**2 - amp**2 + q2) / 2. / amp
c	    write(*,*) 'my nu is ', nu, wp**2, q2, rq2
cdg	    call w1w2(rq2,wp,w1,w2,inp)
c	    call F1F2IN06(1.D0,1.D0,rq2,wp**2,w11,w22)
c	    call F1F2IN09(1.D0,1.D0,rq2,wp**2,w11,w22,r)
	    call ineft(rq2,wp,w11,w22,dble(1.0))
	    if((w11.eq.0).or.(w22.eq.0)) goto 10
c	    w11=w11/amp
c	    w22=w22/nu

	    w1=w11
	    w2=w22

	    if (wp .gt. wp_max) wp_max = wp

	    arg1 = w1 + (1 -u*u)*p*p*w2/2./amp2
	    arg21 = (1. + p*u*q2/amp/anup/q3)**2*anup*anup/anu/anu
c       errorinBR	arg21 = (1. - p*u*q2/amp/anup/q3)*anup*anup/anu/anu
	    arg22 = (1.- u*u)*p*p*q2/2./amp2/q3/q3
	    arg2 = w2*(arg21 + arg22)
	    w1ap1 = w1ap1 + arg1*du
	    w2ap1 = w2ap1 + arg2*du
	    inp = 2
cdg	    call w1w2(rq2,wp,w1,w2,inp)	    
	    w11=0.
	    w22=0.
c	    call F1F2IN06(0.D0,1.D0,rq2,wp**2,w11,w22)
c	    call F1F2IN09(0.D0,1.D0,rq2,wp**2,w11,w22,r)
c	    call ineft(rq2,wp,w11tmp,w22tmp,dble(1.0))
	    call ineft(rq2,wp,w11,w22,dble(0.0))
c	    w11=2.0*w11-w11tmp
c	    w22=2.0*w22-w22tmp
	    if((w11.eq.0).or.(w22.eq.0)) goto 10
c	    w11=w11/amp
c	    w22=w22/nu

	    w1=w11
	    w2=w22

	    arg1 = w1 + (1 -u*u)*p*p*w2/2./amp2
c       errorinBR	arg21 = (1. - p*u*q2/amp/anup/q3)*anup*anup/anu/anu
	    arg21 = (1. + p*u*q2/amp/anup/q3)**2*anup*anup/anu/anu
	    arg22 = (1.- u*u)*p*p*q2/2./amp2/q3/q3
	    arg2 = w2*(arg21 + arg22)
	    w1an1 = w1an1 + arg1*du
	    w2an1 = w2an1 + arg2*du
c       
 10	  continue



	  weight1 = ((exp(-(ag*p)**2)*(f0-bigB)*alpha**2)*(((ag**2)
	1   *(alpha**2+p**2)+1)/(alpha**2+p**2)**2))


	  
	  if (aa.eq.2) then
	    if(p.lt.0) then
	      weight1=weight1-(bg*bigB*exp(-bg*abs(p))/(2*p))
	    else
	      weight1=weight1+(bg*bigB*exp(-bg*abs(p))/(2*p))
	    endif
	  else
	      weight1=weight1+(bg**2*bigB*exp(-(bg*p)**2))
	  endif


	  weight1=weight1*p*p*dp/pi!*2.48176

	  weight2 = (agd*bgd/pi*exp(-bgd*p**2) + cgd*dgd/pi*exp(-dgd*p
	1   **2) )*p*p*dp
	  
	  w1ap = w1ap + w1ap1*weight1
	  w2ap = w2ap + w2ap1*weight1
	  w1an = w1an + w1an1*weight1
	  w2an = w2an + w2an1*weight1
	enddo
	
	fackf = 2.*pi
c	
	w1a = fackf*(zz*w1ap + ann*w1an)
	w2a = fackf*(zz*w2ap + ann*w2an)

	sigm = cos(theta/2.)/(2.*137.*eic*s2)
c        sigm=0.389379*sigm**2  !(GeV**2 mb )
	sigm = (0.1973*sigm)**2
	tt2 = tan(theta/2.)**2
	sigdeep = 1.e7*sigm*(w2a + 2.*w1a*tt2)
c       write(6,*)
c	write(6,*)'End of Bdis, DIS, MOTT',sigdeep, sigm, tt2, w1a, w2a,w1ap,zz,ann
c	write(50,*)x,q2,wp_max

	return
	end
	
