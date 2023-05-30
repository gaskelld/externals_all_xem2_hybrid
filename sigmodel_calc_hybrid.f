	subroutine sigmodel_calc(e1pass,e2pass,thpass,zpass,apass,mpass,sig_dis_pass,sig_qe_pass,sig,xflag,factpass)
C       +______________________________________________________________________________
c	
C       Calculate cross section using Peter's F1F209.f routine
c       
c       ARGUMENTS:
c       
c       E1:		-	Incident energy in GeV.
c       E2:		- Scattered energy in GeV.
c       TH:		- Scattering angle in Degrees.
c       A:		- 'A' of nucleus.
c       Z:		- Number of protons in nucleus.
c       M_TGT:	- Mass of target nucleus in GeV/c2.
c       M_REC:	- Mass of recoiling nucleon in GeV/c2.
c       E_SEP:	- Separation energy for target nucleus in GeV/c2.
c       SIG  :	- Calculated cross section in nb/(MeV-ster).
C       ______________________________________________________________________________

        implicit none
	include 'math_physics.inc'
	include 'include.inc'

C       Declare arguments.

	real*8		e1,e2,th,a,z,m_tgt
	real*4          e1pass,e2pass,thpass,mpass
	real*4          sig,factpass,sig_dis_pass,sig_qe_pass
	integer         zpass,apass
	logical		modelhack

C       Declare locals.

	real*8	 	sig_qe,sig_dis,y,normfac,fact
        real*8		thr,cs,sn,tn,elastic_peak
	real*8          Q2,nu,WSQ, x
	real*8          F1,F2,W1,W2,sigmott,r
	real*8          W1p,W1n,W1D,W2p,W2n,W2D
	real*8          sig_smear,w1smear,w2smear,w1ineft,w2ineft
	real*8          inelastic_it,smear_cor
	real*8          f1mec,f2mec,W1mec,W2mec
	real*8          a1,b1,bigB,eps,f0,alpha1,pmax
	real*8          frac
	integer         xflag !flag for which xsec to calculate 1=both 2=QE only 3=DIS only
	logical         first


	save

        real*8 emc_func_xem
	external emc_func_xem

	real*8 emc_func_slac
	external emc_func_slac

	data first/.true./

	e1=dble(e1pass)
	e2=dble(e2pass)
	th=dble(thpass)
	a=dble(apass)
	z=dble(zpass)
	m_tgt=dble(mpass)


	sig =0.0
	sig_qe=0.0
	sig_dis=0.0

C       If at or beyond elastic peak, return ZERO!

	thr = th*d_r
	cs = cos(thr/2.)
	sn = sin(thr/2.)
	tn = tan(thr/2.)
	elastic_peak = e1/(1.+2.*e1*sn**2/m_tgt)
      	if (e2.ge.elastic_peak) then
       	   sig = 0.0
       	   return
       	endif
c	write(6,*) 'got to 1'

	Q2 = 4.*e1*e2*sn**2
	nu=e1-e2
	WSQ = -Q2 + m_p**2 + 2.0*m_p*nu 
        x = Q2/2/m_p/nu



	F1=0
	F2=0
	r=0
	if((xflag.eq.1).or.(xflag.eq.3)) then
c----------------------------------------------------------------
c       
c       do inelastic stuff
C Use old Bodek fit + SLAC EMC fit for now, b/c F1F2IN09 doesn't like large Q2,W2
	   W1=0.0
	   W2=0.0
	   if(wsq.gt.1.16) then
	      if(A.eq.2.0) then
		 call ineft(Q2,sqrt(wsq),W1D,W2D,dble(2.0))
		 W1ineft=2.0*W1D
		 W2ineft=2.0*W2D
	      else
		 call ineft(Q2,sqrt(wsq),W1p,W2p,dble(1.0))
c		 call ineft(Q2,sqrt(wsq),W1D,W2D,dble(2.0))
c		 W1n=2.0*W1D-W1p
c		 W2n=2.0*W2D-W2p
		 call ineft(Q2,sqrt(wsq),W1n,W2n,dble(0.0))
		 W1ineft=Z*W1p+(A-Z)*W1n ! + W1mec
		 W2ineft=Z*W2p+(A-Z)*W2n ! + W2mec
		 W1ineft=W1ineft*emc_func_slac(x,A)
		 W2ineft=W2ineft*emc_func_slac(x,A)
	      endif
	   endif
	   if(A.gt.1.0) then
	      if(wsq.lt.5.0) then
		 innt=30
		 innp=30
		 pmax=1.0
		 if(A.eq.2.0) then
		    eps=0.0022
		    f0=0.009042
		    bigB=0.0008522
		    a1=0.007727
		    b1=0.009394
		    alpha1= 45.384
		 elseif(A.gt.2.0 .and. A.le. 20.0) then !use carbon
		    eps=0.01727
		    f0=0.003182
		    bigB=0.001359
		    a1=0.003027
		    b1=0.007050
		    alpha1=137.2
		 elseif(A.gt.20.0 .and. A.le.80.0) then ! use copper
		    eps=0.00855
		    f0=0.00287403
		    bigB= 0.0008866
		    a1=0.0030959
		    b1=0.0070945
		    alpha1=132.458
		    elseif(A.gt.80.0) then ! use gold
		    eps=0.00693
		    f0=0.0026424
		    bigB= 0.0007632
		    a1=0.0030654
		    b1=0.0067678
		    alpha1=132.452
		 endif
		 call bdisnew4he3(e1,e2,th,A, Z,A-Z, 
     >  	   eps, pmax, innp, innt,  f0, bigB, a1, 
     >  	   b1,alpha1,sig_smear,w1smear,w2smear)
		 w1smear=w1smear*smear_cor(x,A)
		 w2smear=w2smear*smear_cor(x,A)
		 if(wsq.le.3.8) then
		    W1=w1smear
		    W2=w2smear
		 else
		    frac=(wsq-3.8)/(5.0-3.8)
		    W1=(1-frac)*w1smear + W1ineft*frac
		    W2=(1-frac)*w2smear + W2ineft*frac
		 endif
	      else
		 W1=W1ineft
		 W2=W2ineft
	      endif
	   else
	      W1=W1p
	      W2=W2p
	   endif
C       Mott cross section
	   sigmott=(19732.0/(2.0*137.0388*e1*sn**2))**2*cs**2/1.d6
	   sig_dis = 1d3*sigmott*(W2+2.0*W1*tn**2) ! this is nb I think
	endif


	if((xflag.eq.1).or.(xflag.eq.2)) then
	   call F1F2QE09(Z, A, Q2, WSQ, F1, F2)
C       Convert F1,F2 to W1,W2
	   W1 = F1/m_p
	   W2 = F2/nu
C       Mott cross section
	   sigmott=(19732.0/(2.0*137.0388*e1*sn**2))**2*cs**2/1.d6
	   sig_qe = 1d3*sigmott*(W2+2.0*W1*tn**2)
	endif


	sig = sig_qe + sig_dis !sig is already real*4


	sig_qe_pass = sig_qe ! pass back as real*4
	sig_dis_pass = sig_dis

	return
	end


c-------------------------------------------------------------------------------------------
	real*8 function smear_cor(x,A)
C DJG: Correction to smearing calculation
	implicit none
	real*8 A,x,xpass
	real*8 p1,p2,p3,p4
	real*8 x1,x2,xit


	smear_cor=1.0
	x1=0.60
	x2=0.95
	if(A.gt.1.5.and.A.lt.2.5) x2=1.2
	if(x.lt.x1) then
	   xit=x1
	elseif(x.gt.x2) then
	   xit=x2
	else
	   xit=x
	endif

	if(A.gt.1.5.and.A.lt.2.5) then
	   
	   smear_cor=-34.57898661+225.89600206*xit-563.07790319*xit**2+
     >      689.12634955*xit**3-414.41696003*xit**4+98.18408451*xit**5
	   smear_cor=smear_cor/1.025

	elseif(A.gt.2.5.and.A.lt.3.5) then !3He
c	   smear_cor=-58.1433332 + 327.70804494*xit - 672.44236276*xit**2  
c     >      + 608.1593233*xit**3 - 204.73844846*xit**4
	   smear_cor=1.0
	elseif(A.gt.3.5.and.A.lt.4.5) then ! 4He
c	   smear_cor=-33.37855579 + 197.50227515*xit - 417.40813687*xit**2  
c     >      + 387.18924082*xit**3 - 133.36528243*xit**4
	   smear_cor=1.0
	elseif(A.gt.5.5.and.A.lt.6.5) then ! 6Li
c	   smear_cor=65.23349909 - 486.84210205*xit + 1466.93311236*xit**2 
c     >      - 2186.29697553*xit**3 + 1609.31813952*xit**4 
c     >      - 467.95090794*xit**5
	   smear_cor=1.0
	elseif(A.gt.6.5.and.A.lt.7.5) then ! 7Li
c	   smear_cor=63.10711668 - 471.75998378*xit + 1424.51649584*xit**2 
c     >      - 2127.14621541*xit**3 + 1568.13221276*xit**4 
c     >      - 456.46738585*xit**5
	   smear_cor=1.0
	elseif(A.gt.8.5.and.A.lt.9.5) then
c	   smear_cor=73.3171768 - 545.3096107*xit + 1634.67582811*xit**2
c     >      - 2424.82036345*xit**3 + 1776.73735915*xit**4 
c     >      - 514.20267654*xit**5
	   smear_cor=1.0
	else
	   smear_cor = 1.0	! set to 1 by default
	endif

	return
	end

c-------------------------------------------------------------------------------------------
	real*8 function emc_func_slac(xdum,A)
	real*8 x,A,atemp,xdum
	real*8 alpha,C

	atemp = A
!	if(A.eq.4.or.A.eq.3) then  ! emc effect is more like C for these 2...
!	   atemp = 12
!	endif

C use fixed EMC correction starting at x=0.9
	if(xdum.gt.0.9) then
	   x=0.9
	else
	   x=xdum
	endif
	
	if(A.lt.3) then
	   emc_func_slac=1.0
	else
	   
	   alpha = -0.070+2.189*x - 24.667*x**2 + 145.291*x**3
	1	-497.237*x**4 + 1013.129*x**5 - 1208.393*x**6
	2	+775.767*x**7 - 205.872*x**8

	   C = exp( 0.017 + 0.018*log(x) + 0.005*log(x)**2)
	
	   emc_func_slac = C*atemp**alpha
	endif
	return 
	end      
	   

	
      SUBROUTINE MEC2020(z,a,w2,q2,f1mec)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Subroutine for Transverse Enhancement new the QE and Delta due to meson   CCC
CCC   exchange currents and isobar excitations in the medium.  This is assumed  CCC
CCC   to be due to quasi-deuteron 2-body currents.  Shape is a distorted        CCC
CCC   Gaussian in W^2 with a cut-off near the 2-body threshold near x=2.        CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      
! fit to low q2 dip region: purefly empirical
! assume contribution is pure transverse
      implicit none
      real*8 z,a,q2,w2,mp/0.938272/,mp2,w,nu
      real*8 a1,b1,c1,t1,dw2
      real*8 x, f1mec
C use A=12 values
      real*8 xvalm(45) / 
     & 0.40388E+00,0.28053E+01,0.29253E+00,0.60248E+01,0.85436E+00,
     & 0.94254E-01,0.21100E+01,0.12762E+01,0.40249E+00,-.22073E+01,
     & 0.98110E+00,0.95991E+00,0.10305E+01,0.99548E+00,0.97845E+00,
     & 0.10000E+01,0.97772E+00,0.10062E+01,0.99096E+00,0.99684E+00,
     & 0.10051E+01,0.99674E+00,0.10032E+01,0.10070E+01,0.10054E+01,
     & 0.16075E+01,0.23972E+00,0.26369E+01,0.45830E+00,0.91716E+00,
     & 0.37683E-01,0.00000E+00,0.13844E+00,0.12306E+00,0.15551E+01,
     & 0.99849E+02,0.44784E-02,0.10040E-04,0.18568E-01,0.33686E-01,
     & 0.19603E+00,0.24499E+00,0.32459E+00,0.00000E+00,0.10000E+01 /

      mp2 = mp*mp
      f1mec = 0.0
      if(w2.le.0.0) return
      w  = sqrt(w2)
      nu = (w2 - mp2 + q2)/2./mp
      x  = q2/(2.0*mp*nu )

      if(A.lt.2.5) return

      a1 = A*q2**2*xvalm(1)*exp(-1.0*q2/xvalm(2))/
     &                                (xvalm(3)+q2)**xvalm(4)

      b1 = xvalm(5)+xvalm(6)*q2

      c1 = xvalm(33)+xvalm(34)*q2

      t1 = (w2-b1)**2/c1**2/2.

      dw2 = w2+q2*(1.-1./2.2)-1.0*mp2

      if(dw2.LT.0.0) dw2 = 0.0
      f1mec = a1*(exp(-1.*t1)*sqrt(dw2))

      
       if(f1mec.LE.1.0E-9 ) f1mec=0.0


      return
      end


C DJG NOTE 7/11/2007
c Note that although this is called "emc_func", it is really a fit to the ratio of the emc effect to 
c a pure smearing claculation. So - this function is applied to the smeared n+p cross section
c to reproduce the correct emc effect. If you plot these functions, it will not look like
c the emc effect, so don't panic.
  
	real*8 function emc_func_xem(x,A) ! now compute the emc effect from our own fits.
	implicit none
        real*8 x,A,xtmp, Atmp
	real*8 emc
	
	if(A.gt.9 .and. A.le.20) then
	   Atmp=12.0
	elseif(A.gt.20. and. A.le.80.) then
	   Atmp=64.0
	elseif(A.gt.80.0) then
	   Atmp=197
	else
	   Atmp=A
	endif


	
c	if (x.le.1.0) then
	if(x.le.0.9) then
	   xtmp = x
	else
	   xtmp = 0.9
	endif

	emc =1.0
CDeuterium******************************************************************************
	if(Atmp.eq.2) then
C it 2
c	   if(xtmp.lt.0.2) then
c	      emc=1.06
c	   else
c	      emc = 0.79515 +1.9777*xtmp - 3.9724*xtmp**2 -0.66967*xtmp**3
c	1	   +8.3082*xtmp**4 - 5.5263*xtmp**5
c	   endif
c	   emc = emc*0.96689
c it 3
c	   emc = 0.70443 +2.3742*xtmp - 4.6566*xtmp**2 -0.78540*xtmp**3
c     >   	+9.3838*xtmp**4 - 6.1256*xtmp**5
c it4
	   emc = 0.83818 +1.4201*xtmp - 2.8189*xtmp**2 -0.67289*xtmp**3
     >      +6.0543*xtmp**4 - 3.8575*xtmp**5
CHe3***********************************************************************************
	else if(Atmp.eq.3) then
C it 2
c	      emc = 1.0118 +1.1029*xtmp -2.5081*xtmp**2 - 0.22033*xtmp**3
c	1	   + 4.8120*xtmp**4 - 3.2865*xtmp**5
C it 3
c	   emc = 0.92170 +1.7544*xtmp -3.7324*xtmp**2 - 0.24293*xtmp**3
c     >  	+ 6.7613*xtmp**4 - 4.6089*xtmp**5
c it4
	   emc = 0.95649 +1.2695*xtmp -2.3735*xtmp**2 - 0.62320*xtmp**3
     >        + 4.4999*xtmp**4 - 2.7871*xtmp**5
CHe4***********************************************************************************	      
	else if(Atmp.eq.4) then
C it2
c            emc = 0.84622 + 2.2462*xtmp - 4.7909*xtmp**2
c	1	   + 0.065713*xtmp**3 + 7.6154*xtmp**4 - 5.2029*xtmp**5
C it3
c	   emc = 0.70050 + 3.1241*xtmp - 6.1738*xtmp**2
c     >  	- 0.049988*xtmp**3 + 9.3053*xtmp**4 - 6.1348*xtmp**5
C it4
	   emc = 0.85453 + 2.0067*xtmp - 4.0285*xtmp**2
     >       + 0.20116*xtmp**3 + 5.3353*xtmp**4 - 3.5246*xtmp**5
C Be**********************************************************************************
	else if(Atmp.eq.9) then
C it 2
c	      emc = 0.80887 + 3.9354*xtmp - 8.6056*xtmp**2 -0.16342*xtmp**3
c	1	   + 14.074*xtmp**4 -9.3065*xtmp**5
C it 3
c	   emc = 0.46324 + 6.1220*xtmp - 12.184*xtmp**2 -1.0956*xtmp**3
c     >  	+ 20.316*xtmp**4 -12.899*xtmp**5
C it 4
	emc = 0.67858 + 4.5917*xtmp - 9.0413*xtmp**2 -1.4294*xtmp**3
     > 	   + 15.497*xtmp**4 -9.4333*xtmp**5
C Carbon**********************************************************************************
	else if(Atmp.eq.12) then
C it 2
c         emc = 0.8279 + 3.5070*xtmp -7.5807*xtmp**2 
c	1	   -0.60935*xtmp**3 +13.081*xtmp**4 -8.5083*xtmp**5
C it 3
c	   emc = 0.63653 + 4.6458*xtmp -9.2994*xtmp**2 
c     >  	-1.2226*xtmp**3 +16.157*xtmp**4 -10.236*xtmp**5
C it 4
	   emc = 0.72341 + 3.9547*xtmp -7.7589*xtmp**2 
     >      -1.5734*xtmp**3 +14.090*xtmp**4 -8.6972*xtmp**5
C Al**********************************************************************************
	else if(atmp.eq.27) then
	   emc = 0.98645 + 3.0385*xtmp - 22.072*xtmp**2 + 74.981*xtmp**3
     >  	- 132.97*xtmp**4 + 113.06*xtmp**5 -35.612*xtmp**6
C Copper**********************************************************************************
	else if(Atmp.eq.64) then 
C it 2
c	      emc = 1.1075 + 2.7709*xtmp - 6.5395*xtmp**2 -0.46848 *xtmp**3
c	1	   +10.534*xtmp**4 - 6.6257*xtmp**5
c it 3
c	   emc = 0.58372 + 6.0358*xtmp - 11.988*xtmp**2 -1.0211*xtmp**3
c     >  	+18.567*xtmp**4 - 11.482*xtmp**5
c it 4
	   emc = 0.70732 + 5.0635*xtmp - 9.5274*xtmp**2 -1.7399*xtmp**3
     >     +14.835*xtmp**4 - 8.4185*xtmp**5
C Gold**********************************************************************************	      
	else if(Atmp.eq.197) then
C it 2
c	      emc = 1.1404 + 4.0660*xtmp -10.318*xtmp**2 -1.9036*xtmp**3
c	1	   + 21.969*xtmp**4 - 14.461*xtmp**5	   
C it 3
c	   emc = 0.44132 + 8.1232*xtmp -16.141*xtmp**2 -5.6562*xtmp**3
c     >  	+ 35.606*xtmp**4 - 22.008*xtmp**5
C it4
	      emc = 0.79377 + 5.8630*xtmp -11.149*xtmp**2 -5.8010*xtmp**3
     >	   + 26.139*xtmp**4 - 14.992*xtmp**5	
	      
	else  
	   write(*,*) '** in emc_func_xem, unknown target'
	   stop		
	endif
	
	emc_func_xem= emc
	return
	end


	real*8 function inelastic_it(x,A)
C DJG: Correction to inelastic cross section to F1F209 from XEM
C DJG: 40 degree data. Just a simple one-pass iteration for use
C DJG: to check our model dependence.
	implicit none
	real*8 A,x,xpass
	real*8 p1,p2,p3,p4
	real*8 x1,x2,xit


	inelastic_it = 1.0 ! set to 1 by default

	if(A.lt.2.0) then ! hydrogen
c	   if(x.le.1.0 .and. x.ge.0.22) then 
c	      xit=x
c	   elseif(x.gt.1.0) then
c	      xit=1.0
c	   elseif(x.lt.0.22) then
c	      xit=0.22
c	   endif
c	   inelastic_it = 1.1409-1.6098*xit+5.7926*xit**2-8.2497*xit**3
c     >       +4.1543*xit**4
	   inelastic_it=1.0
	endif

	if(A.gt.1.5.and.A.lt.2.5) then !deuterium
	   x1=0.2
	   x2=1.45
	   if(x.lt.x1) then
	      xit=0.2
	   elseif(x.gt.x2) then
	      xit=1.45
	   else
	      xit=x
	   endif

	   if (xit.le.0.99) then
	      inelastic_it=0.7451958+2.61132496*xit-7.45341781*xit**2
     >           +8.96206613*xit**3-4.02966634*xit**4
	   elseif(xit.gt.0.99) then
	      inelastic_it=201.77189575-729.71020935*xit+995.01149727*xit**2
     >          -604.65771117*xit**3+138.44769913*xit**4
	   endif
	endif

	if(A.gt.2.5 .and. A.lt.3.5) then !3He
	   x1=0.6
	   x2=0.9
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   else
	      xit=x
	   endif
	   
	   inelastic_it=-52.67589987+298.51612558*xit-614.09651574*xit**2+
     >         +556.88720307*xit**3-188.16611293*xit**4

	endif

	if(A.gt.3.5 .and. A.lt.4.5) then !4He
c	   x1=0.3172
c	   x2=1.0927
c	   if(x.lt.x1) then
c	      xit=x1
c	   elseif(x.gt.x2) then
c	      xit=x2
c	   elseif(x.ge.x1 .and. x.le.x2) then
c	      xit=x
c	   endif
c	   inelastic_it=1.505d0 - 4.8103d0*xit + 15.221d0*xit**2
c     >     -19.713d0*xit**3 + 8.9693d0*xit**4
	   inelastic_it=1.0
	endif

	if(A.gt.8.5 .and. A.lt.9.5) then !beryllium
	   x1=0.22
	   x2=1.0
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.85082 + 1.4186*xit -4.9446*xit**2 + 7.4273*xit**3
     >             -3.9151*xit**4
	endif

	if(A.gt.9.5 .and. A.lt.10.5) then !boron-10
	   x1=0.22
	   x2=1.0166
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.98334 + 0.51046*xit - 2.8438*xit**2
     >        +5.2527*xit**3 -3.0320*xit**4
	endif

	if(A.gt.10.5 .and. A.lt.11.5) then !boron-11
	   x1=0.22
	   x2=1.0166
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.93151 + 0.95067*xit - 3.9928*xit**2
     >        +6.4022*xit**3 -3.3918*xit**4
	endif

	if(A.gt.11.5 .and. A.lt.13.5) then ! carbon
	   x1=0.22
	   x2=0.966
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.97414 + 0.303133*xit - 1.8188*xit**2
     >        +3.6032*xit**3 -2.1878*xit**4
	endif

	if(A.gt.61.0 .and. A.lt.66.0) then !copper
	   x1=0.3172
	   x2=1.005
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.77646d0 + 0.90481d0*xit - 0.83815d0*xit**2
	endif


	if(A.gt.196.0 .and. A.lt.198.0) then !gold
	   x1=0.3172
	   x2=1.005
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.70439d0 + 1.0510d0*xit - 0.91679d0*xit**2
	endif

	return
	end
