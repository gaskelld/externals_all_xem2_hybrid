      subroutine get_cc_info(iA, iZ, deltae_cc)
C     DG note
c     subroutine to get deltae for Coulomb Corrections
c     just making a subroutine so we only need to update
c     stuff in one place

c      Ai Daniel note:
c     boost to beam energy and scattered electron energy
c     and do Coulomb correction.  recipe is given in:
c     Aste et al., Eur. Phys. J. A 26, 167 (2005)
c     note that these deltae_cc's are calculated as follows:
c     deltae_cc = (1.5)*(Z/R)*(hbar*c)*alpha*0.775
c     the 0.775 comes from Aste EPJ A 26 167 (2005)
c     also, all deltae_cc's are computed for Z-1, not Z!/-* 

      integer iA, iZ
      real deltae_cc,Re
      real alpha,hbarc

      parameter(alpha=1.0/137.036)
      parameter(hbarc=0.197327) ! GeV fm

      if(iA.eq.1) then          !Hydrogen
         deltae_cc = 0.0
      elseif(iA.eq.2) then      !Deuterium
         deltae_cc = 0.0
      elseif(iA.eq.3) then      !Helium-3
         deltae_cc = 0.00085
      elseif(iA.eq.4) then      !Helium-4
         deltae_cc = 0.0010
      elseif(iA.eq.9) then      !Beryllium	  
         deltae_cc = 0.001875
      elseif(iA.eq.12) then     !Carbon
         deltae_cc = 0.00292
      elseif(iA.eq.27) then     !Aluminum
         deltae_cc = 0.0061
      elseif(iA.eq.64) then     !Copper
         deltae_cc = 0.0102
      elseif(iA.eq.197) then    !Gold	  
         deltae_cc = 0.0199
      else                      ! use empirical fit
         Re=1.1*iA**(1./3.)+0.86*iA**(-1./3.)
         deltae_cc = 1.5*(iZ-1.0)/Re*(hbarc)*alpha*0.775
      endif

        return
        end
