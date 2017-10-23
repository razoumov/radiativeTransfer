c=======================================================================
c//////////////////////  SUBROUTINE CALC_RATES  \\\\\\\\\\\\\\\\\\\\\\\\c
      subroutine calc_rates(nratec, temstart, temend,
     &             ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &             ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &             reHeII2a, reHeIIIa, brema, lineHIa, compa,
     &             hyd01ka, h2k01a, vibha, rotha, rotla, 
     &             gpldl, gphdl, hdltea, hdlowa,
     &             k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &             k11a, k12a, k13a, k13dda, k14a, k15a, k16a, k17a,
     &             k18a, k19a, k20a, k21a, k22a,
     &             k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &             recombinationType
     &                     )
c
c  COMPUTE MULTISPECIES RATE LOOKUP TABLE
c
c  written by: Yu Zhang, Peter Anninos and Tom Abel
c  date:       
c  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
c  modified2: October, 1996 by GB; moved to AMR
c  modified3: February, 2000 by GB; added Toms new collisional rates
c  modified4: 2005 by Alexei Razoumov; adapted to FTTE
c
c  PURPOSE:
c    Construct tables for rate coefficients and some cooling functions 
c    versus the logarithm of temperature.
c
c  UNITS:
c    Most rate coefficients are computed in cgs units and then converted
c    into more convenient units.  This would usually be code units
c    (see units.src) but in many cases the natural code units are not
c    near unity (e.g. particles/box volume), so we include units and
c    constants from the equations themselves.  This is explained for
c    each of two types used: rate coefficients, cooling coefficients.
c    Note also that the densities in the rate equations are true, not
c    comoving densities and so we must include a factor of a^3.
c    NOTE: if the spectrum changes, this routine must be called again.
c
c  PARAMETERS:
c    nnu      = number of frequency bins
c    everg    = ergs per eV (Rydberg constant) in cgs
c    evhz     = everg / h
c    mh       = mass of hydrogen atom
c    temstart = table start temperature (K)
c    temend   = table end temperature (K)
c    tevk     = conversion constant [eV] = [K] / tevk
c    ireco    = Recombination cooling flag (1 = on)
c    ibrco    = Brehmmstrahlung cooling flag (1 = on)
c    icico    = Collisional ionization flag (1 = on)
c    iceco    = Collisional excitation flag (1 = on)
c    iphoto_mol = Molecular photo-dissociation flag (1 = on)
c
c  INPUTS
c    nratec   = number of temperature bins in rate equations
c
c    tbase1   = code time units
c    xbase1   = code length units
c    kunit    = conversion factor xbase1**3/tbase1
c
c  RATES:
c
c  (New rates numbering as of Feb/2000, agrees with original
c   Abel etal 1997)
c    old   -> new              old   -> new
c    k1-10 -> k1-10            k16   -> k14
c    k11   -> k13              k17   -> k15
c    k12   -> k11              k18   -> k16
c    k13   -> removed          k19   -> k17
c    k14   -> k12              k20   -> k18
c    k15   -> removed          k21   -> k19
c
C ------- 1:    HI    + e   -> HII   + 2e
C ------- 2:    HII   + e   -> H     + p
C ------- 3:    HeI   + e   -> HeII  + 2e
C ------- 4:    HeII  + e   -> HeI   + p
C ------- 5:    HeII  + e   -> HeIII + 2e
C ------- 6:    HeIII + e   -> HeII  + p
C ------- 7:    HI    + e   -> HM    + p
C ------- 8:    HM    + HI  -> H2I*  + e
C ------- 9:    HI    + HII -> H2II  + p
C ------- 10:   H2II  + HI  -> H2I*  + HII
c
c
C --------------  old  ---------------------
C ------- 11:   H2I   + H   -> 3H
C ------- 12:   H2I   + HII -> H2II  + H
C ------- 13:   H2I   + e   -> HI    + HM
C ------- 14:   H2I   + e   -> 2HI   + e
C ------- 15:   H2I   + H2I -> H2I   + 2HI
C ------- 16:   HM    + e   -> HI    + 2e
C ------- 17:   HM    + HI  -> 2H    + e
C ------- 18:   HM    + HII -> 2HI
C ------- 19:   HM    + HII -> H2II  + e
C ------- 20:   H2II  + e   -> 2HI
C ------- 21:   H2II  + HM  -> HI    + H2I
C --------------  old  ---------------------
c
c --------------  new  ---------------------
C ---11--       H2I   + HII -> H2II  + H
C ---12--       H2I   + e   -> 2HI   + e
C ---13--       H2I   + H   -> 3H
C ---14--       HM    + e   -> HI    + 2e
C ---15--       HM    + HI  -> 2H    + e
C ---16--       HM    + HII -> 2HI
C ---17--       HM    + HII -> H2II  + e
C ---18--       H2II  + e   -> 2HI
C ---19--       H2II  + HM  -> HI    + H2I
c (20-21) - unused
c --------------  new  ---------------------
c
c
C ------- 22:   2H    + H   -> H2I   + H
C ------- 24:   HI    + p   -> HII   + e
C ------- 25:   HeII  + p   -> HeIII + e
C ------- 26:   HeI   + p   -> HeII  + e
C ------- 27:   HM    + p   -> HI    + e
C ------- 28:   H2II  + p   -> HI    + HII
C ------- 29:   H2I   + p   -> H2II  + e
C ------- 30:   H2II  + p   -> 2HII  + e
C ------- 31:   H2I   + p   -> 2HI
C
c ------- 50-56: deuterium rates (given below)
c
c-----------------------------------------------------------------------
c
c  This define indicates if the Tom Abels new (as of Feb/2000) collisional
c    rates (k1-k19) should be used.
c
      implicit none
      integer nratec, recombinationType, nnu
      double precision temstart, temend, utim
      double precision ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &     ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &     ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &     reHeII2a(nratec), reHeIIIa(nratec), brema(nratec),
     &     lineHIa(nratec)
      double precision hyd01ka(nratec), h2k01a(nratec),
     &     vibha(nratec),
     &     rotha(nratec), rotla(nratec), gpldl(nratec),
     &     gphdl(nratec), hdltea(nratec), hdlowa(nratec)
      double precision compa
      double precision k1a (nratec), k2a(nratec), k3a(nratec),
     &     k4a(nratec),
     &        k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &        k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &        k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &        k17a(nratec), k18a(nratec), k19a(nratec), k20a(nratec), 
     &        k21a(nratec),
     &        k22a(nratec), k50a(nratec), k51a(nratec), k52a(nratec), 
     &        k53a(nratec), k54a(nratec), k55a(nratec), k56a(nratec)
      double precision k13dda(nratec, 7)
      parameter (nnu = 400)
      integer ireco, ibrco, icico, iceco, iphoto_mol
      parameter (ireco = 1, ibrco = 1, icico = 1, iceco = 1, 
     &           iphoto_mol = 0)
      double precision everg, evhz, tevk, mh, pi, tiny_cgs, dhuge, kb
      parameter (everg = 1.60184d-12, evhz  = 2.41838d14, 
     &           tevk = 1.16d+4, mh = 1.67d-24, pi=3.14159,
     &           tiny_cgs = 1.0d-37, dhuge = 1.0d+30,
     &           kb = 1.380658d-16)
c
c         Set various cutoff values in eV.
c
      double precision e24,e25,e26,e27,e28a,e28b,e28c,e29a,e29b,e29c,
     &                 e30a,e30b
        PARAMETER(
     &    e24  = 13.6
     &   ,e25  = 54.4
     &   ,e26  = 24.6
     &   ,e27  = 0.755
     &   ,e28a = 2.65
     &   ,e28b = 11.27
     &   ,e28c = 21.0
     &   ,e29a = 15.42
     &   ,e29b = 16.5
     &   ,e29c = 17.7
     &   ,e30a = 30.0
     &   ,e30b = 70.0)

      integer caseA, caseB
      parameter (caseA = 1, caseB = 2)

c  Locals
c
      integer i,j,n
      double precision nu, delnu, hnu, logttt, ttt, log10ttt,
     &     tev, logtev,
     &     xx, dum, kunit, coolunit, dlogtem, kunit_3bdy
      double precision sigma24(nnu), sigma25(nnu), sigma26(nnu), 
     &     sigma27(nnu), sigma28(nnu), sigma29(nnu), 
     &     sigma30(nnu), sigma31(nnu)
      double precision tm
      double precision HDLR, HDLV, lt, t3
      double precision tempMellema(81), tmp, HIIrecCoolMellema(81)
      double precision tempGnedin(201), tmp1, HeIIrecCoolGnedin(201),
     &     tmp2, HeIIIrecCoolGnedin(201)
c
c  Get conversion units
c
c     t/x/dbase1 is the number (z dependant) that converts from the
c       dimensionless code units to physical units.  Also, in the
c       code aye = 1 at z=zinit, so to convert the usual a (=1 at z=0)
c       to a~ (written in the code as aye), we use a = a~*[a] 
c
c  1) Set the dimensions of the (non-radiative) rate coefficients.  
c    Note that we have included the units that convert density to 
c    number density, so the rate equations should look like 
c    (in dimensionless units, hence the primes):
c
c       d(d0~)/dt~ = k~ * d1~ * d2~ / a~^3
c
c    where k~ is the dimenionless rate coefficients and d0-2~ are three
c     dimensionless densities (i.e. d = [dens]*d~) and a~ is the 
c     dimensionless expansion coefficient (see above).
c
c    rate eqn        : delta(n0)  = k  * n1        * n2        * dt     / a^3
c    rate eqn units  : [dens]/mh  = k  * [dens]/mh * [dens]/mh * [time] / [a]^3
c    rate eqn dimless: delta(n0~) = k~ * n1~       * n2~       * dt~    / a~^3
c    so: k = [k] * k~  where [k] = ( [a]^3 * mh ) / ( [dens] * [time] )  (~)
c    reminder: the number densities here are normalized with [dens] which
c              is not a constant (it has a factor a^3), so the number
c              densities must be converted from comoving to proper.
c
c
c  2) Set the dimension of the cooling coefficients (including constants)
c     (this equation has a rho because e is the specific energy, not
c      energy/unit volume).
c       delta(e)  = L     * n1        * n2        * dt     / dens   / a^3
c       [e]       = L     * [dens]/mh * [dens]/mh * [time] / [dens] / [a]^3
c       delta(e~) = L~    * n1~       * n2~       * dt~    / dens~  / a~^3 [~]
c     so L = [L] * L~ where [L] = [e] * mh**2 * [a]^3 / ([dens] * [time]) [~]
c       but [e] = ([a]*[x])**2 / [time]**2 and ([a] = 1 / (1 + zri) )
c      [L] = ([a]**5 * [x]**2 * mh**2) / ([dens] * [time]**3)
c
c
c    Note: some of the coffiecients have only one power of n.  These
c          do not have the /a^3 factor, also they have units
c          [L1] = ([a]**2 * [x]**2 * mh) / [time]**3
c               = [L] * [dens] * [a]**3 / mh
c          This is done through the dom variable in cool.src
c         (some have three powers of n and they are different by the
c          reciprocal of the above factor multiplying [L]).
c
c  3) the units for the radiative rate coefficients is just 1/[time]
c
c  Compute log spacing in temperature
c
      ttt    = temstart
      logttt = dlog(ttt)
      dlogtem= (log(temend) - log(temstart))/real(nratec-1)
c
c  Initialize constants to 1.0e-20
c
      do i = 1, nratec
c
        k1a (i) = 1.0e-20
        k2a (i) = 1.0e-20
        k3a (i) = 1.0e-20
        k4a (i) = 1.0e-20
        k5a (i) = 1.0e-20
        k6a (i) = 1.0e-20
        k7a (i) = 1.0e-20
        k8a (i) = 1.0e-20
        k9a (i) = 1.0e-20
        k10a(i) = 1.0e-20
        k11a(i) = 1.0e-20
        k12a(i) = 1.0e-20
        k13a(i) = 1.0e-20
        k14a(i) = 1.0e-20
        k15a(i) = 1.0e-20
        k16a(i) = 1.0e-20
        k17a(i) = 1.0e-20
        k18a(i) = 1.0e-20
        k19a(i) = 1.0e-20
        k20a(i) = 1.0e-20
        k21a(i) = 1.0e-20
        k22a(i) = 1.0e-20
        k50a(i) = 1.0e-20
        k51a(i) = 1.0e-20
        k52a(i) = 1.0e-20
        k53a(i) = 1.0e-20
        k54a(i) = 1.0e-20
        k55a(i) = 1.0e-20
        k56a(i) = 1.0e-20
c
        do j = 1, 7
           k13dda(i,j) = 1.0e-20
        enddo
c
        ceHIa(i)    = 1.0e-20
        ceHeIa(i)   = 1.0e-20
        ceHeIIa(i)  = 1.0e-20
        ciHIa(i)    = 1.0e-20
        ciHeIa(i)   = 1.0e-20
        ciHeISa(i)  = 1.0e-20
        ciHeIIa(i)  = 1.0e-20
        reHIIa(i)   = 1.0e-20
        reHeII1a(i) = 1.0e-20
        reHeII2a(i) = 1.0e-20
        reHeIIIa(i) = 1.0e-20
        brema(i)    = 1.0e-20
        lineHIa(i)   = 1.0e-20
        compa       = 1.0e-20
c
        hyd01ka(i)  = 1.0e-20
        h2k01a(i)   = 1.0e-20
        vibha(i)    = 1.0e-20
        rotha(i)    = 1.0e-20
        rotla(i)    = 1.0e-20
        hdltea(i)   = 1.0e-20
        hdlowa(i)   = 1.0e-20
c
      enddo
c
c  Fill in tables over the range temstart to temend
c
c -------------------------------------------------
c  1) rate coefficients (excluding external radiation field)
c
      do i = 1, nratec
c
c       Compute temperature of this bin (in eV)
c
        logttt = log(temstart) + real(i-1)*dlogtem
        ttt = exp(logttt)
        tev = ttt/tevk
        logtev = dlog(tev)
c

c
c       Call Tom Abels routine (from his web page, Feb/2000)
c
        call coll_rates(ttt, k1a(i), k2a(i), k3a(i), k4a(i), k5a(i),
     &                  k6a(i), k7a(i), k8a(i), k9a(i), k10a(i),
     &                  k11a(i), k12a(i), k13a(i), k14a(i), k15a(i), 
     &                  k16a(i), k17a(i), k18a(i), k19a(i),
     &                  recombinationType)
c

c
c       Compute the density-dependant collision H2 dissociation rates.
c          (this givens a 7-variable function that is used to compute
c           the log10 of the rate -- normalize by dividing by kunit).
c
        call colh2diss(ttt, k13dda(i,1), k13dda(i,2), k13dda(i,3), 
     &                 k13dda(i,4), k13dda(i,5), k13dda(i,6), 
     &                 k13dda(i,7))
        k13dda(i,1) = k13dda(i,1)
c
c       ------ 3-body H2 rate ----
c       The first bit is my fit to A.E. Orel 1987, J.Chem.Phys., 87, 
c       314. I then match it to the T^-1 of Palla etal (1983) 
c       Which is then 4 times smaller than the Palla rate ! Thats
c       molecule uncertainties ! :-)
c
        if (ttt .le. 300.0) then
           k22a(i) = 1.3d-32 * (ttt/300.0)**(-0.38)
        else
           k22a(i) = 1.3d-32 * (ttt/300.0)**(-1.0)
        endif
c
c       ------ Deuterium rates -----
c
c       50) H+  +  D  -> H  +  D+
c       51) H   +  D+ -> H+ +  D
c       52) H2  +  D+ -> HD +  H+
c       53) HD  +  H+ -> H2 +  D+
c       54) H2  +  D  -> HD +  H
c       55) HD  +  H  -> H2 +  D
c       56) D   +  H- -> HD +  e-
c      [57) D-  +  H  -> HD +  e-]  included by multiply 56 by 2
c
        k50a(i) = 1.0d-9 *exp(-4.1d1/ttt)
        k51a(i) = 1.0d-9
        k52a(i) = 2.1d-9
        k53a(i) = 1.0d-9 *exp(-4.57d2/ttt)
        k54a(i) = 7.5d-11*exp(-3.82d3/ttt)
        k55a(i) = 7.5d-11*exp(-4.24d3/ttt)
        k56a(i) = 1.5d-9 *(ttt/300.0)**(-0.1)
c
      enddo
c
c  Write out cgs rate coefficients
c
!      open(10, file='k1-6.out', status='unknown')
!      do i = 1, nratec
!        write(10,1010) exp(log(temstart) + real(i-1)*dlogtem),
!     &       k1a(i),k2a(i),k3a(i),k4a(i),k5a(i),k6a(i)
!      enddo
! 1010 format(7(1pe11.4))
!      close(10)
c
c -------------------------------------------------
c  2) Cooling/heating rates (excluding external radiation field)

      open(unit=1,file='HII-ktbetas.tab',action="read",status="old")
      do i = 1, 81
        read(1,*) tempMellema(i), tmp, HIIrecCoolMellema(i)
      enddo
      close(1)

      open(unit=2,file='cratesHe.res',action="read",status="old")
      do i = 1, 201
        read(2,*) tempGnedin(i), tmp1, HeIIrecCoolGnedin(i),
     &       tmp2, HeIIIrecCoolGnedin(i)
        tempGnedin(i) = log10(tempGnedin(i))
        HeIIrecCoolGnedin(i) = log10(HeIIrecCoolGnedin(i))
        HeIIIrecCoolGnedin(i) = log10(HeIIIrecCoolGnedin(i))
      enddo
      close(2)

c
      do i = 1, nratec
c
        logttt = log(temstart) + real(i-1)*dlogtem
        ttt = exp(logttt)
        log10ttt = log10(ttt)
c
c    a) Collisional excitations (Black 1981; Cen 1992)
c
      if ( iceco .eq. 1 ) then
        ceHIa(i)    = 7.5d-19*exp(-dmin1(dlog(dhuge),118348./ttt))
     &              /(1.+dsqrt(ttt/1.0d5))
        ceHeIa(i)   = 9.1d-27*exp(-dmin1(dlog(dhuge),13179./ttt))
     &              *ttt**(-0.1687)/(1.+dsqrt(ttt/1.0d5))
        ceHeIIa(i)  = 5.54d-17*exp(-dmin1(dlog(dhuge),473638./ttt))
     &              *ttt**(-0.397)/(1.+sqrt(ttt/1.0d5))
      else
        ceHIa(i)    = 1.0e-20
        ceHeIa(i)   = 1.0e-20
        ceHeIIa(i)  = 1.0e-20
      endif
c
c    b) Collisional ionizations (Cen 1992 or Abel 1996)
c
      if ( icico .eq. 1 ) then
c
c        ciHIa(i)    = 1.27d-21*sqrt(ttt)/(1.+sqrt(ttt/1.0d5))
c     &              *exp(-dmin1(dlog(dhuge),157809.1/ttt))
c        ciHeIa(i)   = 9.38d-22*sqrt(ttt)/(1.+sqrt(ttt/1.0d5))
c     &              *exp(-dmin1(dlog(dhuge),285335.4/ttt))
c        ciHeIIa(i)  = 4.95d-22*sqrt(ttt)/(1.+sqrt(ttt/1.0d5))
c     &              *exp(-dmin1(dlog(dhuge),631515.0/ttt))
c
        ciHeISa(i)  = 5.01d-27*(ttt)**(-0.1687)/(1.+sqrt(ttt/1.0d5))
     &              *exp(-dmin1(dlog(dhuge),55338./ttt))
c
c       Collisional ionizations (polynomial fits from Tom Abel)
c
        ciHIa(i)    = 2.18d-11*k1a(i)
        ciHeIa(i)   = 3.94d-11*k3a(i)
        ciHeIIa(i)  = 8.72d-11*k5a(i)
      else
        ciHIa(i)    = 1.0e-20
        ciHeIa(i)   = 1.0e-20
        ciHeIIa(i)  = 1.0e-20
        ciHeISa(i)  = 1.0e-20
      endif
c
c    c) Recombinations (Cen 1992)
c

      if ( ireco .eq. 1 ) then

!       HII recombination cooling rates
        if (recombinationType.eq.caseA) then ! case A - from Cen 1992
          reHIIa(i)   = 8.70d-27*sqrt(ttt)*(ttt/1000.0)**(-0.2)
     &         / (1.0 + (ttt/1.0d6)**(0.7))
        else                    ! case B - from Garrelt Mellema
          if (log10ttt.lt.tempMellema(1).or.log10ttt.gt.tempMellema(81))
     &         then
            reHIIa(i) = 0.
          else
            j = 1
            do while (log10ttt.gt.tempMellema(j))
              j = j + 1
            enddo
            reHIIa(i) = 10.**((log10ttt-tempMellema(j-1))/
     &           (tempMellema(j)-tempMellema(j-1))*
     &           (HIIrecCoolMellema(j)-HIIrecCoolMellema(j-1))+
     &           HIIrecCoolMellema(j-1))
          endif
        endif

!       HeII and HeIII recombination cooling rates [erg cm^3/s]
        if (recombinationType.eq.caseA) then ! case A - from Cen 1992
          reHeII1a(i) = 1.55d-26*ttt**0.3647
          reHeIIIa(i) = 3.48d-26*sqrt(ttt)*(ttt/1000.0)**(-0.2)
     &         / (1.0 + (ttt/1.0d6)**(0.7))
        else ! case B - from Nick Gnedin
          if (log10ttt.lt.tempGnedin(1).or.log10ttt.gt.tempGnedin(201))
     &         then
            reHeII1a(i) = 0.
            reHeIIIa(i) = 0.
          else
            j = 1
            do while (log10ttt.gt.tempGnedin(j))
              j = j + 1
            enddo
            reHeII1a(i) = 10.**((log10ttt-tempGnedin(j-1))/
     &           (tempGnedin(j)-tempGnedin(j-1))*
     &           (HeIIrecCoolGnedin(j)-HeIIrecCoolGnedin(j-1))+
     &           HeIIrecCoolGnedin(j-1))
            reHeIIIa(i) = 10.**((log10ttt-tempGnedin(j-1))/
     &           (tempGnedin(j)-tempGnedin(j-1))*
     &           (HeIIIrecCoolGnedin(j)-HeIIIrecCoolGnedin(j-1))+
     &           HeIIIrecCoolGnedin(j-1))
          endif
        endif

!       HeII dielectronic recombination cooling rate [erg cm^3/s]
        reHeII2a(i) = 1.24d-13*ttt**(-1.5)
     &              *exp(-dmin1(dlog(dhuge),470000./ttt))
     &              *(1.+0.3*exp(-dmin1(dlog(dhuge),94000./ttt)))

      else
        reHIIa(i)   = 1.0e-20
        reHeII1a(i) = 1.0e-20
        reHeII2a(i) = 1.0e-20
        reHeIIIa(i) = 1.0e-20
      endif
c
c    d) Bremsstrahlung (Black 1981)(Spitzer & Hart 1979)
c
      if ( ibrco .eq. 1 ) then
        brema(i)    = 1.43d-27*sqrt(ttt)
     &       *(1.1+0.34*exp(-(5.5-log10(ttt))**2/3.0))
      else
        brema(i)    = 1.0e-20
      endif
c
c         Bremsstrahlung (Shapiro & Kang 1987)(Spitzer 1978)
c
c        if(ttt .lt. 1.0*3.2d5) gaunt = 0.79464 + 0.1243*log10(ttt/1.)
c        if(ttt .ge. 1.0*3.2d5) gaunt = 2.13164 - 0.1240*log10(ttt/1.)
c        brem1a(i) = 1.426d-27*sqrt(ttt)*gaunt
c        if(ttt .lt. 4.0*3.2d5) gaunt = 0.79464 + 0.1243*log10(ttt/4.)
c        if(ttt .ge. 4.0*3.2d5) gaunt = 2.13164 - 0.1240*log10(ttt/4.)
c        brem2a(i) = 1.426d-27*sqrt(ttt)*gaunt

!     HI line excitation cooling from hui.97
      tmp = 2. * 13.598 * 1.60217646d-12 / (1.3806503d-16 * ttt)
      lineHIa(i) = 7.5e-19 * exp(-0.75*tmp/2.) / (1.+sqrt(ttt/1.e5)) ! [erg cm^3 / s]

c    e) Molecular hydrogen cooling
c
c       (Which one is used is controlled by a flag in cool1d_mulit.src)
c
c    e part1) - Lepp & Shull rates
c
        xx = log10(ttt/1.0d4)
        vibha(i) = 1.1d-18*exp(-dmin1(dlog(dhuge),6744.0/ttt))
c
        if( ttt .gt. 1635) then
           dum = 1.0d-12*sqrt(ttt)*exp(-1000.0/ttt)
        else
           dum = 1.4d-13*exp((ttt/125.) - (ttt/577.)**2)
        endif
        hyd01ka(i) = dum*exp(-dmin1(dlog(dhuge),
     &                            8.152d-13/(1.38d-16*ttt)))
c
        dum = 8.152d-13*(4.2/(1.38d-16*(ttt+1190))+1./(1.38d-16*ttt))
        h2k01a(i) = 1.45d-12*sqrt(ttt)*exp(-dmin1(dlog(dhuge),dum))
c
        if(ttt .gt. 4031)then
          rotla(i) = 1.38d-22*exp(-9243.0/ttt)
        else
          rotla(i) = 10.0**(-22.9 - 0.553*xx - 1.148*xx*xx)
        endif
c
        if(ttt .gt. 1087)then
          rotha(i) = 3.9d-19*exp(-6118.0/ttt)
        else
          rotha(i) = 10.0**(-19.24 + 0.474*xx - 1.247*xx*xx)
        endif
c
c     e part2) - Galli and Palla (1999) rates as fit by Tom Abel
c
        tm = max(ttt, 13.d0)       ! no cooling below 13 Kelvin ...
        tm = min(tm, 1.d5)      ! fixes numerics
        lt = log10(tm)
c     low density limit from Galli and Palla
        gpldl(i)  = 10.**(-103.0+97.59*lt-48.05*lt**2+10.80*lt*lt*lt
     &       -0.9032*lt*lt*lt*lt)
c     high density limit from HM79      
        t3 = tm/1000.
        HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*
     &            exp(-(0.13/t3)**3)+3.e-24*exp(-0.51/t3))
        HDLV = (7.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3))
        gphdl(i)  = (HDLR + HDLV)
c
c     f) HD cooling
c
c        HD Cooling Function (ergs cm3 /s)
c
c Fit to Lepp and Shull 1984 LTE (ergs/s) -> hdlte (ergs cm3/s)
c
      hdltea(i) = -35.6998d0 + 15.35716d0*dlog10(ttt) -
     &            5.58513d0   * (dlog10(ttt))**2 +
     &            0.8561149d0 * (dlog10(ttt))**3 - 
     &            1.75538d-2  * (dlog10(ttt))**4
      hdltea(i) = (10.0**hdltea(i))
c      hdlte=10**lhdlte/den
c
c Galli and Palla 1998 low density limit (erg cm3 /s)
c  uses HD-He collisional data, so reduced by 1.27
c
      hdlowa(i) = ((3.0d0 * (4.4d-12 + 3.6d-13*ttt**0.77) *
     &            dexp(-128.d0/ttt) * 128.d0 +
     &            (5.0d0/3.0d0) * (4.1d-12+2.1e-13*ttt**0.92) *
     &            dexp(-255.d0/ttt) * 255.d0) * kb/1.27d0)
c      hdcool=hdlte/(1.0d0+hdlte/hdlow)
c
      enddo
c
c     g) Compton cooling
c
      compa   = 5.65d-36
c
c -------------------------------------------------
c  3) External radiative processes
c
c     Initialize to 1.0e-20
c
c
c     Loop over all frequency bins, compute cross-sections
c        nu is in eV (range is 0.74 eV to 7.2 keV if nnu=400)
c
      do n = 1, nnu
        nu = 10**((n-14)*0.01)                    ! 0.74 ev -- 7.24d9 ev
c
c       24) HI photo-ionization cross-section
c
        if (nu .gt. e24) then
           dum = sqrt(nu/e24-1)
           sigma24(n) = 6.3d-18 * (e24/nu)**4
     &               *  exp(4.0-4.0*atan(dum)/dum)
     &               / (1-exp(-2.0*pi/dum))
        else
           sigma24(n) = 0.0
        endif
c
c       25) HeII photo-ionization cross-section
c
        if (nu .gt. e25) then
           dum = sqrt(nu/e25-1)
           sigma25(n) = 1.58d-18 * (e25/nu)**4
     &               * exp(4.0-4.0*atan(dum)/dum)
     &               / (1-exp(-2.0*pi/dum))
        else
           sigma25(n) = 0.0
        endif
c
c       26) HeI photo-ionization cross-section
c
        if (nu .gt. e26) then
           sigma26(n) = 7.42d-18*(1.66*(nu/e26)**(-2.05)
     &                          - 0.66*(nu/e26)**(-3.05))
        else
           sigma26(n) = 0.0
        endif
c
c       27) HM    + p   -> HI    + e
c
        if (nu .gt. e27) then
          sigma27(n) = 2.11d-16*(nu-e27)**1.5/nu**3
        else
          sigma27(n) = 0.0
        endif
c
c       28) H2II  + p   -> HI    + HII
c
        if (nu .gt. e28a .and. nu .le. e28b) then
          sigma28(n) = 10**(-40.97+6.03*nu-0.504*nu**2+1.387d-2*nu**3)
        elseif (nu .gt. e28b .and. nu .lt. e28c) then
          sigma28(n) = 10**(-30.26+2.79*nu-0.184*nu**2+3.535d-3*nu**3)
        else
          sigma28(n) = 0.0
        endif
c
c       29) H2I   + p   -> H2II  + e
c
        if (nu .gt. e29a .and. nu .le. e29b) then
          sigma29(n) = 6.2d-18*nu - 9.4d-17
        elseif (nu .gt. e29b .and. nu .le. e29c) then
          sigma29(n) = 1.4d-18*nu - 1.48d-17
        elseif (nu .gt. e29c) then
          sigma29(n) = 2.5d-14*nu**(-2.71)
        else
          sigma29(n) = 0.0
        endif
c
c       30) H2II  + p   -> 2HII  + e
c
        if (nu .ge. e30a .and. nu .lt. e30b) then
          sigma30(n) = 10**(-16.926-4.528d-2*nu+2.238d-4*nu**2
     &                      +4.245d-7*nu**3)
        else
          sigma30(n) = 0.0
        endif
c
c       31) H2I   + p   -> 2HI
c
c        sigma31(n) = 0.0
        if (nu .ge. e28b .and. nu .lt. e24) THEN
           sigma31(n) = 3.71e-18
        else
           sigma31(n) = 0.0
        endif

      enddo
c
c     Integrate over the frequency spectrum
c
      do n = 2, nnu
          delnu = (10**((n-14)*0.01) - 10**((n-15)*0.01))*evhz
          hnu   = 0.5*(10**((n-14)*0.01) + 10**((n-15)*0.01))*everg

c            Add to HI,HeI,HeII photoionization heating coefficients

          dum = hnu - e24*everg
c
          dum = hnu - e26*everg
c
          dum = hnu - e25*everg

      enddo
c
c     The secondary electrons produced by steep spectral sources
c       are accounted for here using data from Shull & Steenberg.
c       Note that this is only good for ionization fractions < 1e-3 or so
c          (now down in the solvers to cover all ionization fractions)
c
      open(10, file='cool_rates.out', status='unknown')
      do i=1, nratec
         logttt = log(temstart) + real(i-1)*dlogtem
         ttt = exp(logttt)
         write(10,1000) ttt, ceHIa(i), ceHeIa(i), ceHeIIa(i), ciHIa(i),
     &                  ciHeIa(i), ciHeISa(i), ciHeIIa(i), reHIIa(i),
     &                  reHeII1a(i), reHeII2a(i), reHeIIIa(i),
     &                  brema(i), compa       
 1000    format(1p,30(e10.3,1x))
      enddo
      close(10)
c
      open(10, file='rates.out', status='unknown')
      do i=1, nratec
         logttt = log(temstart) + real(i-1)*dlogtem
         ttt = exp(logttt)
         write(10,1000) ttt, k1a(i), k2a(i), k3a(i), k4a(i), k5a(i),
     &                  k6a(i), k7a(i), k8a(i), k9a(i), k10a(i),
     &                  k11a(i), k12a(i), k13a(i), k14a(i), k15a(i), 
     &                  k16a(i), k17a(i), k18a(i), k19a(i), k22a(i)
      enddo
      close(10)

      return
      end
