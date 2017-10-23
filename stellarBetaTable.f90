subroutine stellarBetaTable(nfreq,freqdel,totalIntegral,iSpectrum,coefSpectrum,iMetal,coefMetal)

  use definitions
  use stellarPopulationModule, only: stellarPopulation
  use dust

  implicit none
  integer, intent(in) :: nfreq, iSpectrum, iMetal
  real(kind=RealKind), intent(in) :: freqdel, coefSpectrum, coefMetal
  real(kind=RealKind), intent(out) :: totalIntegral

  integer :: i, irate, ier, iWavelength, idepth1, idepth2, idepth3, idepthDust, ienergy
  real(kind=RealKind) :: delta_nu, dtmp, tmp1, tmp2, tmp3, &
       thisWavelength, thisSpecificLuminosity, freq, dum, &
       denominatorForBeta1, denominatorForBeta2, denominatorForBeta3, coefWavelength, &
       tau1, tau2, tau3, tauDust, atmp, lambda
  real(kind=RealKind), dimension(nfreq) :: sigma24, sigma25, sigma26, &
       sigma27, sigma28, sigma29, sigma30, sigma31, nu, sigmaDust

  real(kind=RealKind), parameter :: &
       e27  = 0.755, &
       e28a = 2.65, &
       e28b = 11.27, &
       e28c = 21.0, &
       e29a = 15.42, &
       e29b = 16.5, &
       e29c = 17.7, &
       e30a = 30.0, &
       e30b = 70.0

  do i = 1, nfreq

     nu(i) = 10.0**(dble(i-1)*dble(freqdel))

     lambda = clight/(nu(i)*eV_to_Hz)*1.e8 ! angstroms
     sigmaDust(i) = dustCrossSection(lambda/1.e4,1) * 1.e-22 ! cm^2

!     print*, '---', nu(i), sigmaDust(i)
!     print*, dustCrossSection(clight/(hydrogenIonization*eV_to_Hz)*1.e8/1.e4,1)
!     print*, dustCrossSection(0.5d0,1) ! 5000A, 10^-22 [cm^-2]

     IF ( nu(i) .GT. hydrogenIonization ) THEN

! i) HI photon-ionization cross section 

        dum = dsqrt(nu(i)/hydrogenIonization-1)
        sigma24(i) = 6.3e-18*(hydrogenIonization/nu(i))**4* &
             exp(4.0-4.0*atan(dum)/dum)/ (1-exp(-2.0*pi/dum))
     ELSE
        sigma24(i) = 0.0
     ENDIF

     IF ( nu(i) .GT. doubleHeliumIonization ) THEN

! iii) HeII photon-ionization cross section 

        dum = dsqrt(nu(i)/doubleHeliumIonization-1)
        sigma25(i) = 1.58e-18*(doubleHeliumIonization/nu(i))**4* &
             exp(4.0-4.0*atan(dum)/dum)/ (1-exp(-2.0*pi/dum))
     ELSE
        sigma25(i) = 0.0
     ENDIF

     IF ( nu(i) .GT. singleHeliumIonization ) THEN

! ii) HeI photon-ionization cross section 

        sigma26(i) = 7.42e-18*(1.66*(nu(i)/singleHeliumIonization)**(-2.05) - &
             0.66*(nu(i)/singleHeliumIonization)**(-3.05))
     ELSE
        sigma26(i) = 0.0
     ENDIF

     IF ( nu(i) .GT. e27 ) THEN
        sigma27(i) = 2.11e-16*(nu(i)-e27)**1.5/nu(i)**3
     ELSE
        sigma27(i) = 0.0
     ENDIF

     IF ( nu(i) .GT. e28a .AND. nu(i) .LE. e28b ) THEN
        sigma28(i) = 10**(-40.97+6.03*nu(i)-0.504*nu(i)**2+1.387e-2*nu(i)**3) 
     ELSEIF ( nu(i) .GT. e28b .AND. nu(i) .LT. e28c ) THEN
        sigma28(i) = 10**(-30.26+2.79*nu(i)-0.184*nu(i)**2+3.535e-3*nu(i)**3) 
     ELSE
        sigma28(i) = 0.0
     ENDIF

     IF ( nu(i) .GT. e29a .AND. nu(i) .LE. e29b ) THEN
        sigma29(i) = 6.2e-18*nu(i) - 9.4e-17
     ELSEIF ( nu(i) .GT. e29b .AND. nu(i) .LE. e29c ) THEN
        sigma29(i) = 1.4e-18*nu(i) - 1.48e-17
     ELSEIF ( nu(i) .GT. e29c ) THEN
        sigma29(i) = 2.5e-14*nu(i)**(-2.71)
     ELSE
        sigma29(i) = 0.0
     ENDIF

     IF ( nu(i) .GE. e30a .AND. nu(i) .LT. e30b ) THEN
        sigma30(i) = 10**(-16.926-4.528e-2*nu(i)+2.238e-4*nu(i)**2+4.245e-7*nu(i)**3)
     ELSE
        sigma30(i) = 0.0
     ENDIF

     IF ( nu(i) .GT. e28b .AND. nu(i) .LT. hydrogenIonization ) THEN
        sigma31(i) = 3.71e-18
     ELSE
        sigma31(i) = 0.0
     ENDIF
!       IF (H2_shield_flag) sigma31(i) = 0.0

!     write(*,1010) lambda, sigmaDust(i), sigma24(i), sigma25(i), sigma26(i)

1010 format(f6.0, 4es15.3)

  enddo

! output frequencies

  do ienergy = 1, nenergy
!     freq = lowerEnergy * exp((float(ienergy)-0.5)/float(nenergy)*(log(upperEnergy)-log(lowerEnergy)))

     freq = lowerEnergy * exp(float(ienergy-1)/float(nenergy-1)*(log(upperEnergy)-log(lowerEnergy)))

     outputFreq(ienergy) = freq
!     print*, ienergy, freq
     lambda = clight/(freq*eV_to_Hz)*1.e8 ! angstroms
     outputSigmaDust(ienergy) = dustCrossSection(lambda/1.e4,1) * 1.e-22 ! cm^2
     if (freq.gt.hydrogenIonization ) then
        dum = dsqrt(freq/hydrogenIonization-1)
        outputSigma24(ienergy) = 6.3e-18*(hydrogenIonization/freq)**4* &
             exp(4.-4.*atan(dum)/dum)/ (1-exp(-2.*pi/dum))
     else
        if (freq.eq.hydrogenIonization) then
           outputSigma24(ienergy) = 6.3e-18
        else
           outputSigma24(ienergy) = 0.
        endif
     endif
     if (freq.gt.doubleHeliumIonization ) then
        dum = dsqrt(freq/doubleHeliumIonization-1)
        outputSigma25(ienergy) = 1.58e-18*(doubleHeliumIonization/freq)**4* &
             exp(4.-4.*atan(dum)/dum)/ (1-exp(-2.*pi/dum))
     else
        outputSigma25(ienergy) = 0.
     endif
     if (freq.gt.singleHeliumIonization ) then
        outputSigma26(ienergy) = 7.42e-18*(1.66*(freq/singleHeliumIonization)**(-2.05) - &
             0.66*(freq/singleHeliumIonization)**(-3.05))
     else
        outputSigma26(ienergy) = 0.
     endif
  enddo








! ... initialize recombination and bremsstrahlung emissivities
!     call init_eta(ier)

  totalIntegral = 0.

  reactionRate1 = 0.
  reactionRate2 = 0.
  reactionRate3 = 0.
  energyRate1 = 0.
  energyRate2 = 0.
  energyRate3 = 0.

!   group1%beta24 = 0.d0
!   group1%beta25 = 0.d0
!   group1%beta26 = 0.d0
!   group1%beta27 = 0.d0
!   group1%beta28 = 0.d0
!   group1%beta29 = 0.d0
!   group1%beta30 = 0.d0
!   group1%beta31 = 0.d0
!   group1%gammaHI = 0.d0
!   group1%gammaHeI = 0.d0
!   group1%gammaHeII = 0.d0

!   group2%beta24 = 0.d0
!   group2%beta25 = 0.d0
!   group2%beta26 = 0.d0
!   group2%beta27 = 0.d0
!   group2%beta28 = 0.d0
!   group2%beta29 = 0.d0
!   group2%beta30 = 0.d0
!   group2%beta31 = 0.d0
!   group2%gammaHI = 0.d0
!   group2%gammaHeI = 0.d0
!   group2%gammaHeII = 0.d0

!   group3%beta24 = 0.d0
!   group3%beta25 = 0.d0
!   group3%beta26 = 0.d0
!   group3%beta27 = 0.d0
!   group3%beta28 = 0.d0
!   group3%beta29 = 0.d0
!   group3%beta30 = 0.d0
!   group3%beta31 = 0.d0
!   group3%gammaHI = 0.d0
!   group3%gammaHeI = 0.d0
!   group3%gammaHeII = 0.d0

!   denominatorForBeta1 = 0.
!   denominatorForBeta2 = 0.
!   denominatorForBeta3 = 0.

!   tmp1 = stellarPopulation(iSpectrum,coefSpectrum,nu1)/(nu1*eV_to_erg) ! [1/sec/Hz]
!   tmp2 = stellarPopulation(iSpectrum,coefSpectrum,nu2)/(nu2*eV_to_erg) ! [1/sec/Hz]
!   tmp3 = stellarPopulation(iSpectrum,coefSpectrum,nu3)/(nu3*eV_to_erg) ! [1/sec/Hz]

  do i = 2, nfreq

     freq = nu(i)
     delta_nu = nu(i) - nu(i-1)

! lookup spectrum in the stellar population model

     thisSpecificLuminosity = stellarPopulation(iSpectrum,coefSpectrum,iMetal,coefMetal,freq) ! [erg/s/Hz]

     dtmp = thisSpecificLuminosity/(freq*eV_to_erg) * delta_nu*eV_to_Hz ! [1/s] - tabulated spectrum
     if (freq.ge.nu1) totalIntegral = totalIntegral + dtmp ! [1/s]

     do idepth1 = 0, ndepth1
        do idepth2 = 0, ndepth2
           do idepth3 = 0, ndepth3
              do idepthDust = 0, ndepthDust

                 ! optical depths at each reaction's threshold

                 tau1 = float(idepth1)/float(ndepth1)*maxOpticalDepth1
                 tau2 = float(idepth2)/float(ndepth2)*maxOpticalDepth2
                 tau3 = float(idepth3)/float(ndepth3)*maxOpticalDepth3

                 ! dust optical depth at the Ly-limit (hydrogenIonization)

                 tauDust = float(idepthDust)/float(ndepthDust)*maxOpticalDepthDust

                 ! convert to frequency-dependent optical depths

                 tau1 = sigma24(i) / 6.3e-18 * tau1
                 tau2 = sigma26(i) / 7.42e-18 * tau2
                 tau3 = sigma25(i) / 1.58e-18 * tau3

                 tauDust = sigmaDust(i) / 5.4116737e-22 * tauDust

                 if (freq.ge.nu1) then
                    atmp = dtmp * exp(-(tau1+tau2+tau3+tauDust)) ! [1/s]
                    reactionRate1(idepth1,idepth2,idepth3,iDepthDust) = &
                         reactionRate1(idepth1,idepth2,idepth3,iDepthDust) + &
                         atmp ! [1/s]
                    energyRate1(idepth1,idepth2,idepth3,iDepthDust) = &
                         energyRate1(idepth1,idepth2,idepth3,iDepthDust) + &
                         (freq-nu1)*eV_to_erg * atmp ! [erg/s]
                 endif

                 if (freq.ge.nu2) then
                    atmp = dtmp * exp(-(tau1+tau2+tau3+tauDust))
                    reactionRate2(idepth1,idepth2,idepth3,iDepthDust) = &
                         reactionRate2(idepth1,idepth2,idepth3,iDepthDust) + &
                         atmp ! [1/s]
                    energyRate2(idepth1,idepth2,idepth3,iDepthDust) = &
                         energyRate2(idepth1,idepth2,idepth3,iDepthDust) + &
                         (freq-nu2)*eV_to_erg * atmp ! [erg/s]
                 endif

                 if (freq.ge.nu3) then
                    atmp = dtmp * exp(-(tau1+tau2+tau3+tauDust))
                    reactionRate3(idepth1,idepth2,idepth3,iDepthDust) = &
                         reactionRate3(idepth1,idepth2,idepth3,iDepthDust) + &
                         atmp ! [1/s]
                    energyRate3(idepth1,idepth2,idepth3,iDepthDust) = &
                         energyRate3(idepth1,idepth2,idepth3,iDepthDust) + &
                         (freq-nu3)*eV_to_erg * atmp ! [erg/s]
                 endif

              enddo
           enddo
        enddo
     enddo




!      if (freq.ge.nu1 .and. freq.le.nu2) then

! !        dtmp = freq**2/(exp(freq*eV_to_erg/(kb*temp))-1.) * delta_nu ! [eV*3] - thermal spectrum
!         dtmp = thisSpecificLuminosity/(freq*eV_to_erg) * delta_nu ! [eV/sec/Hz] - tabulated spectrum

!         totalIntegral = totalIntegral + dtmp ! [eV*3] or [eV/sec/Hz]

!         group1%beta24 = group1%beta24 + dtmp*sigma24(i) ! [cm^2 eV*3] or [cm^2 eV/sec/Hz]
!         group1%beta25 = group1%beta25 + dtmp*sigma25(i)
!         group1%beta26 = group1%beta26 + dtmp*sigma26(i)
!         group1%beta27 = group1%beta27 + dtmp*sigma27(i)
!         group1%beta28 = group1%beta28 + dtmp*sigma28(i)
!         group1%beta29 = group1%beta29 + dtmp*sigma29(i)
!         group1%beta30 = group1%beta30 + dtmp*sigma30(i)
!         group1%beta31 = group1%beta31 + dtmp*sigma31(i)

!         denominatorForBeta1 = denominatorForBeta1 + dtmp ! [eV*3] or [eV/sec/Hz]

!         group1%gammaHI = group1%gammaHI + dtmp*(freq-nu1)*eV_to_erg ! [erg eV*3] or [erg eV/sec/Hz]

!      endif

!      if (freq.ge.nu2 .and. freq.le.nu3) then

! !        dtmp = freq**2/(exp(freq*eV_to_erg/(kb*temp))-1.) * delta_nu ! [eV*3] - thermal spectrum
!         dtmp = thisSpecificLuminosity/(freq*eV_to_erg) * delta_nu ! [eV/sec/Hz] - tabulated spectrum

!         totalIntegral = totalIntegral + dtmp ! [eV*3] or [eV/sec/Hz]

!         group2%beta24 = group2%beta24 + dtmp*sigma24(i) ! [cm^2 eV*3] or [cm^2 eV/sec/Hz]
!         group2%beta25 = group2%beta25 + dtmp*sigma25(i)
!         group2%beta26 = group2%beta26 + dtmp*sigma26(i)
!         group2%beta27 = group2%beta27 + dtmp*sigma27(i)
!         group2%beta28 = group2%beta28 + dtmp*sigma28(i)
!         group2%beta29 = group2%beta29 + dtmp*sigma29(i)
!         group2%beta30 = group2%beta30 + dtmp*sigma30(i)
!         group2%beta31 = group2%beta31 + dtmp*sigma31(i)

!         denominatorForBeta2 = denominatorForBeta2 + dtmp ! [eV*3] or [eV/sec/Hz]

!         group2%gammaHI = group2%gammaHI + dtmp*(freq-nu1)*eV_to_erg ! [erg eV*3] or [erg eV/sec/Hz]
!         group2%gammaHeI = group2%gammaHeI + dtmp*(freq-nu2)*eV_to_erg

!      endif

!      if (freq.ge.nu3) then

! !        dtmp = freq**2/(exp(freq*eV_to_erg/(kb*temp))-1.) * delta_nu ! [eV*3] - thermal spectrum
!         dtmp = thisSpecificLuminosity/(freq*eV_to_erg) * delta_nu ! [eV/sec/Hz] - tabulated spectrum

!         totalIntegral = totalIntegral + dtmp ! [eV*3] or [eV/sec/Hz]

!         group3%beta24 = group3%beta24 + dtmp*sigma24(i) ! [cm^2 eV*3] or [cm^2 eV/sec/Hz]
!         group3%beta25 = group3%beta25 + dtmp*sigma25(i)
!         group3%beta26 = group3%beta26 + dtmp*sigma26(i)
!         group3%beta27 = group3%beta27 + dtmp*sigma27(i)
!         group3%beta28 = group3%beta28 + dtmp*sigma28(i)
!         group3%beta29 = group3%beta29 + dtmp*sigma29(i)
!         group3%beta30 = group3%beta30 + dtmp*sigma30(i)
!         group3%beta31 = group3%beta31 + dtmp*sigma31(i)

!         denominatorForBeta3 = denominatorForBeta3 + dtmp ! [eV*3] or [eV/sec/Hz]

!         group3%gammaHI = group3%gammaHI + dtmp*(freq-nu1)*eV_to_erg ! [erg eV*3] or [erg eV/sec/Hz]
!         group3%gammaHeI = group3%gammaHeI + dtmp*(freq-nu2)*eV_to_erg
!         group3%gammaHeII = group3%gammaHeII + dtmp*(freq-nu3)*eV_to_erg

!      endif

  enddo

!   group1%beta24 = group1%beta24 / denominatorForBeta1 ! [cm^2]
!   group1%beta25 = group1%beta25 / denominatorForBeta1
!   group1%beta26 = group1%beta26 / denominatorForBeta1
!   group1%beta27 = group1%beta27 / denominatorForBeta1
!   group1%beta28 = group1%beta28 / denominatorForBeta1
!   group1%beta29 = group1%beta29 / denominatorForBeta1
!   group1%beta30 = group1%beta30 / denominatorForBeta1
!   group1%beta31 = group1%beta31 / denominatorForBeta1
! !  dtmp = (exp(nu1*eV_to_erg/(kb*temp))-1.) / nu1**2 * eV_to_Hz
! !  group1%ksi = denominatorForBeta1 * dtmp ! [Hz]
! !  group1%gammaHI = group1%gammaHI * dtmp ! [erg Hz]
!   group1%ksi = denominatorForBeta1 / tmp1 * eV_to_Hz ! [Hz]
!   group1%gammaHI = group1%gammaHI / tmp1 * eV_to_Hz ! [erg Hz]

!   group2%beta24 = group2%beta24 / denominatorForBeta2 ! [cm^2]
!   group2%beta25 = group2%beta25 / denominatorForBeta2
!   group2%beta26 = group2%beta26 / denominatorForBeta2
!   group2%beta27 = group2%beta27 / denominatorForBeta2
!   group2%beta28 = group2%beta28 / denominatorForBeta2
!   group2%beta29 = group2%beta29 / denominatorForBeta2
!   group2%beta30 = group2%beta30 / denominatorForBeta2
!   group2%beta31 = group2%beta31 / denominatorForBeta2
! !  dtmp = (exp(nu2*eV_to_erg/(kb*temp))-1.) / nu2**2 * eV_to_Hz
! !  group2%ksi = denominatorForBeta2 * dtmp ! [Hz]
! !  group2%gammaHI = group2%gammaHI * dtmp ! [erg Hz]
! !  group2%gammaHeI = group2%gammaHeI * dtmp
!   group2%ksi = denominatorForBeta2 / tmp2 * eV_to_Hz ! [Hz]
!   group2%gammaHI = group2%gammaHI / tmp2 * eV_to_Hz ! [erg Hz]
!   group2%gammaHeI = group2%gammaHeI / tmp2 * eV_to_Hz ! [erg Hz]

!   group3%beta24 = group3%beta24 / denominatorForBeta3 ! [cm^2]
!   group3%beta25 = group3%beta25 / denominatorForBeta3
!   group3%beta26 = group3%beta26 / denominatorForBeta3
!   group3%beta27 = group3%beta27 / denominatorForBeta3
!   group3%beta28 = group3%beta28 / denominatorForBeta3
!   group3%beta29 = group3%beta29 / denominatorForBeta3
!   group3%beta30 = group3%beta30 / denominatorForBeta3
!   group3%beta31 = group3%beta31 / denominatorForBeta3
! !  dtmp = (exp(nu3*eV_to_erg/(kb*temp))-1.) / nu3**2 * eV_to_Hz
! !  group3%ksi = denominatorForBeta3 * dtmp ! [Hz]
! !  group3%gammaHI = group3%gammaHI * dtmp ! [erg Hz]
! !  group3%gammaHeI = group3%gammaHeI * dtmp
! !  group3%gammaHeII = group3%gammaHeII * dtmp
!   group3%ksi = denominatorForBeta3 / tmp3 * eV_to_Hz ! [Hz]
!   group3%gammaHI = group3%gammaHI / tmp3 * eV_to_Hz ! [erg Hz]
!   group3%gammaHeI = group3%gammaHeI / tmp3 * eV_to_Hz ! [erg Hz]
!   group3%gammaHeII = group3%gammaHeII / tmp3 * eV_to_Hz ! [erg Hz]

!   write(*,*) group1%beta24, group1%beta25, group1%beta26
!   write(*,*) group1%beta27, group1%beta28, group1%beta29
!   write(*,*) group1%beta30, group1%beta31
!   write(*,*) 'g =', group1%gammaHI, group1%gammaHeI, group1%gammaHeII
!   write(*,*) 'g =', group2%gammaHI, group2%gammaHeI, group2%gammaHeII
!   write(*,*) 'g =', group3%gammaHI, group3%gammaHeI, group3%gammaHeII
!   stop

  return

end subroutine stellarBetaTable
