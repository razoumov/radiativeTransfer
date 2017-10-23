subroutine uniformTable(nfreq,freqdel,alphaQuasar,alphaStellar)

  !--- Alexei Razoumov, 2005-
  !--- Supplementary file for the Fully Threaded Transport Engine.

  use definitions

  implicit none
  integer, intent(in) :: nfreq
  real(kind=RealKind), intent(in) :: freqdel, alphaQuasar, alphaStellar

  integer :: i, irate, ier
  real(kind=RealKind) :: delta_nu, dtmp, dtmpOverEnergy, freq, dum
  real(kind=RealKind), dimension(nfreq) :: sigma24, sigma25, sigma26, &
       sigma27, sigma28, sigma29, sigma30, sigma31, nu

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

!        write(*,*) nu(i), sigma24(i), sigma25(i), sigma26(i)

  enddo






! ... initialize recombination and bremsstrahlung emissivities
!     call init_eta(ier)

  quasar%ksi24 = 0.d0
  quasar%ksi25 = 0.d0
  quasar%ksi26 = 0.d0
  quasar%ksi27 = 0.d0
  quasar%ksi28 = 0.d0
  quasar%ksi29 = 0.d0
  quasar%ksi30 = 0.d0
  quasar%ksi31 = 0.d0
  quasar%gammaHI = 0.d0
  quasar%gammaHeI = 0.d0
  quasar%gammaHeII = 0.d0

  stellar%ksi24 = 0.d0
  stellar%ksi25 = 0.d0
  stellar%ksi26 = 0.d0
  stellar%ksi27 = 0.d0
  stellar%ksi28 = 0.d0
  stellar%ksi29 = 0.d0
  stellar%ksi30 = 0.d0
  stellar%ksi31 = 0.d0
  stellar%gammaHI = 0.d0
  stellar%gammaHeI = 0.d0
  stellar%gammaHeII = 0.d0

  do i = 2, nfreq

     freq = nu(i)
     delta_nu = nu(i) - nu(i-1)

!    rates from quasar component

     dtmp = (freq/nu1)**(-alphaQuasar) * delta_nu
     dtmpOverEnergy = dtmp*eV_to_Hz/(freq*eV_to_erg)

     if (freq.ge.nu1) then
        quasar%ksi24 = quasar%ksi24 + dtmpOverEnergy*sigma24(i) ! [cm^2 Hz / erg]
        quasar%ksi25 = quasar%ksi25 + dtmpOverEnergy*sigma25(i)
        quasar%ksi26 = quasar%ksi26 + dtmpOverEnergy*sigma26(i)
        quasar%ksi27 = quasar%ksi27 + dtmpOverEnergy*sigma27(i)
        quasar%ksi28 = quasar%ksi28 + dtmpOverEnergy*sigma28(i)
        quasar%ksi29 = quasar%ksi29 + dtmpOverEnergy*sigma29(i)
        quasar%ksi30 = quasar%ksi30 + dtmpOverEnergy*sigma30(i)
        quasar%ksi31 = quasar%ksi31 + dtmpOverEnergy*sigma31(i)
        quasar%gammaHI = quasar%gammaHI + dtmpOverEnergy*(freq-nu1)*eV_to_erg*sigma24(i) ! [cm^2 Hz]
     endif

     if (freq.ge.nu2) then
        quasar%gammaHeI = quasar%gammaHeI + dtmpOverEnergy*(freq-nu2)*eV_to_erg*sigma26(i)
     endif

     if (freq.ge.nu3) then
        quasar%gammaHeII = quasar%gammaHeII + dtmpOverEnergy*(freq-nu3)*eV_to_erg*sigma25(i)
     endif

!    rates from stellar component

     dtmp = (freq/nu1)**(-alphaStellar) * delta_nu
     dtmpOverEnergy = dtmp*eV_to_Hz/(freq*eV_to_erg)

     if (freq.ge.nu1) then
        stellar%ksi24 = stellar%ksi24 + dtmpOverEnergy*sigma24(i) ! [cm^2 Hz / erg]
        stellar%ksi25 = stellar%ksi25 + dtmpOverEnergy*sigma25(i)
        stellar%ksi26 = stellar%ksi26 + dtmpOverEnergy*sigma26(i)
        stellar%ksi27 = stellar%ksi27 + dtmpOverEnergy*sigma27(i)
        stellar%ksi28 = stellar%ksi28 + dtmpOverEnergy*sigma28(i)
        stellar%ksi29 = stellar%ksi29 + dtmpOverEnergy*sigma29(i)
        stellar%ksi30 = stellar%ksi30 + dtmpOverEnergy*sigma30(i)
        stellar%ksi31 = stellar%ksi31 + dtmpOverEnergy*sigma31(i)
        stellar%gammaHI = stellar%gammaHI + dtmpOverEnergy*(freq-nu1)*eV_to_erg*sigma24(i) ! [cm^2 Hz]
     endif

     if (freq.ge.nu2) then
        stellar%gammaHeI = stellar%gammaHeI + dtmpOverEnergy*(freq-nu2)*eV_to_erg*sigma26(i)
     endif

     if (freq.ge.nu3) then
        stellar%gammaHeII = stellar%gammaHeII + dtmpOverEnergy*(freq-nu3)*eV_to_erg*sigma25(i)
     endif

  enddo

!  write(*,*) quasar%gammaHI, quasar%gammaHeI, quasar%gammaHeII
!  write(*,*) stellar%gammaHI, stellar%gammaHeI, stellar%gammaHeII
!  stop

  return

end subroutine uniformTable
