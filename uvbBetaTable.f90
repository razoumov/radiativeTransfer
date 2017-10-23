subroutine uvbBetaTable(nfreq,freqdel,alpha)

  !--- Alexei Razoumov, 2005-
  !--- Supplementary file for the Fully Threaded Transport Engine.

  use definitions

  implicit none
  integer, intent(in) :: nfreq
  real(kind=RealKind), intent(in) :: freqdel
  real(kind=RealKind), dimension(3), intent(in) :: alpha

  integer :: i, irate, ier
  real(kind=RealKind) :: shape1, shape2, shape3, &
       energyShape1, energyShape2, energyShape3, &
       delta_nu, dtmp, dtmpOverEnergy, freq, dum
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

  group1%beta24 = 0.d0
  group1%beta25 = 0.d0
  group1%beta26 = 0.d0
  group1%beta27 = 0.d0
  group1%beta28 = 0.d0
  group1%beta29 = 0.d0
  group1%beta30 = 0.d0
  group1%beta31 = 0.d0
  group1%ksi24 = 0.d0
  group1%ksi25 = 0.d0
  group1%ksi26 = 0.d0
  group1%ksi27 = 0.d0
  group1%ksi28 = 0.d0
  group1%ksi29 = 0.d0
  group1%ksi30 = 0.d0
  group1%ksi31 = 0.d0
  group1%gammaHI = 0.d0
  group1%gammaHeI = 0.d0
  group1%gammaHeII = 0.d0

  group2%beta24 = 0.d0
  group2%beta25 = 0.d0
  group2%beta26 = 0.d0
  group2%beta27 = 0.d0
  group2%beta28 = 0.d0
  group2%beta29 = 0.d0
  group2%beta30 = 0.d0
  group2%beta31 = 0.d0
  group2%ksi24 = 0.d0
  group2%ksi25 = 0.d0
  group2%ksi26 = 0.d0
  group2%ksi27 = 0.d0
  group2%ksi28 = 0.d0
  group2%ksi29 = 0.d0
  group2%ksi30 = 0.d0
  group2%ksi31 = 0.d0
  group2%gammaHI = 0.d0
  group2%gammaHeI = 0.d0
  group2%gammaHeII = 0.d0

  group3%beta24 = 0.d0
  group3%beta25 = 0.d0
  group3%beta26 = 0.d0
  group3%beta27 = 0.d0
  group3%beta28 = 0.d0
  group3%beta29 = 0.d0
  group3%beta30 = 0.d0
  group3%beta31 = 0.d0
  group3%ksi24 = 0.d0
  group3%ksi25 = 0.d0
  group3%ksi26 = 0.d0
  group3%ksi27 = 0.d0
  group3%ksi28 = 0.d0
  group3%ksi29 = 0.d0
  group3%ksi30 = 0.d0
  group3%ksi31 = 0.d0
  group3%gammaHI = 0.d0
  group3%gammaHeI = 0.d0
  group3%gammaHeII = 0.d0

  do i = 2, nfreq

     freq = nu(i)
     delta_nu = nu(i) - nu(i-1)

     if (freq.ge.nu1 .and. freq.le.nu2) then

        dtmp = (freq/nu1)**(-alpha(1))*delta_nu
        dtmpOverEnergy = dtmp*eV_to_Hz/(freq*eV_to_erg)

        group1%beta24 = group1%beta24 + dtmp*sigma24(i)
        group1%beta25 = group1%beta25 + dtmp*sigma25(i)
        group1%beta26 = group1%beta26 + dtmp*sigma26(i)
        group1%beta27 = group1%beta27 + dtmp*sigma27(i)
        group1%beta28 = group1%beta28 + dtmp*sigma28(i)
        group1%beta29 = group1%beta29 + dtmp*sigma29(i)
        group1%beta30 = group1%beta30 + dtmp*sigma30(i)
        group1%beta31 = group1%beta31 + dtmp*sigma31(i)

        group1%ksi24 = group1%ksi24 + dtmpOverEnergy*sigma24(i) ! [cm^2 Hz / erg]
        group1%ksi25 = group1%ksi25 + dtmpOverEnergy*sigma25(i)
        group1%ksi26 = group1%ksi26 + dtmpOverEnergy*sigma26(i)
        group1%ksi27 = group1%ksi27 + dtmpOverEnergy*sigma27(i)
        group1%ksi28 = group1%ksi28 + dtmpOverEnergy*sigma28(i)
        group1%ksi29 = group1%ksi29 + dtmpOverEnergy*sigma29(i)
        group1%ksi30 = group1%ksi30 + dtmpOverEnergy*sigma30(i)
        group1%ksi31 = group1%ksi31 + dtmpOverEnergy*sigma31(i)

        group1%gammaHI = group1%gammaHI + dtmpOverEnergy*(freq-nu1)*eV_to_erg*sigma24(i) ! [cm^2 Hz]

     endif

     if (freq.ge.nu2 .and. freq.le.nu3) then

        dtmp = (freq/nu2)**(-alpha(2))*delta_nu
        dtmpOverEnergy = dtmp*eV_to_Hz/(freq*eV_to_erg)

        group2%beta24 = group2%beta24 + dtmp*sigma24(i)
        group2%beta25 = group2%beta25 + dtmp*sigma25(i)
        group2%beta26 = group2%beta26 + dtmp*sigma26(i)
        group2%beta27 = group2%beta27 + dtmp*sigma27(i)
        group2%beta28 = group2%beta28 + dtmp*sigma28(i)
        group2%beta29 = group2%beta29 + dtmp*sigma29(i)
        group2%beta30 = group2%beta30 + dtmp*sigma30(i)
        group2%beta31 = group2%beta31 + dtmp*sigma31(i)

        group2%ksi24 = group2%ksi24 + dtmpOverEnergy*sigma24(i) ! [cm^2 Hz / erg]
        group2%ksi25 = group2%ksi25 + dtmpOverEnergy*sigma25(i)
        group2%ksi26 = group2%ksi26 + dtmpOverEnergy*sigma26(i)
        group2%ksi27 = group2%ksi27 + dtmpOverEnergy*sigma27(i)
        group2%ksi28 = group2%ksi28 + dtmpOverEnergy*sigma28(i)
        group2%ksi29 = group2%ksi29 + dtmpOverEnergy*sigma29(i)
        group2%ksi30 = group2%ksi30 + dtmpOverEnergy*sigma30(i)
        group2%ksi31 = group2%ksi31 + dtmpOverEnergy*sigma31(i)

        group2%gammaHI = group2%gammaHI + dtmpOverEnergy*(freq-nu1)*eV_to_erg*sigma24(i) ! [cm^2 Hz]
        group2%gammaHeI = group2%gammaHeI + dtmpOverEnergy*(freq-nu2)*eV_to_erg*sigma26(i)

     endif

     if (freq.ge.nu3) then

        dtmp = (freq/nu3)**(-alpha(3))*delta_nu
        dtmpOverEnergy = dtmp*eV_to_Hz/(freq*eV_to_erg)

        group3%beta24 = group3%beta24 + dtmp*sigma24(i)
        group3%beta25 = group3%beta25 + dtmp*sigma25(i)
        group3%beta26 = group3%beta26 + dtmp*sigma26(i)
        group3%beta27 = group3%beta27 + dtmp*sigma27(i)
        group3%beta28 = group3%beta28 + dtmp*sigma28(i)
        group3%beta29 = group3%beta29 + dtmp*sigma29(i)
        group3%beta30 = group3%beta30 + dtmp*sigma30(i)
        group3%beta31 = group3%beta31 + dtmp*sigma31(i)

        group3%ksi24 = group3%ksi24 + dtmpOverEnergy*sigma24(i) ! [cm^2 Hz / erg]
        group3%ksi25 = group3%ksi25 + dtmpOverEnergy*sigma25(i)
        group3%ksi26 = group3%ksi26 + dtmpOverEnergy*sigma26(i)
        group3%ksi27 = group3%ksi27 + dtmpOverEnergy*sigma27(i)
        group3%ksi28 = group3%ksi28 + dtmpOverEnergy*sigma28(i)
        group3%ksi29 = group3%ksi29 + dtmpOverEnergy*sigma29(i)
        group3%ksi30 = group3%ksi30 + dtmpOverEnergy*sigma30(i)
        group3%ksi31 = group3%ksi31 + dtmpOverEnergy*sigma31(i)

        group3%gammaHI = group3%gammaHI + dtmpOverEnergy*(freq-nu1)*eV_to_erg*sigma24(i) ! [cm^2 Hz]
        group3%gammaHeI = group3%gammaHeI + dtmpOverEnergy*(freq-nu2)*eV_to_erg*sigma26(i)
        group3%gammaHeII = group3%gammaHeII + dtmpOverEnergy*(freq-nu3)*eV_to_erg*sigma25(i)

     endif

  enddo

  shape1 = (1.-(nu2/nu1)**(1.-alpha(1)))/(alpha(1)-1.)
  shape2 = (1.-(nu3/nu2)**(1.-alpha(2)))/(alpha(2)-1.)
  shape3 = 1./(alpha(3)-1.)

  print*, 'shape =', shape1, shape2, shape3

  energyShape1 = shape1*nu1
  group1%beta24 = group1%beta24 / energyShape1 ! [cm^2]
  group1%beta25 = group1%beta25 / energyShape1
  group1%beta26 = group1%beta26 / energyShape1
  group1%beta27 = group1%beta27 / energyShape1
  group1%beta28 = group1%beta28 / energyShape1
  group1%beta29 = group1%beta29 / energyShape1
  group1%beta30 = group1%beta30 / energyShape1
  group1%beta31 = group1%beta31 / energyShape1

  energyShape2 = shape2*nu2
  group2%beta24 = group2%beta24 / energyShape2 ! [cm^2]
  group2%beta25 = group2%beta25 / energyShape2
  group2%beta26 = group2%beta26 / energyShape2
  group2%beta27 = group2%beta27 / energyShape2
  group2%beta28 = group2%beta28 / energyShape2
  group2%beta29 = group2%beta29 / energyShape2
  group2%beta30 = group2%beta30 / energyShape2
  group2%beta31 = group2%beta31 / energyShape2

  energyShape3 = shape3*nu3
  group3%beta24 = group3%beta24 / energyShape3 ! [cm^2]
  group3%beta25 = group3%beta25 / energyShape3
  group3%beta26 = group3%beta26 / energyShape3
  group3%beta27 = group3%beta27 / energyShape3
  group3%beta28 = group3%beta28 / energyShape3
  group3%beta29 = group3%beta29 / energyShape3
  group3%beta30 = group3%beta30 / energyShape3
  group3%beta31 = group3%beta31 / energyShape3

!   write(*,*) group1%beta24, group1%beta25, group1%beta26
!   write(*,*) group1%beta27, group1%beta28, group1%beta29
!   write(*,*) group1%beta30, group1%beta31
!   write(*,*) 'g =', group1%gammaHI, group1%gammaHeI, group1%gammaHeII
!   write(*,*) 'g =', group2%gammaHI, group2%gammaHeI, group2%gammaHeII
!   write(*,*) 'g =', group3%gammaHI, group3%gammaHeI, group3%gammaHeII

  return

end subroutine uvbBetaTable
