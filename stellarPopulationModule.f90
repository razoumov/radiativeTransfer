module stellarPopulationModule

  use definitions

contains

  function stellarPopulation(iSpectrum,coefSpectrum,iMetal,coefMetal,freq) result(SP)

    implicit none
    integer, intent(in) :: iSpectrum, iMetal
    real(kind=RealKind), intent(in) :: coefSpectrum, coefMetal, freq
    integer :: iWavelength
    real(kind=RealKind) :: thisWavelength, coefWavelength, SP, sp1, sp2

    thisWavelength = clight/(freq*eV_to_Hz)
    iWavelength = 1
    do while (thisWavelength.gt.wavelength(iWavelength+1))
       iWavelength = iWavelength + 1
    enddo

    coefWavelength = (thisWavelength-wavelength(iWavelength))/ &
         (wavelength(iWavelength+1)-wavelength(iWavelength))
    coefWavelength = min(max(0.d0,coefWavelength),1.d0)

!    print*, coefMetal, coefSpectrum, coefWavelength

    sp1 = coefSpectrum * &
         ((1.-coefWavelength)*specificLuminosity(iMetal,iSpectrum+1,iWavelength) + &
         coefWavelength*1.*specificLuminosity(iMetal,iSpectrum+1,iWavelength+1)) +  &
         (1.-coefSpectrum) * &
         ((1.-coefWavelength)*specificLuminosity(iMetal,iSpectrum,iWavelength) + &
         coefWavelength*specificLuminosity(iMetal,iSpectrum,iWavelength+1)) ! [log(erg/sec/A)]
    sp2 = coefSpectrum * &
         ((1.-coefWavelength)*specificLuminosity(iMetal+1,iSpectrum+1,iWavelength) + &
         coefWavelength*1.*specificLuminosity(iMetal+1,iSpectrum+1,iWavelength+1)) +  &
         (1.-coefSpectrum) * &
         ((1.-coefWavelength)*specificLuminosity(iMetal+1,iSpectrum,iWavelength) + &
         coefWavelength*specificLuminosity(iMetal+1,iSpectrum,iWavelength+1)) ! [log(erg/sec/A)]
    SP = (1.-coefMetal)*sp1 + coefMetal*sp2

!     SP = coefSpectrum * &
!          ((1.-coefWavelength)*specificLuminosity(iMetal,iSpectrum+1,iWavelength) + &
!          coefWavelength*1.*specificLuminosity(iMetal,iSpectrum+1,iWavelength+1)) +  &
!          (1.-coefSpectrum) * &
!          ((1.-coefWavelength)*specificLuminosity(iMetal,iSpectrum,iWavelength) + &
!          coefWavelength*specificLuminosity(iMetal,iSpectrum,iWavelength+1)) ! [log(erg/sec/A)]

    SP = (10.**SP)/angstrom * clight/(freq*eV_to_Hz)**2 ! [erg/sec/Hz]

  end function stellarPopulation

end module stellarPopulationModule
