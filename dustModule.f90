module dust

  use definitions

  implicit none
  real(kind=RealKind), parameter :: extinctionToTauCorrection = 0.9210340372 ! =ln(10)/2.5

contains

  subroutine dustInitialize

    implicit none
    integer i

    open(14,file='smc_dust_parameters.dat',status='old')
    open(15,file='lmc_dust_parameters.dat',status='old')

    do i = 1, 7
       read(14,*)a_smc(i,:)
    enddo

    do i = 1, 7
       read(15,*)a_lmc(i,:)
    enddo

    return

  end subroutine dustInitialize

  function dustCrossSection(lambda,idust) result(crossSection)

    implicit none

    integer, intent(in) :: idust
    real(kind=RealKind), intent(in) :: lambda
    real(kind=RealKind) :: crossSection, sigma, x
    integer :: i

    if (idust==1) then !SMC dust
      
       sigma=0
      
       do i=1,7
          x=lambda/a_smc(i,1)
          sigma=sigma + a_smc(i,2)/(x**a_smc(i,4)+x**(-a_smc(i,5))+a_smc(i,3))
       enddo

       crossSection = 1.1*sigma*extinctionToTauCorrection

       return

    elseif(idust==2) then

       sigma = 0
      
       do i = 1, 7
          x = lambda/a_lmc(i,1)
          sigma = sigma + a_lmc(i,2)/(x**a_lmc(i,4)+x**(-a_lmc(i,5))+a_lmc(i,3))
       enddo

       crossSection = 3.3*sigma*extinctionToTauCorrection

       return

    else

       print*,'unknown dust-type'

       stop

    endif

  end function dustCrossSection

end module dust
