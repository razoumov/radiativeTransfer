module definitions

  implicit none
  integer, parameter :: SingleKind = kind(1.0)
  integer, parameter :: RealKind = kind(1.0d0)
  integer*8 :: LongIntegerExample
  integer, parameter :: LongInteger = kind(LongIntegerExample)
  real(kind=RealKind), parameter :: pi = 3.141592654
  real(kind=RealKind), parameter :: halfPi = 0.5*pi
  real(kind=RealKind), parameter :: twoPi = 2.*pi
  real(kind=RealKind), parameter :: fourPi = 4.*pi
  real(kind=RealKind), parameter :: quarterPi = 0.25*pi
  real(kind=RealKind), parameter :: fortyFiveDegrees = pi/4.
  real(kind=RealKind), parameter :: ninetyDegrees = pi/2.
  real(kind=RealKind), parameter :: hp = 6.6260693e-27
  real(kind=RealKind), parameter :: kb = 1.3806503e-16
  real(kind=RealKind), parameter :: clight = 2.99792458e10
  real(kind=RealKind), parameter :: yr = 31557600
  real(kind=RealKind), parameter :: kyr = 1.e3*yr
  real(kind=RealKind), parameter :: Myr = 1.e6*yr
  real(kind=RealKind), parameter :: pc = 3.08568025e18
  real(kind=RealKind), parameter :: kpc = 1.e3*pc
  real(kind=RealKind), parameter :: Mpc = 1.e6*pc
  real(kind=RealKind), parameter :: angstrom = 1.e-8
  real(kind=RealKind), parameter :: mp = 1.6726231e-24
  real(kind=RealKind), parameter :: mn = 1.67492728e-24
  real(kind=RealKind), parameter :: mh = mp
  real(kind=RealKind), parameter :: mhe = 2.*(mp+mn)
  real(kind=RealKind), parameter :: msun = 1.98892e33
  real(kind=RealKind), parameter :: hydrogenIonization = 13.598
  real(kind=RealKind), parameter :: singleHeliumIonization = 24.587
  real(kind=RealKind), parameter :: doubleHeliumIonization = 54.418
  real(kind=RealKind), parameter :: nu1 = hydrogenIonization
  real(kind=RealKind), parameter :: nu2 = singleHeliumIonization
  real(kind=RealKind), parameter :: nu3 = doubleHeliumIonization
  real(kind=RealKind), parameter :: eV_to_erg = 1.60217646d-12
  real(kind=RealKind), parameter :: eV = 1.60217646d-12
  real(kind=RealKind), parameter :: eV_to_Hz = eV_to_erg / hp
  real(kind=RealKind), parameter :: gamma = 1.6667
  real(kind=RealKind), parameter :: nu_alpha = 2.466e15
  integer, parameter :: nAngularLevel = 3 ! (1 gives 12 angles, 2 gives 48 angles, 3 gives 192, 4 gives 768, 5 gives 3072, and so on)

  real(kind=RealKind), dimension(7,5) :: a_smc, a_lmc

  real(kind=RealKind) :: compa, piHI, piHeI, piHeII, currentRedshift

  integer, parameter :: caseA = 1, caseB = 2
  integer, parameter :: recombinationType = caseB

  integer, parameter :: ksi = 3 ! defines accuracy/speed of the tree walk

  integer, parameter :: mapHIColumnDensity = 1, mapTemperature = 2, checkStarHostCells = 3, &
       mapTotalGasDensity = 4, mapMaxGasDensity = 5, mapMetals = 6

  real(kind=RealKind) :: uvbStellar1, uvbStellar2, uvbStellar3, &
       uvbQuasar1, uvbQuasar2, uvbQuasar3, uvb1, uvb2, uvb3, &
       stellar99, quasar99, pascal02, component1, component2, step, &
       quasar02, stellar02, gaussian, newQuasar02, newStellar02, total02, &
       tmplambda, density_units, vel_units, &
       baseDensity, cosmicJmean1, fcosmic, meanHostDensity, selfShieldingThreshold, &
       totalMass, volSum, volSumDensity, volSumDensity2, totalPath, &
       uniformQuasar, uniformStellar

  integer, parameter :: npdf = 50
  real(kind=RealKind), parameter :: apdf = - 8., bpdf = 3.
  integer, dimension(npdf) :: pdfStar
  real(kind=RealKind), dimension(npdf) :: pdfGas
  integer :: pdfStarOutside
  real(kind=RealKind) :: pdfGasOutside

  ! hydrogen/helium table
  integer, parameter :: ndepth1 = 10, ndepth2 = 10, ndepth3 = 10, ndepthDust = 10
  real(kind=RealKind), parameter :: maxOpticalDepth1 = 10., &
       maxOpticalDepth2 = 10., maxOpticalDepth3 = 10., maxOpticalDepthDust = 10.
  real(kind=RealKind), dimension(0:ndepth1,0:ndepth2,0:ndepth3,0:ndepthDust) :: &
       reactionRate1, reactionRate2, reactionRate3, &
       energyRate1, energyRate2, energyRate3

  integer :: icosmic, ncosmic, nSources, nNodes, dustApproximation
  real(kind=RealKind) :: neutralHydrogenMass, totalHydrogenMass, finalRadius, densityCoefficient, &
       sourceTotalHydrogenDensity, ddd, totalVolume, oxygenMass, gasMass, &
       minOxygenMass, maxOxygenMass
  integer, parameter :: nabun = 1000
  real(kind=RealKind), parameter :: minAbun = -10., maxAbun = 2.
  real(kind=RealKind), dimension(nabun) :: abunArray, massAbove
  logical, parameter :: expansionFlag = .false.
  integer, parameter :: noDust = 0, completeSublimation = 1, noSublimation = 2
  logical :: readMetals, readKinematics
  integer :: massStellarParticle
  integer, parameter :: normal = 1, hiRes = 2, superHiRes = 3, massive = 10, hiResHeavy = 4, &
       crazyHiRes = 5, light = 6, lyAlpha = 7
  integer, parameter :: maxNumberReadVariables = 10

  type :: normCrossSectionType
     real(kind=RealKind) :: &
          ! group-averaged absorption cross-sections [cm^2]:
          beta24, beta25, beta26, beta27, beta28, beta29, beta30, beta31, &
          ! spectrum-averaged photo-reaction rates [cm^2 Hz / erg]:
          ksi24, ksi25, ksi26, ksi27, ksi28, ksi29, ksi30, ksi31, &
          ! spectrum-averaged photo-heating rates [cm^2 Hz]:
          gammaHI, gammaHeI, gammaHeII
  end type normCrossSectionType

  type(normCrossSectionType) :: quasar, stellar, group1, group2, group3

  type :: coolingRatesType
     ! photo-ionization cooling/heating (no temperature dependence)
     real(kind=RealKind) :: gammaHI, gammaHeI, gammaHeII
  end type coolingRatesType

  type :: xyRayType
     real(kind=RealKind) :: Iout1, Iout2, Iout3
  end type xyRayType

  type :: xzRayType
     real(kind=RealKind) :: Iout1, Iout2, Iout3
  end type xzRayType

  type :: yzRayType
     real(kind=RealKind) :: Iout1, Iout2, Iout3
  end type yzRayType

  type :: xyRayPatternType
     real(kind=RealKind) :: x0, y0, len
  end type xyRayPatternType

  type :: xzRayPatternType
     real(kind=RealKind) :: x0, z0, len
  end type xzRayPatternType

  type :: yzRayPatternType
     real(kind=RealKind) :: y0, z0, len
  end type yzRayPatternType

  type :: radiativeTransportType
     type(xyRayType) :: xyRay
     type(xzRayType) :: xzRay
     type(yzRayType) :: yzRay
  end type radiativeTransportType

  type :: patternType
     type(xyRayPatternType) :: xyRay
     type(xzRayPatternType) :: xzRay
     type(yzRayPatternType) :: yzRay
     logical(1) :: xzRayActive, yzRayActive
     logical(1) :: refined
     integer(1) :: level ! maximum 127 levels of refinement
     integer(1) :: xyTop ! xyTop shows which ray hits xy-plane: 1=xy, 2=yz, and 3=xz
     integer(1) :: xzTop ! xzTop shows which ray hits xz-plane: 1=xy, 2=yz, and 3=xz
     integer(1) :: yzTop ! yzTop shows which ray hits yz-plane: 1=xy, 2=yz, and 3=xz
     type(patternType), dimension(:,:,:), pointer :: cell
  end type patternType

  type :: angleType
     real(kind=RealKind) :: phi, theta
     integer(1) :: izone
  end type angleType

  integer, parameter :: xyEnd = 1, yzEnd = 2, xzEnd = 3
  integer, parameter :: xy = 1, yz = 2, xz = 3
  integer, parameter :: up = 1, down = -1

  type :: zoneType
     type(radiativeTransportType) :: rt
     real(kind=RealKind) :: rho, tgas, HI, HeI, HeII, &
          krate24, krate25, krate26, crate24, crate25, crate26, &
          HIprev, HeIprev, HeIIprev, hydroHeating, rhoCoef, &
          velx, vely, velz, abun2
     real(kind=RealKind) :: kappa1, kappa2, kappa3
! tgasprev, abun1, abun3, abun4
     real(kind=RealKind) :: Jmean1, Jmean2, Jmean3
     logical(1) :: refined
     logical(1) :: xyNeighbourPresent, xzNeighbourPresent, yzNeighbourPresent
     integer(1) :: level ! maximum 127 levels of cell refinement
     type(zoneType), pointer :: parent
     type(zoneType), pointer :: xyNeighbour, xzNeighbour, yzNeighbour
     type(zoneType), dimension(:,:,:), pointer :: cell
     type(segmentType), pointer :: segment, lastSegment
     type(patternType), pointer :: pattern
  end type zoneType

  type(zoneType), target :: baseGrid
  type(zoneType), pointer :: localCell
  integer, dimension(33) :: localCallSequence ! maximum 10 levels of refinement

  type :: pointType
     real(kind=RealKind) :: x, y, z
  end type pointType

  type :: pixelType
     logical(1) :: refined
     integer :: level
     real(kind=RealKind) :: phi, theta
     type(pixelType), dimension(:), pointer :: pixel
     type(pixelType), pointer :: parent
  end type pixelType

  type :: readLevelType
     integer :: ncell
     real*4, dimension(:,:), pointer :: pos, vel, abun
     real*4, dimension(:), pointer :: lT, lnH, lx
  end type readLevelType

  type :: starType
     integer :: level, weight
     integer, dimension(:), pointer :: position
     type(zoneType), pointer :: hostCell
     real(kind=RealKind) :: age, xpos, ypos, zpos, rvir, mvir, vcirc, ndm
!     real(kind=RealKind) :: luminosity, mass
     integer*8 :: location ! unique integer identifying position
  end type starType

  type :: segmentType
     real(kind=RealKind) :: len, ndot1
     type(segmentType), pointer :: segment
  end type segmentType

  type :: nodeType
     integer :: n, weight ! n = # of non-degenerate particles per node, weight = # of degenerate particles per node
     integer :: level ! not level of the node, but the level of the cell hosting that node
     real(kind=RealKind) :: xpos, ypos, zpos ! first in units of current node, then in absolute units
     real(kind=RealKind) :: size ! can be (1) physical size of the node box, or (2) distance from the center of the node to its furtherst particle, in absolute units
     type(nodeType), dimension(:,:,:), pointer :: node
     type(nodeType), pointer :: parent
     integer, dimension(:), pointer :: position
     type(zoneType), pointer :: hostCell
  end type nodeType

  type(zoneType), pointer :: densestCell

  integer, parameter :: xyPlane = 1, yzPlane = 2, xzPlane = 3
  integer, parameter :: proceed = 1, split = 2, boundary = 3

  real(kind=RealKind), parameter :: crossSection = 5. ! absorption cross-section

  real(kind=RealKind), parameter :: temstart = 1.
  real(kind=RealKind), parameter :: temend = 1.e8
  integer, parameter :: nratec = 5000 ! # of temperature bins
  integer, parameter :: nfbins = 400 ! # of frequency bins
  integer, parameter :: itmax = 5000 ! max # of chemistry substeps
  real(kind=RealKind), parameter :: frequencyBinWidth = 0.02 ! in log10(eV)
  real(kind=RealKind), parameter :: comp_xraya = 0. ! X-ray compton rates
  real(kind=RealKind), parameter :: comp_temp = 0.

  real(kind=RealKind), dimension(nratec) :: k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, &
       k9a, k10a, k11a, k12a, k13a, k14a, k15a, k16a, &
       k17a, k18a, k19a, k20a, k21a, k22a, k50a, &
       k51a, k52a, k53a, k54a, k55a, k56a, &
       ceHIa, ceHeIa, ceHeIIa, &
       ciHIa, ciHeIa, ciHeISa, ciHeIIa, &
       reHIIa, reHeII1a, reHeII2a, reHeIIIa, brema, lineHIa, &
       hyd01ka, h2k01a, vibha, rotha, rotla, gpldl, gphdl, hdltea, hdlowa

  real(kind=RealKind), dimension(nratec,7) :: k13dda

  real(kind=RealKind) :: logtem0, logtem9, dlogtem, physicalBoxSize

  integer, parameter :: nrmax = 30
  real(kind=RealKind), dimension(nrmax) :: rmax

  real(kind=RealKind), parameter :: psi = 0.76
  real(kind=RealKind), parameter :: starParticleMass = 1.13e6 * msun
  integer, parameter :: temp = 1.e5
! for 1e5 K specific energy (per Hz) peaks at 24.3 eV
!           specific energy (per Angstrom) peaks at 42.8 eV
!           specific photon number (per Hz) peaks at 13.7 eV
  integer, parameter :: nMetallicity = 5, nSpectra = 37, nWavelengths = 1221
!  real(kind=RealKind), dimension(nStars,nSpectra,nWavelengths) :: specificLuminosity
!  real(kind=RealKind), dimension(nSpectra,nWavelengths) :: specificLuminosity
  real(kind=RealKind), dimension(nMetallicity,nSpectra,nWavelengths) :: specificLuminosity
  real(kind=RealKind), dimension(nMetallicity) :: metallicity
  real(kind=RealKind), dimension(nSpectra) :: spectrumTime
  real(kind=RealKind), dimension(nWavelengths) :: wavelength

  integer :: highestPixelLevel, highestNodeLevel
  type(zoneType), pointer :: neighbourCell, splitContinuationCell
  real(kind=RealKind) :: xneighbour, yneighbour, zneighbour
  integer, dimension(33) :: neighbourCallSequence, splitContinuationCallSequence ! maximum 10 levels of cell refinement
  real(kind=RealKind) :: xbase, ybase, zbase, xsrc, ysrc, zsrc
  type(pointType) :: splitContinuation

  type(pointType) :: localPoint

! >3) combined full spectrum at 100 kpc from all sources.
! I think we want 3), including the intensity of the UV light (to estimate
! escape-fractions), We can for one frame calculate average escape fractions
! at, say, d~10, 50 and 100 kpc, just to demonstrate that when the photons
! get out, they get out more or less all the way ...

  integer, parameter :: nenergy = 300
  real(kind=RealKind), parameter :: lowerEnergy = hydrogenIonization
  real(kind=RealKind), parameter :: upperEnergy = 10.*hydrogenIonization
  real(kind=RealKind), dimension(nenergy) :: ndotSpectrum, cosmicSpectrum, outputFreq, &
       outputSigma24, outputSigma25, outputSigma26, outputSigmaDust

!  real(kind=RealKind), dimension(:,:), pointer :: ndotOriginal

  ! HDF variables
  integer sd_id, sds_id, status, n_datasets, n_file_attrs, &
       dfacc_write, dfacc_read, dfacc_create, rank, &
       max_nc_name, max_var_dims, nset
  character*80 fileread
  parameter (dfacc_read = 1, dfacc_write = 2, dfacc_create = 4)
  parameter (max_nc_name = 256, max_var_dims = 5)
  integer start(max_var_dims), edges(max_var_dims), stride(max_var_dims)
  integer n_attrs, data_type
  integer dim_sizes(max_var_dims)
  character sds_name *(max_nc_name)
  integer dfnt_char8, dfnt_uchar8, dfnt_int8, dfnt_uint8, &
       dfnt_int16, dfnt_uint16, dfnt_int32, dfnt_uint32, &
       dfnt_float32, dfnt_float64
  parameter (dfnt_char8 = 4) ! 8-bit character type
  parameter (dfnt_uchar8 = 3) ! 8-bit unsigned character type
  parameter (dfnt_int8 = 20) ! 8-bit integer type
  parameter (dfnt_uint8 = 21) ! 8-bit unsigned integer type
  parameter (dfnt_int16 = 22) ! 16-bit integer type
  parameter (dfnt_uint16 = 23) ! 16-bit unsigned integer type
  parameter (dfnt_int32 = 24) ! 32-bit integer type
  parameter (dfnt_uint32 = 25) ! 32-bit unsigned integer type
  parameter (dfnt_float32 = 5) ! 32-bit floating-point type
  parameter (dfnt_float64 = 6) ! 64-bit floating-point type

  integer*4, dimension(:), pointer :: cellArrayLevel
  real*4, dimension(:), pointer :: cellArrayHI, cellArrayHeI, cellArrayHeII, &
       cellArrayTemp, cellArrayDensity, cellArrayAbun2, &
       cellArrayVelx, cellArrayVely, cellArrayVelz, cellArrayAbun

end module definitions
