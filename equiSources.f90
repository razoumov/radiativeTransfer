!--- Alexei Razoumov, 2005-
!--- Main driver for the Fully Threaded Transport Engine: diffuse and point-source solvers.
!--- Diffuse solver: please cite Razoumov A. and Cardall C.Y., 2005, MNRAS, 362, 1413.
!--- Point-source solver: please cite Razoumov A. and Sommer-Larsen J., 2006, ApJ, 651, L89.

module localDefinitions

  use definitions
  integer, parameter :: nradius = 7, maxPixelLevel = 6
  real(kind=RealKind), dimension(nradius) :: outputRadius = (/ 0.1, 0.3, 1., 3., 10., 30., 100. /) ! [kpc]
  real(kind=RealKind), dimension(nradius) :: ndotRemaining, ndotBoundary, fraction
  real(kind=RealKind) :: ndotDust, maxDepth
  real(kind=RealKind), dimension(:), pointer :: lymanFraction

end module localDefinitions

program pointTransfer

  use definitions
  use localDefinitions
  use stellarPopulationModule, only: stellarPopulation
  use utilities, only: sortStellarParticles
  use dust
  use rotateIndicesModule
  use transportRoutinesModule

  implicit none
  integer :: nlevels, icell, level, i, j, k, nx, ny, nz, itmp, i0, j0, k0, &
       nside, iAngularLevel, imax, ir, itime, iStar, jStar, iSpectrum, iWavelength, &
       iMetal, restart, &
       nStars, nStarsSpecificAge, iradius, ienergy, ilambda, iBadSource, &
       iRedshift, iDensity, iSpecificEnergy, iref, mode, nvariables, counter, &
       reionizationModel, nxtransfer, nytransfer, nztransfer, imean
  integer*8 :: longIntegerArgument, iray
  real(kind=RealKind) :: tmp, tmp1, tmp2, tmp3, x0, y0, z0, xnew, ynew, znew, &
       xa, xb, ya, yb, za, zb, phi, theta, weight, totalIntegral, &
       ndot1, xneu, nh, nhe, coefSpectrum, coefMetal, &
       depth1, depth2, depth3, depthDust, timeReadTable, &
       upperAgeLimit, phiLarge, thetaLarge, phi1, theta1, &
       xmean, ymean, zmean, localTemperature, maxOxygenAbundance, &
       totalOxygen, rho, uvbCoefficient, oldNeutralFraction, actualRate, rateCoef
  real(kind=RealKind) :: Jmean1, Jmean2, Jmean3, cellSizeAbsoluteUnits, Iin1, Iin2, Iin3, &
       dpath, tau1, tau2, tau3, tmpemi1, tmpemi2, tmpemi3, nemi1, nemi2, nemi3, &
       tmpabs1, tmpabs2, tmpabs3
  type(readLevelType), dimension(:), pointer :: readLevel
  real*4, dimension(maxNumberReadVariables) :: readArray
  type(zoneType), pointer :: currentCell
  type(nodeType), pointer :: currentNode, newNode
  type(pixelType), target :: sphere
  type(pixelType), pointer :: currentPixel, tmpPixel, leafPixel
  type(starType), dimension(:), pointer :: star
  type(starType), pointer :: currentStar, nextStar
  type(pointType) :: startingPoint
  type(segmentType), pointer :: newSegment
  type(nodeType), target :: box
  integer :: sfstart, sffinfo, sfselect, sfginfo, sfrdata, sfendacc, sfend, sfcreate, sfwdata
  integer, dimension(3) :: baseGridSize
  logical :: readNewSpectrum
  character(60) :: synthesisFile, synthfile(5), sphDir, synthesisDir, grid, sources
  character(90) :: scratchString, restartCellArrayName, filename
  real(kind=RealKind), parameter :: alphaQuasar = 1.8
  real(kind=RealKind), parameter :: alphaStellar = 5.
  real(kind=RealKind), parameter :: contributionQuasar = 1.
  real(kind=RealKind), parameter :: contributionStellar = 1.
  integer, parameter :: stellarTransferThinUVB = 1, plotPDFs = 2, initialConfiguration = 3, &
       printNumberOfCells = 4, noStarsThinUVB = 6, clumpingFactor = 7, &
       bothStellarUVBTransfer = 8, UVBTransferOnly = 9
  logical :: runStellarTransfer, runUVBTransfer
  integer :: nxmap, nymap, kstart, kend, jcell, kcell
  integer(1) :: izone
  integer, dimension(2,2,2) :: is, js, ks
  real*4, dimension(:,:), pointer :: mapArray
  real(kind=RealKind), dimension(2) :: startPoint, endPoint
  real(kind=RealKind), dimension(3) :: centerPoint
  real(kind=RealKind) :: zoomFactor, zstart, zend
  character(60) :: dirname, mapname
  real(kind=RealKind), dimension(:,:,:), pointer :: unigrid, tmpgrid
  real(kind=RealKind), dimension(17) :: z, rate
  type(angleType) :: angles
  real(kind=RealKind), dimension(3) :: alpha
  type(patternType), dimension(:), pointer :: pattern

  synthfile(1) = 'model41-salpeter-burst34/spectrum.out'
  synthfile(2) = 'model42-salpeter-burst34/spectrum.out'
  synthfile(3) = 'model43-salpeter-burst34/spectrum.out'
  synthfile(4) = 'model44-salpeter-burst34/spectrum.out'
  synthfile(5) = 'model45-salpeter-burst34/spectrum.out'

  dustApproximation = 0
  selfShieldingThreshold = 1.*kpc
  currentRedshift = 3.
  massStellarParticle = 1 ! normal(8x) = 1, hiRes(64x) = 2, superHiRes(512x) = 3, massive = 10, hiResHeavy(64x_heavy) = 4, crazyHiRes(512x8x) = 5
  grid = ''
  sources = ''
  restartCellArrayName = ''
  mode = 1
  upperAgeLimit = 10.*Myr      !   34.*Myr   10.*Myr   3.4*Myr
  sphDir = ''
  synthesisDir = ''
  open(unit=10, file='inputParameters', action="read", status="old")
  status = 0
  restart = 0
  uvbCoefficient = 1.
  reionizationModel = 0
  do while (status.eq.0)
     read(unit=10, fmt="(a90)", IOSTAT=status) scratchString
     if (scratchString(1:20).eq.'dustApproximation = ') read (scratchString(21:90),*) dustApproximation
     if (scratchString(1:25).eq.'selfShieldingThreshold = ') then
        read (scratchString(26:90),*) selfShieldingThreshold
        selfShieldingThreshold = selfShieldingThreshold*kpc
     endif
     if (scratchString(1:18).eq.'currentRedshift = ') read (scratchString(19:90),*) currentRedshift
     if (scratchString(1:17).eq.'uvbCoefficient = ') read (scratchString(18:90),*) uvbCoefficient
     if (scratchString(1:22).eq.'massStellarParticle = ') read (scratchString(23:90),*) massStellarParticle
     if (scratchString(1:7).eq.'grid = ') read (scratchString(8:90),*) grid
     if (scratchString(1:10).eq.'sources = ') read (scratchString(11:90),*) sources
     if (scratchString(1:7).eq.'mode = ') read (scratchString(8:90),*) mode
     if (scratchString(1:15).eq.'upperAgeLimit = ') then
        read (scratchString(16:90),*) upperAgeLimit
        upperAgeLimit = upperAgeLimit * Myr
     endif
     if (scratchString(1:9).eq.'sphDir = ') read (scratchString(10:90),*) sphDir
     if (scratchString(1:15).eq.'synthesisDir = ') read (scratchString(16:90),*) synthesisDir
     if (scratchString(1:10).eq.'restart = ') read (scratchString(11:90),*) restart
     if (scratchString(1:23).eq.'restartCellArrayName = ') read (scratchString(24:90),*) restartCellArrayName
     if (scratchString(1:20).eq.'reionizationModel = ') read (scratchString(21:90),*) reionizationModel
  enddo
  close(10)
  print*, 'dustApproximation =', dustApproximation
  print*, 'selfShieldingThreshold =', selfShieldingThreshold/kpc
  print*, 'currentRedshift =', currentRedshift
  print*, 'massStellarParticle =', massStellarParticle
  print*, 'grid = ', trim(grid)
  print*, 'sources = ', trim(sources)
  print*, 'mode = ', mode
  print*, 'upperAgeLimit = ', upperAgeLimit/Myr
  print*, 'sphDir = ', trim(sphDir)
  print*, 'synthesisDir = ', trim(synthesisDir)
  print*, 'restart = ', restart
  print*, 'restartCellArrayName = ', trim(restartCellArrayName)
  print*, 'uvbCoefficient = ', uvbCoefficient
  print*, 'reionizationModel =', reionizationModel

  readKinematics = .false.
  itmp = len(trim(grid))
  do i = 1, itmp
     if (grid(i:i).eq.'v'.and.i.lt.itmp-1) then
        if (grid(i:i+2).eq.'vel') readKinematics = .true.
     endif
  enddo

  readMetals = .false.
  itmp = len(trim(grid))
  do i = 1, itmp
     if (grid(i:i).eq.'m'.and.i.lt.itmp-1) then
        if (grid(i:i+2).eq.'met') readMetals = .true.
     endif
  enddo

  if (mode.eq.stellarTransferThinUVB .or. mode.eq.bothStellarUVBTransfer) then
     runStellarTransfer = .true.
  else
     runStellarTransfer = .false.
  endif

  if (mode.eq.UVBTransferOnly .or. mode.eq.bothStellarUVBTransfer) then
     runUVBTransfer = .true.
  else
     runUVBTransfer = .false.
  endif

! initialize collisional/recombination rate data

  logtem0 = log(temstart)
  logtem9 = log(temend)
  dlogtem = (log(temend) - log(temstart))/real(nratec-1)

  call calc_rates(nratec, temstart, temend, &
       ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, &
       ciHeISa, ciHeIIa, reHIIa, reHeII1a, &
       reHeII2a, reHeIIIa, brema, lineHIa, compa, &
       hyd01ka, h2k01a, vibha, rotha, rotla, &
       gpldl, gphdl, hdltea, hdlowa, &
       k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a, &
       k11a, k12a, k13a, k13dda, k14a, k15a, k16a, k17a, &
       k18a, k19a, k20a, k21a, k22a, &
       k50a, k51a, k52a, k53a, k54a, k55a, k56a, &
       recombinationType)

! compute radiative rates for the uniform fiducial background

  call uniformTable(nfbins,frequencyBinWidth,alphaQuasar,alphaStellar)

! set background radiation

! Abel & Haehnelt 99 stellar component

  stellar99 = 1./(1.+(7./(1.+currentRedshift))**4)*exp(-(currentRedshift/4.)**3)

! Abel & Haehnelt 99 quasar component

  quasar99 = 10./(1.+(7./(1.+currentRedshift))**4)*exp(-(currentRedshift/2.5)**3)

! Paschos 02 total

  pascal02 = 0.0188*exp(-(currentRedshift-0.5)**2/(1.+0.0625*(currentRedshift+2.09)**2.075))*(1.+currentRedshift)**3.35

! Razoumov 02 stellar component

  component1 = stellar99
  component2 = pascal02
  step = 0.5*(tanh((currentRedshift-4.2)*1.5)+1.)
  stellar02 = (1.-step)*component1 + step*component2

! Razoumov 02 quasar component

  quasar02 = 10./(1.+(7./(1.+currentRedshift))**4)*exp(-(currentRedshift/2.5)**3)

  gaussian = exp(-((currentRedshift-4.5)/2.)**2)*0.3

  newQuasar02 = gaussian*stellar02 + (1.-gaussian)*quasar02
  newStellar02 = (1.-gaussian)*stellar02 + gaussian*quasar02

  step = 0.5*(tanh((currentRedshift-14.)*0.5)+1.)
  newStellar02 = step*0. + (1.-step)*newStellar02

! Razoumov 02 total

  total02 = newStellar02 + newQuasar02

  uniformQuasar = newQuasar02 * 1.d-21 * contributionQuasar * uvbCoefficient
  uniformStellar = newStellar02 * 1.d-21 * contributionStellar * uvbCoefficient

  if (runUVBTransfer) then

     uvbStellar1 = newStellar02 * 1.d-21 * contributionStellar * uvbCoefficient
     uvbStellar2 = uvbStellar1 * (nu2/nu1)**(-alphaStellar)
     uvbStellar3 = uvbStellar2 * (nu3/nu2)**(-alphaStellar)

     uvbQuasar1 = newQuasar02 * 1.d-21 * contributionQuasar * uvbCoefficient
     uvbQuasar2 = uvbQuasar1 * (nu2/nu1)**(-alphaQuasar)
     uvbQuasar3 = uvbQuasar2 * (nu3/nu2)**(-alphaQuasar)

     call powerSpectrumIndex(uvbStellar1,alphaStellar,uvbQuasar1,alphaQuasar,uvb1,alpha(1),nu1,nu2,.true.)
     call powerSpectrumIndex(uvbStellar2,alphaStellar,uvbQuasar2,alphaQuasar,uvb2,alpha(2),nu2,nu3,.true.)
     call powerSpectrumIndex(uvbStellar3,alphaStellar,uvbQuasar3,alphaQuasar,uvb3,alpha(3),nu3,nu3,.false.)

     write(*,'("band1: ",3f15.10)') uvbStellar1/1.d-21, uvbQuasar1/1.d-21, uvb1/1.d-21
     write(*,'("band2: ",3f15.10)') uvbStellar2/1.d-21, uvbQuasar2/1.d-21, uvb2/1.d-21
     write(*,'("band3: ",3f15.10)') uvbStellar3/1.d-21, uvbQuasar3/1.d-21, uvb3/1.d-21
     print*, 'alpha =', alpha

     call uvbBetaTable(nfbins,frequencyBinWidth,alpha)

  endif

! hydrogen ionization rates for models z_re=6 and z_re=10

  if (reionizationModel.ne.0) then
     print*, 'assuming reionization at z=', reionizationModel
     select case (reionizationModel)
     case (6)
        z = (/ 0., 0.316, 0.697, 1.187, 1.513, 2.343, 2.547, 2.765, &
             3.024, 3.296, 3.772, 4.316, 4.657, 4.997, 5.302, 5.609, 100. /)
        rate = (/ 0.0045, 0.0100, 0.0248, 0.0585, 0.0968, 0.1594, 0.1621, 0.1564, &
             0.1403, 0.1159, 0.0683, 0.0248, 0.0112, 0.0058, 0.0017, 0.0004, 0.0000 /) * 1.d-11
     case (10)
        z = (/ 0., 0.316, 0.697, 1.187, 1.513, 2.343, 2.547, 2.972, &
             3.432, 3.976, 5.065, 6.221, 6.902, 7.650, 8.331, 9.419, 100. /)
        rate = (/ 0.0045, 0.0100, 0.0248, 0.0585, 0.0968, 0.1594, 0.1621, 0.1570, &
             0.1444, 0.1240, 0.0710, 0.0262, 0.0128, 0.0058, 0.0014, 0.0003, 0.0000 /) * 1.d-11
     case default
        print*, 'select proper reionization model'
        stop
     endselect
     i = 1
     do while (currentRedshift.gt.z(i))
        i = i + 1
     enddo
     actualRate = (currentRedshift-z(i-1)) / (z(i)-z(i-1)) * (rate(i)-rate(i-1)) + rate(i-1)
     rateCoef = actualRate / (4.*pi * (uniformQuasar*quasar%ksi24 + uniformStellar*stellar%ksi24))
     uniformQuasar = uniformQuasar * rateCoef
     uniformStellar = uniformStellar * rateCoef
     if (runUVBTransfer) then
        uvb1 = uvb1 * rateCoef
        uvb2 = uvb2 * rateCoef
        uvb3 = uvb3 * rateCoef
     endif
  endif

  call dustInitialize ! initialize dust data

  rmax(1)  = 1.984
  rmax(2)  = 3.264
  rmax(3)  = 5.739
  rmax(4)  = 10.65
  rmax(5)  = 20.46
  rmax(6)  = 40.05
  rmax(7)  = 79.25
  rmax(8)  = 157.6
  rmax(9)  = 314.4
  rmax(10) = 627.9

  do ir = 1, nrmax
     rmax(ir) = sqrt(3.)*(sqrt(0.5*4.**(ir-1)-1./12.)+0.5)
     if (mode.ne.printNumberOfCells) print*, ir, rmax(ir)
  enddo

  rmax = rmax/2.

  print*, 'maxPixelLevel =', maxPixelLevel

  sphere%level = 0
  sphere%refined = .false.

  sd_id = sfstart(trim(sphDir)//trim(grid)//'.h4', dfacc_read)

  status = sffinfo(sd_id, n_datasets, n_file_attrs)

  start = 0
  edges(1) = 1
  stride = 1

  sds_id = sfselect(sd_id, 0)
  status = sfrdata(sds_id, start, stride, edges, nlevels)
  status = sfendacc(sds_id)

  nvariables = 4
  if (readKinematics) nvariables = nvariables + 1
  if (readMetals) nvariables = nvariables + 1

  if (nlevels.ne.(n_datasets-1)/nvariables) then
     print*, trim(sphDir)//trim(grid)//'.h4'
     print*, 'error in nlevels', nlevels, n_datasets
     stop
  endif

  allocate(readLevel(nlevels))

  if (readMetals) then
     maxOxygenAbundance = 0.
     totalOxygen = 0.
  endif

  do level = 1, nlevels

     sds_id = sfselect(sd_id, nvariables*(level-1)+1)
     status = sfginfo(sds_id,sds_name,rank,dim_sizes,data_type,n_attrs)
!     if (mode.ne.printNumberOfCells) print*, trim(sds_name), rank, dim_sizes(1), dim_sizes(2)
     if (status.ne.0 .or. rank.ne.2) then
        write(*,*) 'error reading ...', level, status, rank, nlevels
        stop
     endif

     readLevel(level)%ncell = dim_sizes(1)
     allocate(readLevel(level)%pos(readLevel(level)%ncell,3))
     allocate(readLevel(level)%lT(readLevel(level)%ncell))
     allocate(readLevel(level)%lnH(readLevel(level)%ncell))
     allocate(readLevel(level)%lx(readLevel(level)%ncell))
     if (readKinematics) allocate(readLevel(level)%vel(readLevel(level)%ncell,3))
     if (readMetals) allocate(readLevel(level)%abun(readLevel(level)%ncell,4))

!     print*, 'level =', level, '    ncells =', readLevel(level)%ncell

     edges(1:2) = (/ readLevel(level)%ncell, 3 /)
     status = sfrdata(sds_id, start, stride, edges, readLevel(level)%pos)
     status = sfendacc(sds_id)

     sds_id = sfselect(sd_id, nvariables*(level-1)+2)
     edges(1) = readLevel(level)%ncell
     status = sfrdata(sds_id, start, stride, edges, readLevel(level)%lT)
     status = sfendacc(sds_id)

     sds_id = sfselect(sd_id, nvariables*(level-1)+3)
     edges(1) = readLevel(level)%ncell
     status = sfrdata(sds_id, start, stride, edges, readLevel(level)%lnH)
     status = sfendacc(sds_id)

     if (mode.eq.printNumberOfCells) then
        tmp = 0.
        do icell = 1, readLevel(level)%ncell
           tmp = max(tmp,readLevel(level)%lnH(icell))
        enddo
        write(*,*) 'level =', level, readLevel(level)%ncell, tmp
     endif

     sds_id = sfselect(sd_id, nvariables*(level-1)+4)
     edges(1) = readLevel(level)%ncell
     status = sfrdata(sds_id, start, stride, edges, readLevel(level)%lx)
     status = sfendacc(sds_id)

     counter = 5

     if (readMetals) then
        sds_id = sfselect(sd_id, nvariables*(level-1)+counter)
        edges(1:2) = (/ readLevel(level)%ncell, 4 /)
        status = sfrdata(sds_id, start, stride, edges, readLevel(level)%abun)
        status = sfendacc(sds_id)
        counter = counter + 1
     endif

     if (readKinematics) then
        sds_id = sfselect(sd_id, nvariables*(level-1)+counter)
        edges(1:2) = (/ readLevel(level)%ncell, 3 /)
        status = sfrdata(sds_id, start, stride, edges, readLevel(level)%vel)
        status = sfendacc(sds_id)
        counter = counter + 1
     endif

     if (readMetals) then
        do icell = 1, readLevel(level)%ncell
           rho = 10.**readLevel(level)%lnH(icell)
           if (readLevel(level)%abun(icell,2).gt.maxOxygenAbundance) &
                maxOxygenAbundance = readLevel(level)%abun(icell,2)
           totalOxygen = totalOxygen + readLevel(level)%abun(icell,2)*rho/(2.**(3*(level-1)))
        enddo
     endif

  enddo

  if (readMetals) print*, 'maximum oxygen abundance =', maxOxygenAbundance, totalOxygen

  status = sfend(sd_id)

  if (mode.eq.printNumberOfCells) stop

  i = 0
  itmp = 0
  do while (itmp.lt.readLevel(1)%ncell)
     i = i + 1
     itmp = i**3
  enddo
  if (itmp.ne.readLevel(1)%ncell) then
     write(*,*) 'base grid needs to be of size n^3', itmp, readLevel(1)%ncell
     stop
  endif

  nx = i
  ny = i
  nz = i
  write(*,*) 'grid:', nx,'^3'

  print*, 'nlevels =', nlevels
  do i = 1, nlevels
     print*, 'level ', i, ' has', readLevel(i)%ncell, ' cells'
  enddo

  level = 1
  xa = 1.d10
  xb = - 1.d10
  ya = 1.d10
  yb = - 1.d10
  za = 1.d10
  zb = - 1.d10
  do icell = 1, readLevel(level)%ncell
     xa = min(xa,dble(readLevel(level)%pos(icell,1)))
     xb = max(xb,dble(readLevel(level)%pos(icell,1)))
     ya = min(ya,dble(readLevel(level)%pos(icell,2)))
     yb = max(yb,dble(readLevel(level)%pos(icell,2)))
     za = min(za,dble(readLevel(level)%pos(icell,3)))
     zb = max(zb,dble(readLevel(level)%pos(icell,3)))
  enddo

  tmp1 = 0.5*(xa+xb)
  tmp2 = 0.5*(xb-xa)*float(nx)/float(nx-1)
  xa = tmp1 - tmp2
  xb = tmp1 + tmp2

  tmp1 = 0.5*(ya+yb)
  tmp2 = 0.5*(yb-ya)*float(ny)/float(ny-1)
  ya = tmp1 - tmp2
  yb = tmp1 + tmp2

  tmp1 = 0.5*(za+zb)
  tmp2 = 0.5*(zb-za)*float(nz)/float(nz-1)
  za = tmp1 - tmp2
  zb = tmp1 + tmp2

  write(*,*) 'edges of computational grid in kpc:'
  write(*,*) xa, xb
  write(*,*) ya, yb
  write(*,*) za, zb
  physicalBoxSize = abs(xa-xb)*kpc

  do level = 1, nlevels
     do icell = 1, readLevel(level)%ncell
        readLevel(level)%pos(icell,1) = (readLevel(level)%pos(icell,1)-xa)/(xb-xa)
        readLevel(level)%pos(icell,2) = (readLevel(level)%pos(icell,2)-ya)/(yb-ya)
        readLevel(level)%pos(icell,3) = (readLevel(level)%pos(icell,3)-za)/(zb-za)
     enddo
  enddo

  print*, 'setting up a fully nested 3D grid'

  allocate(baseGrid%cell(nx,ny,nz))

  do i = 1, nx
     do j = 1, ny
        do k = 1, nz

           baseGrid%cell(i,j,k)%refined = .false.
           baseGrid%cell(i,j,k)%tgas = 0.d0
           baseGrid%cell(i,j,k)%rho = 0.d0
           baseGrid%cell(i,j,k)%HI = 0.d0
           baseGrid%cell(i,j,k)%HeI = 0.d0
           baseGrid%cell(i,j,k)%HeII = 0.d0
           baseGrid%cell(i,j,k)%rhoCoef = 1.
           baseGrid%cell(i,j,k)%velx = 0.d0
           baseGrid%cell(i,j,k)%vely = 0.d0
           baseGrid%cell(i,j,k)%velz = 0.d0
           baseGrid%cell(i,j,k)%krate24 = 0.d0
           baseGrid%cell(i,j,k)%krate25 = 0.d0
           baseGrid%cell(i,j,k)%krate26 = 0.d0
           baseGrid%cell(i,j,k)%crate24 = 0.d0
           baseGrid%cell(i,j,k)%crate25 = 0.d0
           baseGrid%cell(i,j,k)%crate26 = 0.d0
           baseGrid%cell(i,j,k)%level = 0
           baseGrid%cell(i,j,k)%parent => baseGrid
           nullify(baseGrid%cell(i,j,k)%cell)
           nullify(baseGrid%cell(i,j,k)%segment)
           nullify(baseGrid%cell(i,j,k)%lastSegment)

        enddo
     enddo
  enddo

! smoothing level 1 abundancies
  level = 1
  allocate(unigrid(nx,ny,nz))
  allocate(tmpgrid(nx,ny,nz))
  do icell = 1, readLevel(level)%ncell
     i0 = int(readLevel(level)%pos(icell,1)*nx) + 1
     j0 = int(readLevel(level)%pos(icell,2)*ny) + 1
     k0 = int(readLevel(level)%pos(icell,3)*nz) + 1
     unigrid(i0,j0,k0) = readLevel(level)%abun(icell,2)
  enddo
  do itmp = 1, 2
     tmpgrid = 0.
     do i = 1, nx
        do j = 1, ny
           do k = 1, nz
              tmpgrid(i,j,k) = tmpgrid(i,j,k) + 0.5*unigrid(i,j,k)
              if (i.gt.1)  tmpgrid(i-1,j,k) = tmpgrid(i-1,j,k) + 0.25*unigrid(i,j,k)
              if (i.lt.nx) tmpgrid(i+1,j,k) = tmpgrid(i+1,j,k) + 0.25*unigrid(i,j,k)
           enddo
        enddo
     enddo
     unigrid = tmpgrid
     tmpgrid = 0.
     do i = 1, nx
        do j = 1, ny
           do k = 1, nz
              tmpgrid(i,j,k) = tmpgrid(i,j,k) + 0.5*unigrid(i,j,k)
              if (j.gt.1)  tmpgrid(i,j-1,k) = tmpgrid(i,j-1,k) + 0.25*unigrid(i,j,k)
              if (j.lt.ny) tmpgrid(i,j+1,k) = tmpgrid(i,j+1,k) + 0.25*unigrid(i,j,k)
           enddo
        enddo
     enddo
     unigrid = tmpgrid
     tmpgrid = 0.
     do i = 1, nx
        do j = 1, ny
           do k = 1, nz
              tmpgrid(i,j,k) = tmpgrid(i,j,k) + 0.5*unigrid(i,j,k)
              if (k.gt.1)  tmpgrid(i,j,k-1) = tmpgrid(i,j,k-1) + 0.25*unigrid(i,j,k)
              if (k.lt.nz) tmpgrid(i,j,k+1) = tmpgrid(i,j,k+1) + 0.25*unigrid(i,j,k)
           enddo
        enddo
     enddo
     unigrid = tmpgrid
  enddo
  do icell = 1, readLevel(level)%ncell
     i0 = int(readLevel(level)%pos(icell,1)*nx) + 1
     j0 = int(readLevel(level)%pos(icell,2)*ny) + 1
     k0 = int(readLevel(level)%pos(icell,3)*nz) + 1
     readLevel(level)%abun(icell,2) = unigrid(i0,j0,k0)
  enddo
  deallocate(unigrid,tmpgrid)

  do level = 1, nlevels
     do icell = 1, readLevel(level)%ncell

        x0 = readLevel(level)%pos(icell,1)
        y0 = readLevel(level)%pos(icell,2)
        z0 = readLevel(level)%pos(icell,3)

        i0 = int(x0*nx) + 1
        j0 = int(y0*ny) + 1
        k0 = int(z0*nz) + 1

        xnew = x0*float(nx) - float(i0-1)
        ynew = y0*float(ny) - float(j0-1)
        znew = z0*float(nz) - float(k0-1)

        readArray(1)  = readLevel(level)%lT(icell)
        readArray(2)  = readLevel(level)%lnH(icell)
        readArray(3)  = readLevel(level)%lx(icell)
        nvariables = 3

        if (readKinematics) then
           readArray(nvariables+1) = readLevel(level)%vel(icell,1)
           readArray(nvariables+2) = readLevel(level)%vel(icell,2)
           readArray(nvariables+3) = readLevel(level)%vel(icell,3)
           nvariables = nvariables + 3
        endif

        if (readMetals) then
           readArray(nvariables+1) = readLevel(level)%abun(icell,1)
           readArray(nvariables+2) = readLevel(level)%abun(icell,2)
           readArray(nvariables+3) = readLevel(level)%abun(icell,3)
           readArray(nvariables+4) = readLevel(level)%abun(icell,4)
           nvariables = nvariables + 4
        endif

        call placeCellProjectWithVelocity(baseGrid%cell(i0,j0,k0),level,xnew,ynew,znew,readArray)

     enddo
  enddo

  do level = 1, nlevels
     deallocate(readLevel(level)%pos)
     deallocate(readLevel(level)%lT)
     deallocate(readLevel(level)%lnH)
     deallocate(readLevel(level)%lx)
     if (readKinematics) deallocate(readLevel(level)%vel)
     if (readMetals) deallocate(readLevel(level)%abun)
  enddo
  deallocate(readLevel)

!   densestCell => baseGrid%cell(1,1,1)
!   do i = 1, nx
!      do j = 1, ny
!         do k = 1, nz
!            call printCell(baseGrid%cell(i,j,k),0,(/i,j,k/))
!         enddo
!      enddo
!   enddo

!   densestCell => baseGrid%cell(1,1,1)
!   write(*,*) densestCell%nh, densestCell%level
!   do i = 1, nx
!      do j = 1, ny
!         do k = 1, nz
!            call countCells(baseGrid%cell(i,j,k))
!         enddo
!      enddo
!   enddo
!   write(*,*) 'found', icosmic, ' cells'
!   write(*,*) densestCell%nh, densestCell%level

!   currentCell => baseGrid%cell(1,1,1)
!   if (associated(currentCell%lastSegment)) then
!      allocate(currentCell%lastSegment%segment)
!      currentCell%lastSegment => currentCell%lastSegment%segment
!   else
!      allocate(currentCell%segment)
!      currentCell%lastSegment => currentCell%segment
!   endif
!   newSegment => currentCell%lastSegment

  if (mode.eq.clumpingFactor) then
     volSum = 0.
     volSumDensity = 0.
     volSumDensity2 = 0.
     do i = 1, nx
        do j = 1, ny
           do k = 1, nz
              call computeClumping(baseGrid%cell(i,j,k))
           enddo
        enddo
     enddo
     volSumDensity = volSumDensity / volSum
     volSumDensity2 = volSumDensity2 / volSum
     print*, 'clumping =', volSumDensity2 / volSumDensity**2, volSum
     stop
  endif

  if (mode.eq.initialConfiguration) then
     izone = 3
     nxmap = 1024
     nymap = nxmap
     centerPoint = (/ 0.5, 0.5, 0.5 /)
     zoomFactor = 2.
     mapname = 'abun-S33-06-64-10.4.h4'
     startPoint = (/ max(centerPoint(1) - 0.5/zoomFactor,0.), &
          max(centerPoint(2) - 0.5/zoomFactor,0.) /)
     zstart = max(centerPoint(3) - 0.5/zoomFactor,0.)
     zend = min(centerPoint(3) + 0.5/zoomFactor,1.)
     endPoint(1:2) = min(startPoint(1:2)+1./zoomFactor,1.)
     do i = 1, 2
        do j = 1, 2
           do k = 1, 2
              call rotateIndices(i,j,k,2,2,2,izone,is(i,j,k),js(i,j,k),ks(i,j,k))
           enddo
        enddo
     enddo
     allocate(mapArray(nxmap,nymap))
     endPoint(1:2) = min(startPoint(1:2)+1./zoomFactor,1.)
     kstart = max(int(zstart*nz)+1,1)
     kend = min(int(zend*nz)+1,nz)
     do i = 1, nxmap
        print*, 'doing', i, ' out of', nxmap
        x0 = (endPoint(1)-startPoint(1))*(float(i)-0.5)/float(nxmap) + startPoint(1)
        i0 = int(x0*nx) + 1
        do j = 1, nymap
           y0 = (endPoint(2)-startPoint(2))*(float(j)-0.5)/float(nymap) + startPoint(2)
           j0 = int(y0*ny) + 1
           xnew = x0*float(nx) - float(i0-1)
           ynew = y0*float(ny) - float(j0-1)
           mapArray(i,j) = 0.
           totalMass = 0.
           do k0 = kstart, kend
              call rotateIndices(i0,j0,k0,nx,ny,nz,izone,icell,jcell,kcell)
              call projectVariableToMap(baseGrid%cell(icell,jcell,kcell),mapArray(i,j), &
                   xnew,ynew,physicalBoxSize/dfloat(nx),is,js,ks)
           enddo
           mapArray(i,j) = mapArray(i,j) / totalMass
        enddo
     enddo
     sd_id = sfstart(trim(mapname),dfacc_create)
     edges(1) = nxmap
     edges(2) = nymap
     sds_id = sfcreate(sd_id,'map',dfnt_float32,2,edges)
     start = 0
     stride = 1
     status = sfwdata(sds_id,start,stride,edges,mapArray)
     status = sfendacc(sds_id)
     status = sfend(sd_id)
     deallocate(mapArray)
     stop
  endif

  print*, 'read positions and ages of all stars'
  open(unit=15,file=trim(sphDir)//trim(sources),action='read',status='old')
  status = 0
  iStar = 0
  do while (status.eq.0)
     read(15,*,iostat=status) j
     iStar = iStar + 1
  enddo
  nStars = iStar - 1
  close(15)
  allocate(star(nStars))
  open(unit=15,file=trim(sphDir)//trim(sources),action='read',status='old')
  maxOxygenAbundance = 0.
  do iStar = 1, nStars
     currentStar => star(iStar)
     read(15,*,iostat=status) currentStar%level, tmp1, tmp2, tmp3, currentStar%age
     currentStar%age = currentStar%age * Myr
     xbase = (tmp1-xa)/(xb-xa)
     ybase = (tmp2-ya)/(yb-ya)
     zbase = (tmp3-za)/(zb-za)
     call localizeSplitContinuationCell(xbase,ybase,zbase,0,nx,ny,nz)
!     print*, iStar, currentStar%level-1, splitContinuationCell%level
     currentStar%level = splitContinuationCell%level  ! overwrite the levels
     allocate(currentStar%position(3*currentStar%level+3))
     currentStar%position = splitContinuationCallSequence(1:3*currentStar%level+3)
     startingPoint = splitContinuation
!     print*, '--------'
!     print*, 'source #', iStar
!     print*, 'cube =', xbase, ybase, zbase
!     print*, 'level =', splitContinuationCell%level, '    rho =', splitContinuationCell%rho
!     print*, 'stellar metallicity =', splitContinuationCell%abun2
!     print*, iStar, '      level =', currentStar%level
!     print*, 'position =', currentStar%position
!     print*, 'point =', startingPoint
     if (splitContinuationCell%abun2.gt.maxOxygenAbundance) &
          maxOxygenAbundance = splitContinuationCell%abun2
  enddo
  print*, 'maximum stellar oxygen abundance =', maxOxygenAbundance
  close(15)

  nStarsSpecificAge = 0
  do iStar = 1, nStars
     currentStar => star(iStar)
     if (currentStar%age.le.upperAgeLimit) then
        nStarsSpecificAge = nStarsSpecificAge + 1
        currentStar%weight = 1
!        print*, iStar, currentStar%level, currentStar%age/Myr
     else
        currentStar%weight = 0
     endif
  enddo

  if (mode.eq.plotPDFs) then
     print*, 'computing stellar and gas PDF'
     meanHostDensity = 0.
     pdfStar = 0
     pdfStarOutside = 0
     do iStar = 1, nStars
        currentStar => star(iStar)
        if (currentStar%weight.gt.0) then
           startingPoint%x = 0.5
           startingPoint%y = 0.5
           startingPoint%z = 0.5
           i = currentStar%position(1)
           j = currentStar%position(2)
           k = currentStar%position(3)
           call localizeCellFromStar(baseGrid%cell(i,j,k),currentStar)
           currentCell => currentStar%hostCell
           write(*,2015) iStar, currentCell%level, psi*currentCell%rho/mh
           meanHostDensity = meanHostDensity + currentCell%HI
           tmp = dlog10(currentCell%rho/msun*pc**3)
           if (tmp.gt.apdf .and. tmp.lt.bpdf) then
              itmp = int((tmp-apdf)/(bpdf-apdf)*float(npdf)) + 1
              pdfStar(itmp) = pdfStar(itmp) + 1
           else
              pdfStarOutside = pdfStarOutside + 1
           endif
2015       format(2i5,es18.8)
        endif
     enddo
     print*, 'meanHostDensity =', meanHostDensity/float(nStarsSpecificAge), nStarsSpecificAge
     pdfGas = 0.
     pdfGasOutside = 0.
     do i = 1, nx
        do j = 1, ny
           do k = 1, nz
              call computeGasPDF(baseGrid%cell(i,j,k))
           enddo
        enddo
     enddo
     itmp = pdfStarOutside
     tmp = pdfGasOutside
     do i = 1, npdf
        itmp = itmp + pdfStar(i)
        tmp = tmp + pdfGas(i)
     enddo
     do i = 1, npdf
        tmp1 = max(1.d-38,pdfGas(i)/tmp)
        tmp2 = max(1.d-38,float(pdfStar(i))/float(itmp))
        write(*,2025) (float(i)-0.5)/float(npdf)*(bpdf-apdf)+apdf, dlog10(tmp1), dlog10(tmp2)
     enddo
2025 format(3f15.6)
     stop
  endif

! initialize stellar population synthesis data

  if (nMetallicity.ne.5) then
     print*, 'please update the list of metallicities below'
     stop
  endif
  metallicity(1:5) = (/ 0.0004, 0.004, 0.008, 0.020, 0.050 /)
  metallicity = dlog10(metallicity)

  do iMetal = 1, nMetallicity

     filename = synthfile(iMetal)

     open (unit=11,file=trim(synthesisDir)//trim(filename),action='read',status='old')

     status = 0
     i = 0
     readNewSpectrum = .false.
     iSpectrum = 0
     do while (status.eq.0)
        i = i + 1
        read(unit=11, fmt='(a65)', IOSTAT=status) scratchString
        if (status.eq.0) then
!     print*, i, scratchString(1:65)
           if (scratchString(2:10).eq.'TIME [YR]') then
              read(unit=11, fmt='(a65)', IOSTAT=status) scratchString
              read(unit=11, fmt='(a65)', IOSTAT=status) scratchString
              readNewSpectrum = .true.
              iWavelength = 0
           else
              if (readNewSpectrum .and. scratchString(2:6).ne.'MODEL') then
                 if (iWavelength.eq.0) iSpectrum = iSpectrum + 1
                 iWavelength = iWavelength + 1
                 read (scratchString(2:13),*) spectrumTime(iSpectrum)
                 spectrumTime(iSpectrum) = spectrumTime(iSpectrum) * yr
                 read (scratchString(14:28),*) wavelength(iWavelength)
                 wavelength(iWavelength) = wavelength(iWavelength) * angstrom
                 read (scratchString(29:41),*) specificLuminosity(iMetal,iSpectrum,iWavelength)
              else
                 readNewSpectrum = .false.
              endif
           endif
        endif
     enddo
     close(11)

  enddo

!    specificLuminosity was computed for the SF rate of 11.6 msun/yr distributed among 34 star particles
!    specificLuminosity * 34 -- luminosity of 11.6 msun/yr constant SF rate
!    nStars / 347 * 11.6 Msun/yr -- SF rate in the volume
!    nStars / 347 * 10.**specificLuminosity * 34 -- luminosity of the volume
!    nStars / 347 * 10.**specificLuminosity * 34 / nStarsSpecificAge -- luminosity per star particle

  specificLuminosity = specificLuminosity + log10(float(nStars) / 347. * 34. / float(nStarsSpecificAge)) ! per star particle

  select case (massStellarParticle)
  case (hiRes)
     print*, 'high-res model: stellar sources are 1/8th the mass of normal model'
     specificLuminosity = specificLuminosity - log10(8.)
  case (superHiRes)
     print*, 'high-res model: stellar sources are 1/64th the mass of normal model'
     specificLuminosity = specificLuminosity - log10(64.)
  case (crazyHiRes)
     print*, 'high-res model: stellar sources are 1/512th the mass of normal model'
     specificLuminosity = specificLuminosity - log10(512.)
  case (massive)
     print*, 'massive star particles: stellar sources are 2.7818X the mass of normal model'
     specificLuminosity = specificLuminosity + log10(2.7818)
  case (hiResHeavy)
     print*, 'high-res model: stellar sources are (5.832/8)X the mass of normal model'
     specificLuminosity = specificLuminosity + log10(5.832/8.)
  case (light)
     print*, 'high-res model: stellar sources are 0.6^3/512th the mass of normal model'
     specificLuminosity = specificLuminosity + 3.*log10(0.6) - log10(512.)
  case (lyAlpha)
     print*, 'lyAlpha model: stellar sources are 0.9286/8th the mass of normal model'
     specificLuminosity = specificLuminosity + log10(65./(70.*8.))
  endselect

!   ! specify sources by physical coordinates

!   nStars = 4
!   allocate(star(nStars))

!   star(1)%data = (/ 0.9765625, 0.0000000, 0.0000000, 151. /)
!   star(2)%data = (/ -4.882812, 22.46094, 30.27344, 21. /)
!   star(3)%data = (/ -2.929688, -0.9765625, -2.929688, 35. /)
!   star(4)%data = (/ -0.9765625, 0.0000000, -0.9765625, 31. /)

!   do iStar = 1, nStars
!      currentStar => star(iStar)
!      xbase = (currentStar%data(1)-xa)/(xb-xa)
!      ybase = (currentStar%data(2)-ya)/(yb-ya)
!      zbase = (currentStar%data(3)-za)/(zb-za)
!      call localizeSplitContinuationCell(xbase,ybase,zbase,0,nx,ny,nz)
!      currentStar%level = splitContinuationCell%level
!      allocate(currentStar%position(3*currentStar%level+3))
!      currentStar%position = splitContinuationCallSequence(1:3*currentStar%level+3)
!      startingPoint = splitContinuation
!      currentStar%luminosity = currentStar%data(4) * 1.d54
!      print*, '--------'
!      print*, 'source #', iStar
!      print*, 'cube =', xbase, ybase, zbase
!      print*, 'level =', currentStar%level
!      print*, 'position =', currentStar%position
!      print*, 'point =', startingPoint
!      print*, 'luminosity =', currentStar%luminosity
!   enddo

!  currentStar%position = (/65,64,65,1,2,1,1,2,1,2,2,2,1,2,2,1,2,1/)

!   ! specify sources by host cells

!   allocate(star(nStars))

!   do iStar = 1, nStars
!      currentStar => star(iStar)
!      select case (iStar)
!      case(1)
!         currentStar%level = 5
!         allocate(currentStar%position(3*currentStar%level+3))
!         currentStar%position = (/65,64,64, 1,2,2, 1,2,2, 1,1,2, 2,1,2, 2,2,2/) ! rho = 3.95296746662E-22
!         currentStar%luminosity = 5.d55 ! total photon # luminosity in all bands [1/s]
!         currentStar%mass = 151. * starParticleMass
!      case(2)
!         currentStar%level = 5
!         allocate(currentStar%position(3*currentStar%level+3))
!         currentStar%position = (/63,64,63, 2,1,2, 2,2,2, 1,2,2, 1,1,1, 1,1,1/) ! rho = 8.34276624451E-22
!         currentStar%luminosity = 1.d30 ! total photon # luminosity in all bands [1/s]
!         currentStar%mass = 35. * starParticleMass
!      case(3)
!         currentStar%level = 4
!         allocate(currentStar%position(3*currentStar%level+3))
!         currentStar%position = (/64,64,64, 2,2,2, 2,2,1, 2,2,2, 2,1,2/) ! rho = 1.39107558854E-22
!         currentStar%luminosity = 1.d30 ! total photon # luminosity in all bands [1/s]
!         currentStar%mass = 31. * starParticleMass
!      case(4)
!         currentStar%level = 4
!         allocate(currentStar%position(3*currentStar%level+3))
!         currentStar%position = (/62,76,80, 1,2,1, 1,1,2, 2,1,1, 1,1,1/) ! rho = 7.53792642535E-23
!         currentStar%luminosity = 1.d30 ! total photon # luminosity in all bands [1/s]
!         currentStar%mass = 21. * starParticleMass
!      end select
!      startingPoint%x = 0.5
!      startingPoint%y = 0.5
!      startingPoint%z = 0.5
!      call absoluteCoordinates(currentStar%level,currentStar%position,startingPoint,nx,ny,nz)
!      print*, '--------'
!      print*, 'source #', iStar
!      print*, 'cube =', xbase, ybase, zbase
!      print*, 'level =', currentStar%level
!      print*, 'position =', currentStar%position
!      print*, 'point =', startingPoint
! !     print*, 'luminosity =', currentStar%luminosity
!      call localizeSplitContinuationCell(xbase,ybase,zbase,0,nx,ny,nz)
!      print*, 'level,rho =', splitContinuationCell%level, splitContinuationCell%rho
!   enddo

  print*, 'counting cells'
  icosmic = 0
  do i = 1, nx
     do j = 1, ny
        do k = 1, nz
           call countCells(baseGrid%cell(i,j,k))
        enddo
     enddo
  enddo
  ncosmic = icosmic

  print*, 'computing ionization equilibrium'
  neutralHydrogenMass = 0.
  totalHydrogenMass = 0.
  icosmic = 0
  do i = 1, nx
     do j = 1, ny
        do k = 1, nz
           do ilambda = 1, 2 ! do it twice: tau=1kpc surface might have shifted
              call initialIonizationEquilibrium(baseGrid%cell(i,j,k))
           enddo
           call computeMass(baseGrid%cell(i,j,k),nx)
        enddo
     enddo
  enddo
  write(*,1040) neutralHydrogenMass/totalHydrogenMass
1040 format('ionization equilibrium:', es18.8)

  print*, 'computing hydrodynamical and local SF heating terms'
  icosmic = 0
  do i = 1, nx
     do j = 1, ny
        do k = 1, nz
           call thermalEquilibrium(baseGrid%cell(i,j,k))
        enddo
     enddo
  enddo

  if (expansionFlag) then
     print*, 'computing expansion of HII regions'
     ! find the lowest density drop factor for each cell
     do iStar = 1, nStars
        currentStar => star(iStar)
        if (currentStar%weight.gt.0) then
           startingPoint%x = 0.5
           startingPoint%y = 0.5
           startingPoint%z = 0.5
           i = currentStar%position(1)
           j = currentStar%position(2)
           k = currentStar%position(3)
           call localizeCellFromStar(baseGrid%cell(i,j,k),currentStar)
           call computeExpansionParameters(psi*currentStar%hostCell%rho/mh)
!           print*, iStar, finalRadius/kpc, densityCoefficient
           call absoluteCoordinates(currentStar%level,currentStar%position,startingPoint,nx,ny,nz)
!           print*, iStar, xbase, ybase, zbase
           do i = 1, nx
              do j = 1, ny
                 do k = 1, nz
                    call findExpansion(baseGrid%cell(i,j,k), &
                         (dfloat(i)-0.5)/dfloat(nx),(dfloat(j)-0.5)/dfloat(ny),(dfloat(k)-0.5)/dfloat(nz),nx)
                 enddo
              enddo
           enddo
        endif
     enddo
     do i = 1, nx
        do j = 1, ny
           do k = 1, nz
              call applyExpansion(baseGrid%cell(i,j,k))
           enddo
        enddo
     enddo
  endif

  if (restart.eq.0) then

     itime = 0

  else

     filename = restartCellArrayName

     itmp = len(trim(filename))
     read (filename(itmp-6:itmp-3), FMT = '(i4)') itime
     print*, 'reading last write from ', trim(filename)

     open(unit=51, file='time', action="write", access="append")
     write(51,*) 'itime =', itime
     close(51)

     sd_id = sfstart(trim(filename), dfacc_read)

     status = sffinfo(sd_id, n_datasets, n_file_attrs)

     start = 0
     edges(1) = 3
     stride = 1

     sds_id = sfselect(sd_id, 0)
     status = sfrdata(sds_id, start, stride, edges, baseGridSize)
     status = sfendacc(sds_id)

     if (nx.ne.baseGridSize(1) .or. ny.ne.baseGridSize(2) .or. &
          nz.ne.baseGridSize(3)) then
        write(*,*) 'error in dimensions'
        print*, nx, ny, nz
        print*, baseGridSize(1), baseGridSize(2), baseGridSize(3)
        stop
     endif

     sds_id = sfselect(sd_id, 1)
     status = sfginfo(sds_id,sds_name,rank,dim_sizes,data_type,n_attrs)
     if (status.ne.0 .or. rank.ne.1) then
        write(*,*) 'error reading ...'
        stop
     endif
     ncosmic = dim_sizes(1)
     edges(1) = dim_sizes(1)

     icosmic = 0
     do i = 1, nx
        do j = 1, ny
           do k = 1, nz
              call countCells(baseGrid%cell(i,j,k))
           enddo
        enddo
     enddo
     if (icosmic.ne.ncosmic) then
        write(*,*) 'error in icosmic/ncosmic ...', icosmic, ncosmic
        stop
     endif

     allocate (cellArrayLevel(icosmic))
     allocate (cellArrayHI(icosmic))
     allocate (cellArrayHeI(icosmic))
     allocate (cellArrayHeII(icosmic))
     allocate (cellArrayTemp(icosmic))

     status = sfrdata(sds_id, start, stride, edges, cellArrayLevel)
     status = sfendacc(sds_id)

     sds_id = sfselect(sd_id, 2)
     status = sfrdata(sds_id, start, stride, edges, cellArrayHI)
     status = sfendacc(sds_id)

     sds_id = sfselect(sd_id, 3)
     status = sfrdata(sds_id, start, stride, edges, cellArrayHeI)
     status = sfendacc(sds_id)

     sds_id = sfselect(sd_id, 4)
     status = sfrdata(sds_id, start, stride, edges, cellArrayHeII)
     status = sfendacc(sds_id)

     sds_id = sfselect(sd_id, 5)
     status = sfrdata(sds_id, start, stride, edges, cellArrayTemp)
     status = sfendacc(sds_id)

     status = sfend(sd_id)

     icosmic = 0
     do i = 1, nx
        do j = 1, ny
           do k = 1, nz
              call readLatestIonization(baseGrid%cell(i,j,k),0)
           enddo
        enddo
     enddo

     deallocate(cellArrayLevel,cellArrayHI,cellArrayHeI,cellArrayHeII,cellArrayTemp)

  endif

  print*, 'finding degenerate (same host cell) particles'

  ! compute unique integer describing location of each particle
  do iStar = 1, nStars
     currentStar => star(iStar)
     i = currentStar%position(1)
     j = currentStar%position(2)
     k = currentStar%position(3)
     currentStar%location = ((i-1)*ny+(j-1))*nz + k
     iref = size(currentStar%position)/3 - 1
     if (iref.ge.1) then
        do itmp = 1, iref
           i = currentStar%position(3*itmp+1)
           j = currentStar%position(3*itmp+2)
           k = currentStar%position(3*itmp+3)
           if (i.eq.1.and.j.eq.1.and.k.eq.1) currentStar%location = 10*currentStar%location + 1
           if (i.eq.1.and.j.eq.1.and.k.eq.2) currentStar%location = 10*currentStar%location + 2
           if (i.eq.1.and.j.eq.2.and.k.eq.1) currentStar%location = 10*currentStar%location + 3
           if (i.eq.1.and.j.eq.2.and.k.eq.2) currentStar%location = 10*currentStar%location + 4
           if (i.eq.2.and.j.eq.1.and.k.eq.1) currentStar%location = 10*currentStar%location + 5
           if (i.eq.2.and.j.eq.1.and.k.eq.2) currentStar%location = 10*currentStar%location + 6
           if (i.eq.2.and.j.eq.2.and.k.eq.1) currentStar%location = 10*currentStar%location + 7
           if (i.eq.2.and.j.eq.2.and.k.eq.2) currentStar%location = 10*currentStar%location + 8
        enddo
     endif
  enddo

  call sortStellarParticles(star)

  ! find degenerate (same host cell) particles
  do iStar = nStars, 2, -1
     currentStar => star(iStar)
     nextStar => star(iStar-1)
     if (currentStar%location.eq.nextStar%location) then
        nextStar%weight = nextStar%weight + currentStar%weight
        currentStar%weight = 0
     endif
  enddo

  nSources = 0
  do iStar = 1, nStars
     if (star(iStar)%weight.gt.0) nSources = nSources + 1
  enddo
  print*, 'nStars/nStarsSpecificAge/non-degenerate =', nStars, nStarsSpecificAge, nSources

  open(25,file='weight',status='replace',form='formatted')
  do iStar = 1, nStars
     currentStar => star(iStar)
     i = currentStar%position(1)
     j = currentStar%position(2)
     k = currentStar%position(3)
     call localizeCellFromStar(baseGrid%cell(i,j,k),currentStar)
     if (currentStar%weight.gt.0) write(25,1058) iStar, currentStar%weight, currentStar%hostCell%abun2
1058 format(i10, ' ==>  ', i10,es16.4)
  enddo
  close(25)

! start iteration to equilibrium

  oldNeutralFraction = 1.d0

  do while (1.lt.2)

     itime = itime + 1

! select stellar population time slice

     timeReadTable = 10.*Myr
     iSpectrum = 1
     do while (timeReadTable.gt.spectrumTime(iSpectrum+1))
        iSpectrum = iSpectrum + 1
     enddo
     coefSpectrum = (timeReadTable-spectrumTime(iSpectrum))/(spectrumTime(iSpectrum+1)-spectrumTime(iSpectrum))
     coefSpectrum = min(max(0.d0,coefSpectrum),1.d0)

! transport

     do i = 1, nx
        do j = 1, ny
           do k = 1, nz
              call setZeroRates(baseGrid%cell(i,j,k))
           enddo
        enddo
     enddo

     maxDepth = 0.

     if (runStellarTransfer) then

        cosmicSpectrum = 0.

        do iStar = 1, nStars

           currentStar => star(iStar)

           if (currentStar%weight.gt.0) then

              highestPixelLevel = 0
              ndotRemaining = 0.
              ndotDust = 0.
              ndotBoundary = 0.
              ndotSpectrum = 0.

              startingPoint%x = 0.5
              startingPoint%y = 0.5
              startingPoint%z = 0.5

              i = currentStar%position(1)
              j = currentStar%position(2)
              k = currentStar%position(3)

              call localizeCellFromStar(baseGrid%cell(i,j,k),currentStar)

              if (currentStar%hostCell%abun2.gt.1.e-20) then
                 tmp = dlog10(currentStar%hostCell%abun2)
              else
                 tmp = - 20.
              endif
              iMetal = 1
              do while (tmp.gt.metallicity(iMetal+1))
                 iMetal = iMetal + 1
                 if (iMetal+1.eq.nMetallicity) exit
              enddo
              coefMetal = (tmp-metallicity(iMetal))/(metallicity(iMetal+1)-metallicity(iMetal))
              coefMetal = min(max(0.d0,coefMetal),1.d0)

              ! compute radiative rates for given stellarPopulation
              ! - this in general depends on iStar and its metallicity

              call stellarBetaTable(nfbins,frequencyBinWidth,totalIntegral, &
                   iSpectrum,coefSpectrum,iMetal,coefMetal)

              if (.not.sphere%refined) then
                 allocate(sphere%pixel(12))
                 sphere%refined = .true.
              endif

              ndot1 = float(currentStar%weight) ! in units of stellarPopulation(iStar,iSpectrum,coefSpectrum,nu1)

              do iray = 1, 12

                 leafPixel => sphere%pixel(iray)
                 if (leafPixel%level.ne.1) then
                    leafPixel%level = 1
                    leafPixel%refined = .false.
                    leafPixel%parent => sphere
                    nside = 2**(leafPixel%level-1)
                    longIntegerArgument = iray - 1
                    call pix2ang_nest(nside,longIntegerArgument,leafPixel%phi,leafPixel%theta)
                 endif

                 depth1 = 0.
                 depth2 = 0.
                 depth3 = 0.
                 depthDust = 0.

                 call startNewLongRay(currentStar%hostCell,startingPoint,leafPixel,iray, &
                      currentStar%level,currentStar%position,0.d0,ndot1/12.d0, &
                      depth1,depth2,depth3,depthDust,nx,ny,nz)

              enddo

!         ! clean up segments
!         icosmic = 0
!         do i = 1, nx
!            do j = 1, ny
!               do k = 1, nz
!                  call cleanUpSegments(baseGrid%cell(i,j,k))
!               enddo
!            enddo
!         enddo
!         stop

              do iradius = 1, nradius
                 if (ndotBoundary(iradius).lt.1.) then
                    fraction(iradius) = ndotRemaining(iradius)/(ndot1-ndotBoundary(iradius))
                 else
                    fraction(iradius) = 0.
                 endif
              enddo

              cosmicSpectrum = cosmicSpectrum + &
                   float(currentStar%weight) * ndotSpectrum/(ndot1-ndotBoundary(nradius))

              write(*,1015) iStar, currentStar%level, &
                   currentStar%hostCell%HI * mh / (psi * currentStar%hostCell%rho), &
                   highestPixelLevel, fraction, currentStar%weight

1015          format('src: ', i5, i3, es13.5, i3, 7f9.5, i8)

!                  if (ndotBoundary(nradius).lt.1.) then
!                     print*, 'dust =', ndotDust/(ndot1-ndotBoundary(nradius))
!                  endif

           endif
        enddo

        cosmicSpectrum = cosmicSpectrum / float(nStarsSpecificAge)

1025    format('esc: ', i4, f11.3, es15.5)

     endif ! runStellarTransfer

     if (runUVBTransfer) then

        print*, 'computing opacities'
        do i = 1, nx
           do j = 1, ny
              do k = 1, nz
                 call computeOpacities(baseGrid%cell(i,j,k))
              enddo
           enddo
        enddo

        ! set up the rays

        nside = 2**(nAngularLevel-1)
        weight = 1./float(12*4**(nAngularLevel-1))

        print*, 'starting transport'
        do iray = 0, 12*4**(nAngularLevel-1)-1

           call pix2ang_nest(nside,iray,phiLarge,thetaLarge)

           ! phi goes from 0 to 2*pi; theta goes from -pi/2 to pi/2

           angles%izone = 1
           if (phiLarge.gt.0. .and. phiLarge.lt.0.5*pi) then
              phi1 = phiLarge
              angles%izone = angles%izone + 0
           else
              if (phiLarge.gt.0.5*pi .and. phiLarge.lt.pi) then
                 phi1 = phiLarge - 0.5*pi
                 angles%izone = angles%izone + 3
              else
                 if (phiLarge.gt.pi .and. phiLarge.lt.1.5*pi) then
                    phi1 = phiLarge - pi
                    angles%izone = angles%izone + 6
                 else
                    if (phiLarge.gt.1.5*pi .and. phiLarge.lt.2.*pi) then
                       phi1 = phiLarge - 1.5*pi
                       angles%izone = angles%izone + 9
                    else
                       write(*,*) 'error in phi', phiLarge, thetaLarge
                       stop
                    endif
                 endif
              endif
           endif
           if (thetaLarge.gt.0. .and. thetaLarge.lt.0.5*pi) then
              theta1 = thetaLarge
              angles%izone = angles%izone + 0
           else
              if (thetaLarge.gt.-0.5*pi .and. thetaLarge.lt.0.) then
                 theta1 = - thetaLarge
                 angles%izone = angles%izone + 12
              else
                 write(*,*) 'error in theta', phiLarge, thetaLarge
                 stop
              endif
           endif

           tmp1 = 1./sin(theta1)
           tmp2 = 1./(cos(phi1)*cos(theta1))
           tmp3 = 1./(sin(phi1)*cos(theta1))

           if (tmp1.lt.min(tmp2,tmp3)) then
              angles%theta = theta1
              angles%phi = phi1
              angles%izone = angles%izone + 0
           else
              if (tmp2.lt.min(tmp1,tmp3)) then
                 angles%theta = arcsin(cos(theta1)*cos(phi1))
                 angles%phi = arcsin(sin(theta1)/cos(angles%theta))
                 angles%izone = angles%izone + 1
              else
                 if (tmp3.lt.min(tmp1,tmp2)) then
                    angles%theta = arcsin(cos(theta1)*sin(phi1))
                    angles%phi = acos(sin(theta1)/cos(angles%theta))
                    angles%izone = angles%izone + 2
                 else
                    write(*,*) 'error in theta or phi'
                    stop
                 endif
              endif
           endif

!           write(*,*) 'izone =', angles%izone, angles%phi, angles%theta

           select case (angles%izone)
           case(1,7,13,19)
              nxtransfer = nx
              nytransfer = ny
              nztransfer = nz
           case(2,8,14,20)
              nxtransfer = ny
              nytransfer = nz
              nztransfer = nx
           case(3,9,15,21)
              nxtransfer = nz
              nytransfer = nx
              nztransfer = ny
           case(4,10,16,22)
              nxtransfer = nx
              nytransfer = nz
              nztransfer = ny
           case(5,11,17,23)
              nxtransfer = ny
              nytransfer = nx
              nztransfer = nz
           case(6,12,18,24)
              nxtransfer = nz
              nytransfer = ny
              nztransfer = nx
           end select

           do i = 1, 2
              do j = 1, 2
                 do k = 1, 2
                    call rotateIndices(i,j,k,2,2,2,angles%izone,is(i,j,k),js(i,j,k),ks(i,j,k))
                 enddo
              enddo
           enddo

           allocate(pattern(nxtransfer))

           do i = 1, nxtransfer

              if (i.eq.1) then

                 pattern(i)%xyRay%x0 = 0.5
                 pattern(i)%xyRay%y0 = 0.5
                 call setPattern(pattern(i),angles%phi,angles%theta)

!           call checkPattern(pattern(i),angles%phi,angles%theta)

              else

                 select case (pattern(i-1)%xyTop)
                 case(xyEnd) ! use pattern(i-1)%xyRay
                    pattern(i)%xyRay%x0 = pattern(i-1)%xyRay%x0 + cos(angles%phi)/tan(angles%theta)
                    pattern(i)%xyRay%y0 = pattern(i-1)%xyRay%y0 + sin(angles%phi)/tan(angles%theta)
                 case(xzEnd) ! use pattern(i-1)%xzRay
                    pattern(i)%xyRay%x0 = pattern(i-1)%xzRay%x0 + &
                         pattern(i-1)%xzRay%len*cos(angles%theta)*cos(angles%phi)
                    pattern(i)%xyRay%y0 = pattern(i-1)%xzRay%len*cos(angles%theta)*sin(angles%phi)
                 case(yzEnd) ! use pattern(i-1)%yzRay
                    pattern(i)%xyRay%x0 = pattern(i-1)%yzRay%len*cos(angles%theta)*cos(angles%phi)
                    pattern(i)%xyRay%y0 = pattern(i-1)%yzRay%y0 + &
                         pattern(i-1)%yzRay%len*cos(angles%theta)*sin(angles%phi)
                 case(0)
                    write(*,*) 'error in xyTop'
                    stop
                 end select
                 if (pattern(i)%xyRay%x0.gt.1. .or. pattern(i)%xyRay%y0.gt.1.) then
                    write(*,*) '1) Error: xyRay%x0, y0', pattern(i)%xyRay%x0, pattern(i)%xyRay%y0
                    write(*,*) pattern(i-1)%xyRay%x0
                    stop
                 endif
                 call setPattern(pattern(i),angles%phi,angles%theta)

!           call checkPattern(pattern(i),angles%phi,angles%theta)

              endif

              pattern(i)%refined = .false.

              do j = 1, nytransfer
                 do k = 1, nztransfer

                    call rotateIndices(i,j,k,nx,ny,nz,angles%izone,icell,jcell,kcell)

                    currentCell => baseGrid%cell(icell,jcell,kcell)

                    ! call copyPatternToCell(pattern(i),currentCell)
                    currentCell%pattern => pattern(i)
                    if (currentCell%refined) call setRaysRefined(currentCell,pattern(i), &
                         is,js,ks,angles%phi,angles%theta)

                    currentCell%parent => baseGrid

                 enddo
              enddo

           enddo

! find neighbours

           do i = 1, nxtransfer
              do j = 1, nytransfer
                 do k = 1, nztransfer
                    call rotateIndices(i,j,k,nx,ny,nz,angles%izone,icell,jcell,kcell)
!              write(*,'("entered cell",3i4,"   of level",i5)') icell, jcell, kcell, 0
                    call localizeCellFindNeighbours(baseGrid%cell(icell,jcell,kcell), &
                         0,(/i,j,k/),is,js,ks,nx,ny,nz,angles%izone)
                 enddo
              enddo
           enddo

! transport

           cellSizeAbsoluteUnits = physicalBoxSize/dfloat(nx)

           do i = 1, nxtransfer
              do j = 1, nytransfer
                 do k = 1, nztransfer

                    call rotateIndices(i,j,k,nx,ny,nz,angles%izone,icell,jcell,kcell)

                    currentCell => baseGrid%cell(icell,jcell,kcell)

                    if (.not.currentCell%refined) then

                       Jmean1 = 0.
                       Jmean2 = 0.
                       Jmean3 = 0.
                       imean = 0

                       if (.not.currentCell%xyNeighbourPresent) then
                          Iin1 = uvb1
                          Iin2 = uvb2
                          Iin3 = uvb3
                       else
                          select case (currentCell%xyNeighbour%pattern%xyTop)
                          case(xyEnd)
                             Iin1 = currentCell%xyNeighbour%rt%xyRay%Iout1
                             Iin2 = currentCell%xyNeighbour%rt%xyRay%Iout2
                             Iin3 = currentCell%xyNeighbour%rt%xyRay%Iout3
                          case(xzEnd)
                             Iin1 = currentCell%xyNeighbour%rt%xzRay%Iout1
                             Iin2 = currentCell%xyNeighbour%rt%xzRay%Iout2
                             Iin3 = currentCell%xyNeighbour%rt%xzRay%Iout3
                          case(yzEnd)
                             Iin1 = currentCell%xyNeighbour%rt%yzRay%Iout1
                             Iin2 = currentCell%xyNeighbour%rt%yzRay%Iout2
                             Iin3 = currentCell%xyNeighbour%rt%yzRay%Iout3
                          case(0)
                             write(*,*) 'error in xyTop'
                             stop
                          end select
                       endif

                       dpath = cellSizeAbsoluteUnits * currentCell%pattern%xyRay%len ! [cm]
                       tau1 = currentCell%kappa1 * dpath
                       tau2 = currentCell%kappa2 * dpath
                       tau3 = currentCell%kappa3 * dpath
                       tmpabs1 = exp(-tau1)
                       tmpabs2 = exp(-tau2)
                       tmpabs3 = exp(-tau3)
                       if (tau1.gt.1.e-10) then
                          tmpemi1 = (1.-tmpabs1)/currentCell%kappa1
                       else
                          tmpemi1 = dpath
                       endif
                       if (tau2.gt.1.e-10) then
                          tmpemi2 = (1.-tmpabs2)/currentCell%kappa2
                       else
                          tmpemi2 = dpath
                       endif
                       if (tau3.gt.1.e-10) then
                          tmpemi3 = (1.-tmpabs3)/currentCell%kappa3
                       else
                          tmpemi3 = dpath
                       endif
                       nemi1 = 0. ! emissivity [erg/s/cm^2/Hz/ster]
                       nemi2 = 0.
                       nemi3 = 0.
                       currentCell%rt%xyRay%Iout1 = Iin1*tmpabs1 + nemi1*tmpemi1/dpath
                       currentCell%rt%xyRay%Iout2 = Iin2*tmpabs2 + nemi2*tmpemi2/dpath
                       currentCell%rt%xyRay%Iout3 = Iin3*tmpabs3 + nemi3*tmpemi3/dpath

                       call computeCellIntensity(Jmean1,Iin1,currentCell%rt%xyRay%Iout1)
                       call computeCellIntensity(Jmean2,Iin2,currentCell%rt%xyRay%Iout2)
                       call computeCellIntensity(Jmean3,Iin3,currentCell%rt%xyRay%Iout3)
                       imean = imean + 1

                       if (currentCell%pattern%xzRayActive) then
                          if (.not.currentCell%xzNeighbourPresent) then
                             Iin1 = uvb1
                             Iin2 = uvb2
                             Iin3 = uvb3
                          else
                             select case (currentCell%xzNeighbour%pattern%xzTop)
                             case(xyEnd)
                                Iin1 = currentCell%xzNeighbour%rt%xyRay%Iout1
                                Iin2 = currentCell%xzNeighbour%rt%xyRay%Iout2
                                Iin3 = currentCell%xzNeighbour%rt%xyRay%Iout3
                             case(xzEnd)
                                if (.not.currentCell%xzNeighbour%pattern%xzRayActive) then
                                   write(*,*) '1) Error: xzRay should be active ...', i, j, k
                                   stop
                                endif
                                Iin1 = currentCell%xzNeighbour%rt%xzRay%Iout1
                                Iin2 = currentCell%xzNeighbour%rt%xzRay%Iout2
                                Iin3 = currentCell%xzNeighbour%rt%xzRay%Iout3
                             case(yzEnd)
                                if (.not.currentCell%xzNeighbour%pattern%yzRayActive) then
                                   write(*,*) '1) Error: yzRay should be active ...', i, j, k
                                   stop
                                endif
                                Iin1 = currentCell%xzNeighbour%rt%yzRay%Iout1
                                Iin2 = currentCell%xzNeighbour%rt%yzRay%Iout2
                                Iin3 = currentCell%xzNeighbour%rt%yzRay%Iout3
                             case(0)
                                write(*,*) 'error in xzTop'
                                stop
                             end select
                          endif

                          dpath = cellSizeAbsoluteUnits * currentCell%pattern%xzRay%len ! [cm]
                          tau1 = currentCell%kappa1 * dpath
                          tau2 = currentCell%kappa2 * dpath
                          tau3 = currentCell%kappa3 * dpath
                          tmpabs1 = exp(-tau1)
                          tmpabs2 = exp(-tau2)
                          tmpabs3 = exp(-tau3)
                          if (tau1.gt.1.e-10) then
                             tmpemi1 = (1.-tmpabs1)/currentCell%kappa1
                          else
                             tmpemi1 = dpath
                          endif
                          if (tau2.gt.1.e-10) then
                             tmpemi2 = (1.-tmpabs2)/currentCell%kappa2
                          else
                             tmpemi2 = dpath
                          endif
                          if (tau3.gt.1.e-10) then
                             tmpemi3 = (1.-tmpabs3)/currentCell%kappa3
                          else
                             tmpemi3 = dpath
                          endif
                          nemi1 = 0. ! emissivity [erg/s/cm^2/Hz/ster]
                          nemi2 = 0.
                          nemi3 = 0.
                          currentCell%rt%xzRay%Iout1 = Iin1*tmpabs1 + nemi1*tmpemi1/dpath
                          currentCell%rt%xzRay%Iout2 = Iin2*tmpabs2 + nemi2*tmpemi2/dpath
                          currentCell%rt%xzRay%Iout3 = Iin3*tmpabs3 + nemi3*tmpemi3/dpath

                          call computeCellIntensity(Jmean1,Iin1,currentCell%rt%xzRay%Iout1)
                          call computeCellIntensity(Jmean2,Iin2,currentCell%rt%xzRay%Iout2)
                          call computeCellIntensity(Jmean3,Iin3,currentCell%rt%xzRay%Iout3)
                          imean = imean + 1

                       endif

                       if (currentCell%pattern%yzRayActive) then
                          if (.not.currentCell%yzNeighbourPresent) then
                             Iin1 = uvb1
                             Iin2 = uvb2
                             Iin3 = uvb3
                          else
                             select case (currentCell%yzNeighbour%pattern%yzTop)
                             case(xyEnd)
                                Iin1 = currentCell%yzNeighbour%rt%xyRay%Iout1
                                Iin2 = currentCell%yzNeighbour%rt%xyRay%Iout2
                                Iin3 = currentCell%yzNeighbour%rt%xyRay%Iout3
                             case(xzEnd)
                                if (.not.currentCell%yzNeighbour%pattern%xzRayActive) then
                                   write(*,*) '2) Error: xzRay should be active ...', i, j, k
                                   stop
                                endif
                                Iin1 = currentCell%yzNeighbour%rt%xzRay%Iout1
                                Iin2 = currentCell%yzNeighbour%rt%xzRay%Iout2
                                Iin3 = currentCell%yzNeighbour%rt%xzRay%Iout3
                             case(yzEnd)
                                if (.not.currentCell%yzNeighbour%pattern%yzRayActive) then
                                   write(*,*) '2) Error: yzRay should be active ...', i, j, k
                                   stop
                                endif
                                Iin1 = currentCell%yzNeighbour%rt%yzRay%Iout1
                                Iin2 = currentCell%yzNeighbour%rt%yzRay%Iout2
                                Iin3 = currentCell%yzNeighbour%rt%yzRay%Iout3
                             case(0)
                                write(*,*) '1) error in yzTop'
                                write(*,*) currentCell%pattern%xyRay%x0, currentCell%pattern%xyRay%y0
                                write(*,*) currentCell%yzNeighbour%pattern%xyRay%x0, &
                                     currentCell%yzNeighbour%pattern%xyRay%y0
                                stop
                             end select
                          endif

                          dpath = cellSizeAbsoluteUnits * currentCell%pattern%yzRay%len ! [cm]
                          tau1 = currentCell%kappa1 * dpath
                          tau2 = currentCell%kappa2 * dpath
                          tau3 = currentCell%kappa3 * dpath
                          tmpabs1 = exp(-tau1)
                          tmpabs2 = exp(-tau2)
                          tmpabs3 = exp(-tau3)
                          if (tau1.gt.1.e-10) then
                             tmpemi1 = (1.-tmpabs1)/currentCell%kappa1
                          else
                             tmpemi1 = dpath
                          endif
                          if (tau2.gt.1.e-10) then
                             tmpemi2 = (1.-tmpabs2)/currentCell%kappa2
                          else
                             tmpemi2 = dpath
                          endif
                          if (tau3.gt.1.e-10) then
                             tmpemi3 = (1.-tmpabs3)/currentCell%kappa3
                          else
                             tmpemi3 = dpath
                          endif
                          nemi1 = 0. ! emissivity [erg/s/cm^2/Hz/ster]
                          nemi2 = 0.
                          nemi3 = 0.
                          currentCell%rt%yzRay%Iout1 = Iin1*tmpabs1 + nemi1*tmpemi1/dpath
                          currentCell%rt%yzRay%Iout2 = Iin2*tmpabs2 + nemi2*tmpemi2/dpath
                          currentCell%rt%yzRay%Iout3 = Iin3*tmpabs3 + nemi3*tmpemi3/dpath

                          call computeCellIntensity(Jmean1,Iin1,currentCell%rt%yzRay%Iout1)
                          call computeCellIntensity(Jmean2,Iin2,currentCell%rt%yzRay%Iout2)
                          call computeCellIntensity(Jmean3,Iin3,currentCell%rt%yzRay%Iout3)
                          imean = imean + 1

                       endif

                       currentCell%Jmean1 = currentCell%Jmean1 + Jmean1/float(imean)*weight
                       currentCell%Jmean2 = currentCell%Jmean2 + Jmean2/float(imean)*weight
                       currentCell%Jmean3 = currentCell%Jmean3 + Jmean3/float(imean)*weight

                    else
                       call transport(currentCell,weight,is,js,ks,cellSizeAbsoluteUnits)
                    endif

                 enddo
              enddo
           enddo

           do i = 1, nxtransfer
              if (pattern(i)%refined) call patternNullify(pattern(i))
           enddo
           deallocate(pattern)

           write(*,*) '>>>', float(iray+1)/float(12*4**(nAngularLevel-1)), &
                baseGrid%cell(1,1,1)%Jmean1

        enddo ! end of loop over all rays

     endif ! runUVBTransfer

! chemistry part

     do i = 1, nx
        do j = 1, ny
           do k = 1, nz
              call saveRestoreOriginalFields(baseGrid%cell(i,j,k))
           enddo
        enddo
     enddo

     icosmic = 0
     neutralHydrogenMass = 0.
     totalHydrogenMass = 0.

     do i = 1, nx
        do j = 1, ny
           do k = 1, nz
              call solveRateEquations(baseGrid%cell(i,j,k),nx,runUVBTransfer)
              call computeMass(baseGrid%cell(i,j,k),nx)
           enddo
        enddo
     enddo

     open(unit=51, file='time', action="write", access="append")
     write(51,1034) itime, neutralHydrogenMass/totalHydrogenMass
     close(51)
1034 format('itime =', i5, f18.10)

!     currentCell => baseGrid%cell(64,63,61)
!     print*, 'rad', currentCell%HI * mh / (psi * currentCell%rho)

     call writeIonization(itime,nx,ny,nz)

  enddo

  do i = 1, nx
     do j = 1, ny
        do k = 1, nz
           call scanCube(baseGrid%cell(i,j,k),0,(/i,j,k/),nx,ny,nz)
        enddo
     enddo
  enddo

  ! currentCell => baseGrid%cell(1,1,1)
  ! call analyzeCell(currentCell,nx)
  ! print*, 'level =', currentCell%level

  do i = 1, nx
     do j = 1, nx
        do k = 1, nx
           call analyzeCell(baseGrid%cell(i,j,k),nx)
        enddo
     enddo
  enddo

  deallocate(baseGrid%cell)
  stop

contains

  recursive subroutine placeCellProjectWithVelocity(parentCell,level,x0,y0,z0,readArray)

    use definitions

    implicit none
    integer, intent(in) :: level
    integer :: inew, jnew, knew, i, j, k, nvariables
    real(kind=RealKind), intent(in) :: x0, y0, z0
    real*4, dimension(maxNumberReadVariables), intent(in) :: readArray
    real(kind=RealKind) :: xnew, ynew, znew, nh, xneu, nhe
    type(zoneType), target :: parentCell

    if (level.gt.1) then

       if (.not.parentCell%refined) then
          allocate(parentCell%cell(2,2,2))
          parentCell%refined = .true.
          do i = 1, 2
             do j = 1, 2
                do k = 1, 2
                   parentCell%cell(i,j,k)%refined = .false.
                   parentCell%cell(i,j,k)%level = parentCell%level + 1
                   parentCell%cell(i,j,k)%tgas = parentCell%tgas
                   parentCell%cell(i,j,k)%rho = parentCell%rho
                   parentCell%cell(i,j,k)%HI = parentCell%HI
                   parentCell%cell(i,j,k)%HeI = parentCell%HeI
                   parentCell%cell(i,j,k)%HeII = parentCell%HeII
                   parentCell%cell(i,j,k)%rhoCoef = 1.
                   parentCell%cell(i,j,k)%velx = 0.d0
                   parentCell%cell(i,j,k)%vely = 0.d0
                   parentCell%cell(i,j,k)%velz = 0.d0
                   parentCell%cell(i,j,k)%abun2 = 0.d0
                   parentCell%cell(i,j,k)%parent => parentCell
                   nullify(parentCell%cell(i,j,k)%cell)
                   nullify(parentCell%cell(i,j,k)%segment)
                   nullify(parentCell%cell(i,j,k)%lastSegment)
                enddo
             enddo
          enddo
       endif
       if (x0.lt.0.5) then
          inew = 1
          xnew = 2.*x0
       else
          inew = 2
          xnew = 2.*x0 - 1.
       endif
       if (y0.lt.0.5) then
          jnew = 1
          ynew = 2.*y0
       else
          jnew = 2
          ynew = 2.*y0 - 1.
       endif
       if (z0.lt.0.5) then
          knew = 1
          znew = 2.*z0
       else
          knew = 2
          znew = 2.*z0 - 1.
       endif
       call placeCellProjectWithVelocity(parentCell%cell(inew,jnew,knew),level-1, &
            xnew,ynew,znew,readArray)
    else

       parentCell%tgas = 10.**readArray(1)
       nh = 10.**readArray(2)
       !     nh = 1.e-5 ! uniform for testing
       xneu = 10.**readArray(3)
       parentCell%rho = nh * mh/psi
       parentCell%HI = nh * xneu
       nhe = (1.-psi) * parentCell%rho / mhe
       parentCell%HeI = nhe * 1.
       parentCell%HeII = nhe * 0.
       parentCell%rhoCoef = 1.
       nvariables = 3

       if (readKinematics) then
          parentCell%velx = readArray(nvariables+1)
          parentCell%vely = readArray(nvariables+2)
          parentCell%velz = readArray(nvariables+3)
          nvariables = nvariables + 3
       endif

       if (readMetals) then
          parentCell%abun2 = readArray(nvariables+2)
          nvariables = nvariables + 4
       else
          parentCell%abun2 = 0.02
       endif

       !      print*, 'rho,temp =', parentCell%rho, parentCell%tgas
       !      print*, 'HI,xneu =', parentCell%HI, xneu
       !      print*, 'HeI,HeII =', parentCell%HeI, parentCell%HeII
       !      stop

       parentCell%krate24 = 0.
       parentCell%krate25 = 0.
       parentCell%krate26 = 0.
       parentCell%crate24 = 0.
       parentCell%crate25 = 0.
       parentCell%crate26 = 0.
    endif

  end subroutine placeCellProjectWithVelocity

  recursive subroutine placeCellProject(parentCell,level,x0,y0,z0,logtgas,lognh,logxneu)

    use definitions

    implicit none
    integer, intent(in) :: level
    integer :: inew, jnew, knew, i, j, k
    real(kind=RealKind), intent(in) :: x0, y0, z0
    real*4, intent(in) :: logtgas, lognh, logxneu
    real(kind=RealKind) :: xnew, ynew, znew, nh, xneu, nhe
    type(zoneType), target :: parentCell

    !  write(*,*) x0, y0, z0, level

    if (level.gt.1) then

       if (.not.parentCell%refined) then
          allocate(parentCell%cell(2,2,2))
          parentCell%refined = .true.
          do i = 1, 2
             do j = 1, 2
                do k = 1, 2
                   parentCell%cell(i,j,k)%refined = .false.
                   parentCell%cell(i,j,k)%level = parentCell%level + 1
                   parentCell%cell(i,j,k)%tgas = parentCell%tgas
                   parentCell%cell(i,j,k)%rho = parentCell%rho
                   parentCell%cell(i,j,k)%HI = parentCell%HI
                   parentCell%cell(i,j,k)%HeI = parentCell%HeI
                   parentCell%cell(i,j,k)%HeII = parentCell%HeII
                   parentCell%cell(i,j,k)%rhoCoef = 1.
                   parentCell%cell(i,j,k)%parent => parentCell
                   nullify(parentCell%cell(i,j,k)%cell)
                   nullify(parentCell%cell(i,j,k)%segment)
                   nullify(parentCell%cell(i,j,k)%lastSegment)
                enddo
             enddo
          enddo
       endif
       if (x0.lt.0.5) then
          inew = 1
          xnew = 2.*x0
       else
          inew = 2
          xnew = 2.*x0 - 1.
       endif
       if (y0.lt.0.5) then
          jnew = 1
          ynew = 2.*y0
       else
          jnew = 2
          ynew = 2.*y0 - 1.
       endif
       if (z0.lt.0.5) then
          knew = 1
          znew = 2.*z0
       else
          knew = 2
          znew = 2.*z0 - 1.
       endif
       call placeCellProject(parentCell%cell(inew,jnew,knew),level-1, &
            xnew,ynew,znew,logtgas,lognh,logxneu)
    else
       parentCell%tgas = 10.**logtgas
       nh = 10.**lognh
       !     nh = 1.e-5 ! uniform for testing
       xneu = 10.**logxneu
       parentCell%rho = nh * mh/psi
       parentCell%HI = nh * xneu
       nhe = (1.-psi) * parentCell%rho / mhe
       parentCell%HeI = nhe * 1.
       parentCell%HeII = nhe * 0.
       parentCell%rhoCoef = 1.

       !      print*, 'rho,temp =', parentCell%rho, parentCell%tgas
       !      print*, 'HI,xneu =', parentCell%HI, xneu
       !      print*, 'HeI,HeII =', parentCell%HeI, parentCell%HeII
       !      stop

       parentCell%krate24 = 0.
       parentCell%krate25 = 0.
       parentCell%krate26 = 0.
       parentCell%crate24 = 0.
       parentCell%crate25 = 0.
       parentCell%crate26 = 0.
    endif

  end subroutine placeCellProject

  recursive subroutine printCell(currentCell,level,callSequence)

    use definitions

    implicit none
    type(zoneType), target :: currentCell
    integer, intent(in) :: level
    integer, intent(in) :: callSequence(3*level+3)
    integer :: i, j, k, callLength
    integer :: newCallSequence(3*level+6)
    real(kind=RealKind) :: xneu

    callLength = 3*level + 3

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                newCallSequence(1:callLength) = callSequence
                newCallSequence(callLength+1:callLength+3) = (/i,j,k/)
                call printCell(currentCell%cell(i,j,k),level+1,newCallSequence)
             enddo
          enddo
       enddo
    else
       xneu = currentCell%HI * mh / (psi * currentCell%rho)
       select case (level)
       case(0)
          write(*,'(3i4,es18.11,f18.11,"     l=", i1)') callSequence, currentCell%rho, xneu, level
       case(1)
          write(*,'(6i4,es18.11,f18.11,"     l=", i1)') callSequence, currentCell%rho, xneu, level
       case(2)
          write(*,'(9i4,es18.11,f18.11,"     l=", i1)') callSequence, currentCell%rho, xneu, level
       case(3)
          write(*,'(12i4,es18.11,f18.11,"     l=", i1)') callSequence, currentCell%rho, xneu, level
       case(4)
          write(*,'(15i4,es18.11,f18.11,"     l=", i1)') callSequence, currentCell%rho, xneu, level
       case(5)
          write(*,'(18i4,es18.11,f18.11,"     l=", i1)') callSequence, currentCell%rho, xneu, level
       case(6)
          write(*,'(21i4,es18.11,f18.11,"     l=", i1)') callSequence, currentCell%rho, xneu, level
       case(7)
          write(*,'(24i4,es18.11,f18.11,"     l=", i1)') callSequence, currentCell%rho, xneu, level
       case(8)
          write(*,'(27i4,es18.11,f18.11,"     l=", i1)') callSequence, currentCell%rho, xneu, level
       end select
       if (currentCell%rho.gt.densestCell%rho) then
          densestCell => currentCell
          write(*,*) 'new densest cell found'
       endif
    endif

  end subroutine printCell

  subroutine pix2ang_nest(nside,ipix,phi,theta)

    !=======================================================================
    !     renders theta and phi coordinates of the nominal pixel center
    !     for the pixel number ipix (NESTED scheme)  
    !     given the map resolution parameter nside
    !=======================================================================

    use definitions

    implicit none
    integer, intent(in) :: nside
    integer*8, intent(in) :: ipix
    real(kind=RealKind) :: theta, phi, rotationAngle

    integer*8 :: nsideLong, npix, npface, ipf, ip_low, ip_trunc, ip_med
    integer ip_hi, &
         jrt, jr, nr, jpt, jp, kshift, nl4
    real(kind=RealKind) :: z, fn, fact1, fact2
    integer, parameter :: ns_max = 8192*4
    integer pix2x(0:1023), pix2y(0:1023)
    common/nested_block/ pix2x, pix2y

    integer ix, iy, face_num
    common /xy_nest/ ix, iy, face_num ! can be useful to calling routine

    integer jrll(12), jpll(12) ! coordinate of the lowest corner of each face
    data jrll/2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4/ ! in unit of nside
    data jpll/1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7/ ! in unit of nside/2

    !  print*, '>>>>>>', ipix

    nsideLong = nside

    if (nside.lt.1 .or. nside.gt.ns_max) then
       print*, 'nside out of range', nside
       stop
    endif
    npix = 12 * nsideLong**2
    if (ipix .lt.0 .or. ipix.gt.npix-1) then
       print*, 'ipix out of range', ipix, npix, nside
       stop
    endif

    ! initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) .le. 0) call mk_pix2xy

    fn = float(nside)
    fact1 = 1./(3.*fn*fn)
    fact2 = 2./(3.*fn)
    nl4   = 4*nside

    ! finds the face, and the number in the face
    npface = nsideLong**2

    face_num = ipix/npface    ! face number in {0,11}
    ipf = MOD(ipix,npface)    ! pixel number in the face {0,npface-1}

    ! finds the x,y on the face (starting from the lowest corner)
    ! from the pixel number
    ip_low = MOD(ipf,1024_LongInteger)    ! content of the last 10 bits
    ip_trunc = ipf/1024_LongInteger     ! truncation of the last 10 bits
    ip_med = MOD(ip_trunc,1024_LongInteger) ! content of the next 10 bits
    ip_hi = ip_trunc/1024_LongInteger ! content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)

    ! transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy             ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy             ! 'horizontal' in {-nside+1,nside-1}

    ! computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1 ! ring number in {1,4*nside-1}

    nr = nside                ! equatorial region (the most frequent)
    z  = float(2*nside-jr)*fact2
    kshift = MOD(jr - nside, 2)
    if (jr .lt. nside) then   ! north pole region
       nr = jr
       z = 1. - float(nr)*float(nr)*fact1
       kshift = 0
    else if (jr .gt. 3*nside) then ! south pole region
       nr = nl4 - jr
       z = - 1. + float(nr)*float(nr)*fact1
       kshift = 0
    endif
    theta = acos(z) - halfPi

    ! computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2 ! 'phi' number in the ring in {1,4*nr}
    if (jp .gt. nl4) jp = jp - nl4
    if (jp .lt. 1)   jp = jp + nl4

    phi = (float(jp)-float(kshift+1)*0.5)*halfPi/float(nr)

    do while (phi.gt.twoPi)
       phi = phi - twoPi
    end do

    do while (phi.lt.0.)
       phi = phi + twoPi
    end do

    call rotateAngles(phi,theta)

    if (phi.gt.2.*Pi) then
       write(*,*) 'angle too large'
       stop
    endif

    return

  end subroutine pix2ang_nest

  subroutine mk_pix2xy

    !=======================================================================
    !     constructs the array giving x and y in the face from pixel number
    !     for the nested (quad-cube like) ordering of pixels
    !
    !     the bits corresponding to x and y are interleaved in the pixel number
    !     one breaks up the pixel number by even and odd bits
    !=======================================================================

    IMPLICIT none
    INTEGER kpix, jpix, ix, iy, ip, id
    integer pix2x(0:1023), pix2y(0:1023)
    common/nested_block/ pix2x, pix2y

    !cc cf block data      data      pix2x(1023) /0/
    !-----------------------------------------------------------------------
    !      print *, 'initiate pix2xy'

    do kpix=0,1023            ! pixel number
       jpix = kpix
       ix = 0
       iy = 0
       ip = 1                  ! bit position (in x and y)
       do while (jpix.ne.0)    ! go through all the bits
          id = mod(jpix,2)      ! bit value (in kpix), goes in ix
          jpix = jpix/2
          ix = id*ip+ix

          id = mod(jpix,2)      ! bit value (in kpix), goes in iy
          jpix = jpix/2
          iy = id*ip+iy

          ip = 2*ip             ! next bit (in x and y)
       enddo
       pix2x(kpix) = ix        ! in 0,31
       pix2y(kpix) = iy        ! in 0,31

    enddo

    return

  end subroutine mk_pix2xy

  function arcsin(x) result(angle)

    use definitions

    implicit none
    real(kind=RealKind), intent(in) :: x
    real(kind=RealKind) :: angle

    if (x.gt.1.d0) then
       angle = halfPi
    else
       if (x.lt.-1.d0) then
          angle = - halfPi
       else
          angle = asin(x)
       endif
    endif

  end function arcsin

  subroutine rotateAngles(phi,theta)

    use definitions

    implicit none
    real(kind=RealKind), intent(inout) :: phi, theta
    real(kind=RealKind) :: phi0, theta0, cosphi, sinphi, rotationAngle

    ! rotation around x-axis

    phi0 = phi
    theta0 = theta
    rotationAngle = 0.111
    theta = arcsin(cos(theta0)*sin(phi0)*sin(rotationAngle) + sin(theta0)*cos(rotationAngle))
    cosphi = cos(theta0)*cos(phi0)/cos(theta)
    sinphi = (cos(theta0)*sin(phi0)*cos(rotationAngle) - sin(theta0)*sin(rotationAngle))/cos(theta)
    phi = getAngle(cosphi,sinphi)

    ! rotation around y-axis

    phi0 = phi
    theta0 = theta
    rotationAngle = 0.222
    theta = arcsin(cos(theta0)*cos(phi0)*sin(rotationAngle) + sin(theta0)*cos(rotationAngle))
    cosphi = (cos(theta0)*cos(phi0)*cos(rotationAngle) - sin(theta0)*sin(rotationAngle))/cos(theta)
    sinphi = cos(theta0)*sin(phi0)/cos(theta)
    phi = getAngle(cosphi,sinphi)

    ! rotation around z-axis

    !   phi0 = phi
    !   theta0 = theta
    !   rotationAngle = 0.333
    !   theta = theta0
    !   cosphi = (cos(theta0)*cos(phi0)*cos(rotationAngle) - cos(theta0)*sin(phi0)*sin(rotationAngle))/cos(theta)
    !   sinphi = (cos(theta0)*cos(phi0)*sin(rotationAngle) + cos(theta0)*sin(phi0)*cos(rotationAngle))/cos(theta)
    !   phi = getAngle(cosphi,sinphi)

  end subroutine rotateAngles

  function getAngle(cosphi,sinphi) result(phi)

    use definitions

    implicit none
    real(kind=RealKind), intent(in) :: cosphi, sinphi
    real(kind=RealKind) :: phi

    phi = arcsin(sinphi)

    if (cosphi.gt.0.) then
       if (sinphi.gt.0.) then
          phi = phi
       else
          phi = twoPi + phi
       endif
    else
       if (sinphi.gt.0.) then
          phi = pi - phi
       else
          phi = pi - phi
       endif
    endif

  end function getAngle

  recursive subroutine countCells(currentCell)

    use definitions

    implicit none
    type(zoneType), target :: currentCell
    integer :: i, j, k

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call countCells(currentCell%cell(i,j,k))
             enddo
          enddo
       enddo
    else
       icosmic = icosmic + 1
       !     if (currentCell%rho.gt.densestCell%rho) densestCell => currentCell
    endif

  end subroutine countCells

  ! recursive subroutine pixelizeSphere(currentPixel,iAngularLevel,iray)

  !   use definitions

  !   implicit none
  !   type(pixelType), target :: currentPixel
  !   integer, intent(in) :: iAngularLevel, iray
  !   integer :: i, nside
  !   integer*8 :: longIntegerArgument

  !   nside = 2**(iAngularLevel-1)
  !   longIntegerArgument = iray - 1
  ! !  print*, '>>>', longIntegerArgument
  !   call pix2ang_nest(nside,longIntegerArgument,currentPixel%phi,currentPixel%theta)
  !   currentPixel%level = iAngularLevel

  !   if (iAngularLevel.lt.nAngularLevel) then
  !      allocate(currentPixel%pixel(4))
  !      do i = 1, 4
  !         currentPixel%pixel(i)%parent => currentPixel
  !         call pixelizeSphere(currentPixel%pixel(i),iAngularLevel+1,4*iray+i-4)
  !      enddo
  !   endif

  ! end subroutine pixelizeSphere

  subroutine drawSegment(currentCell,startingPoint,currentPixel,level,callSequence,radius,strategy,length,nx,ny,nz)

    use definitions
    use localDefinitions

    implicit none
    type(zoneType) :: currentCell
    type(pointType), intent(inout) :: startingPoint
    type(pointType) :: newStartingPoint, endPoint
    type(pixelType) :: currentPixel
    integer, intent(in) :: level, nx, ny, nz
    integer, intent(in) :: callSequence(3*level+3)
    real(kind=RealKind), intent(inout) :: radius ! always measured in units of the base grid cell size
    integer, intent(inout) :: strategy
    real(kind=RealKind), intent(out) :: length
    integer :: i, j, k, segmentDirection, side, iray
    real(kind=RealKind) :: tmp1, tmp2, tmp3, tmp, prox, proy, proz
    type(segmentType), pointer :: newSegment

    if (currentCell%level.ne.level) then
       print*, 'error in level'
       stop
    endif

    !  print*, 'entered cell:', currentCell%rho, currentCell%level
    !  print*, callSequence
    !  print*, 'via angle', currentPixel%level

    prox = cos(currentPixel%phi)*cos(currentPixel%theta)
    proy = sin(currentPixel%phi)*cos(currentPixel%theta)
    proz = sin(currentPixel%theta)

    if (proz.gt.0.) then  ! xy
       tmp1 = (1.-startingPoint%z)/proz
    else
       tmp1 = - startingPoint%z/proz
    endif

    if (prox.gt.0.) then  ! yz
       tmp2 = (1.-startingPoint%x)/prox
    else
       tmp2 = - startingPoint%x/prox
    endif

    if (proy.gt.0.) then  ! xz
       tmp3 = (1.-startingPoint%y)/proy
    else
       tmp3 = - startingPoint%y/proy
    endif

    !  write(*,*) 'tmp =', tmp1, tmp2, tmp3

    if (tmp1.lt.min(tmp2,tmp3)) then
       segmentDirection = xyPlane
       tmp = tmp1
    else
       if (tmp2.lt.min(tmp1,tmp3)) then
          segmentDirection = yzPlane
          tmp = tmp2
       else
          segmentDirection = xzPlane
          tmp = tmp3
       endif
    endif

    !  print*, 'dir =', segmentDirection, radius*float(2**level)+tmp, rmax(currentPixel%level)
    !  print*, 'from =', startingPoint

    ! create segment from startingPoint to min(tmp,rmax(currentPixel%level))

    !  if (associated(currentCell%lastSegment)) then
    !     allocate(currentCell%lastSegment%segment)
    !     currentCell%lastSegment => currentCell%lastSegment%segment
    !  else
    !     allocate(currentCell%segment)
    !     currentCell%lastSegment => currentCell%segment
    !  endif
    !  newSegment => currentCell%lastSegment

    if ((radius*float(2**level)+tmp.lt.rmax(currentPixel%level)) &
         .or. currentPixel%level.eq.maxPixelLevel) then ! do not split above maxPixelLevel

       !     if (currentPixel%level.eq.maxPixelLevel) print*, radius, level, currentPixel%level

       strategy = proceed
       !     newSegment%len = tmp
       length = tmp
       radius = radius + tmp/float(2**level)
       endPoint%x = startingPoint%x + tmp*prox
       endPoint%y = startingPoint%y + tmp*proy
       endPoint%z = startingPoint%z + tmp*proz
       !     print*, '  to =', endPoint
       select case (segmentDirection)
       case(xyPlane)
          if (proz.lt.0.) then
             side = 0
          else
             side = 1
          endif
          call findXYNeighbour(currentCell,level,callSequence,endPoint%x,endPoint%y,side,strategy,nx,ny,nz)
          if (strategy.ne.boundary) then
             select case (side)
             case(0)
                startingPoint%z = 1.
             case(1)
                startingPoint%z = 0.
             end select
             startingPoint%x = xneighbour
             startingPoint%y = yneighbour
             call checkPoint(startingPoint)
          endif
       case(yzPlane)
          if (prox.lt.0.) then
             side = 0
          else
             side = 1
          endif
          call findYZNeighbour(currentCell,level,callSequence,endPoint%y,endPoint%z,side,strategy,nx,ny,nz)
          if (strategy.ne.boundary) then
             select case (side)
             case(0)
                startingPoint%x = 1.
             case(1)
                startingPoint%x = 0.
             end select
             startingPoint%y = yneighbour
             startingPoint%z = zneighbour
             call checkPoint(startingPoint)
          endif
       case(xzPlane)
          if (proy.lt.0.) then
             side = 0
          else
             side = 1
          endif
          call findXZNeighbour(currentCell,level,callSequence,endPoint%x,endPoint%z,side,strategy,nx,ny,nz)
          if (strategy.ne.boundary) then
             select case (side)
             case(0)
                startingPoint%y = 1.
             case(1)
                startingPoint%y = 0.
             end select
             startingPoint%x = xneighbour
             startingPoint%z = zneighbour
             call checkPoint(startingPoint)
          endif
       end select

    else

       if (radius*float(2**level).ge.rmax(currentPixel%level)) then

          ! refine while entering the cell

          strategy = split
          !        newSegment%len = 0.
          length = 0.

          !        print*, '  to =', startingPoint

       else

          ! refine inside the cell

          strategy = split
          tmp = rmax(currentPixel%level) - radius*float(2**level)
          !        newSegment%len = tmp
          length = tmp
          radius = radius + tmp/float(2**level)

          endPoint%x = startingPoint%x + tmp*prox
          endPoint%y = startingPoint%y + tmp*proy
          endPoint%z = startingPoint%z + tmp*proz

          !        print*, '  to =', endPoint

          startingPoint = endPoint

       endif

    endif

  end subroutine drawSegment

  recursive subroutine localizeCellFromStar(currentCell,currentStar)

    use definitions

    implicit none
    type(zoneType), target :: currentCell
    type(starType) :: currentStar
    integer :: i, j, k, ipos, ir

    if (currentStar%level.gt.currentCell%level) then
       if (.not.currentCell%refined) then
          write(*,*) 'error in star particle position: cell not refined'
          stop
       endif
       ipos = 3*currentCell%level + 4
       i = currentStar%position(ipos)
       j = currentStar%position(ipos+1)
       k = currentStar%position(ipos+2)
       call localizeCellFromStar(currentCell%cell(i,j,k),currentStar)
    else
       currentStar%hostCell => currentCell
    endif

  end subroutine localizeCellFromStar

  recursive subroutine localizeCellFromNode(currentCell,currentNode)

    use definitions

    implicit none
    type(zoneType), target :: currentCell
    type(nodeType) :: currentNode
    integer :: i, j, k, ipos, ir

    if (currentNode%level.gt.currentCell%level) then
       if (.not.currentCell%refined) then
          write(*,*) 'error in node particle position: cell not refined'
          stop
       endif
       ipos = 3*currentCell%level + 4
       i = currentNode%position(ipos)
       j = currentNode%position(ipos+1)
       k = currentNode%position(ipos+2)
       call localizeCellFromNode(currentCell%cell(i,j,k),currentNode)
    else
       currentNode%hostCell => currentCell
    endif

  end subroutine localizeCellFromNode

  recursive subroutine findXYNeighbour(currentCell,level,callSequence,x0,y0,side,strategy,nx,ny,nz)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    integer, intent(in) :: level, side, nx, ny, nz
    integer, intent(in) :: callSequence(3*level+3)
    real(kind=RealKind), intent(in) :: x0, y0
    integer, intent(inout) :: strategy
    integer :: i, j, k
    real(kind=RealKind) :: xnew, ynew
    type(zoneType), pointer :: neighbourContainerCell
    integer :: newCallSequence(3*level+3)

    if (level.gt.0) then
       if ((side.eq.0 .and. callSequence(3*level+3).eq.1) .or. &
            (side.eq.1 .and. callSequence(3*level+3).eq.2)) then
          if (callSequence(3*level+1).eq.1) then
             xnew = 0.5*x0
          else
             xnew = 0.5*x0 + 0.5
          endif
          if (callSequence(3*level+2).eq.1) then
             ynew = 0.5*y0
          else
             ynew = 0.5*y0 + 0.5
          endif
          call findXYNeighbour(currentCell%parent,level-1,callSequence(1:3*level),xnew,ynew,side,strategy,nx,ny,nz)
       else
          newCallSequence = callSequence
          if (side.eq.0) then
             newCallSequence(3*level+3) = 1
          else
             newCallSequence(3*level+3) = 2
          endif
          neighbourContainerCell => currentCell%parent%cell(newCallSequence(3*level+1), &
               newCallSequence(3*level+2),newCallSequence(3*level+3))
          call zoomXYNeighbour(neighbourContainerCell,level,newCallSequence,x0,y0,side)
       endif
    else
       if ((side.eq.0 .and. callSequence(3).eq.1) .or. &
            (side.eq.1 .and. callSequence(3).eq.nz)) then
          !        write(*,*) 'reached boundary'
          strategy = boundary
       else
          newCallSequence = callSequence
          if (side.eq.0) then
             newCallSequence(3) = callSequence(3)-1
          else
             newCallSequence(3) = callSequence(3)+1
          endif
          neighbourContainerCell => currentCell%parent%cell(newCallSequence(1), &
               newCallSequence(2),newCallSequence(3))
          call zoomXYNeighbour(neighbourContainerCell,level,newCallSequence,x0,y0,side)
       endif
    endif

  end subroutine findXYNeighbour

  recursive subroutine findYZNeighbour(currentCell,level,callSequence,y0,z0,side,strategy,nx,ny,nz)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    integer, intent(in) :: level, side, nx, ny, nz
    integer, intent(in) :: callSequence(3*level+3)
    real(kind=RealKind), intent(in) :: y0, z0
    integer, intent(inout) :: strategy
    integer :: i, j, k
    real(kind=RealKind) :: ynew, znew
    type(zoneType), pointer :: neighbourContainerCell
    integer :: newCallSequence(3*level+3)

    if (level.gt.0) then
       if ((side.eq.0 .and. callSequence(3*level+1).eq.1) .or. &
            (side.eq.1 .and. callSequence(3*level+1).eq.2)) then
          if (callSequence(3*level+2).eq.1) then
             ynew = 0.5*y0
          else
             ynew = 0.5*y0 + 0.5
          endif
          if (callSequence(3*level+3).eq.1) then
             znew = 0.5*z0
          else
             znew = 0.5*z0 + 0.5
          endif
          call findYZNeighbour(currentCell%parent,level-1,callSequence(1:3*level),ynew,znew,side,strategy,nx,ny,nz)
       else
          newCallSequence = callSequence
          if (side.eq.0) then
             newCallSequence(3*level+1) = 1
          else
             newCallSequence(3*level+1) = 2
          endif
          neighbourContainerCell => currentCell%parent%cell(newCallSequence(3*level+1), &
               newCallSequence(3*level+2),newCallSequence(3*level+3))
          call zoomYZNeighbour(neighbourContainerCell,level,newCallSequence,y0,z0,side)
       endif
    else
       if ((side.eq.0 .and. callSequence(1).eq.1) .or. &
            (side.eq.1 .and. callSequence(1).eq.nx)) then
          !        write(*,*) 'reached boundary'
          strategy = boundary
       else
          newCallSequence = callSequence
          if (side.eq.0) then
             newCallSequence(1) = callSequence(1)-1
          else
             newCallSequence(1) = callSequence(1)+1
          endif
          neighbourContainerCell => currentCell%parent%cell(newCallSequence(1), &
               newCallSequence(2),newCallSequence(3))
          call zoomYZNeighbour(neighbourContainerCell,level,newCallSequence,y0,z0,side)
       endif
    endif

  end subroutine findYZNeighbour

  recursive subroutine findXZNeighbour(currentCell,level,callSequence,x0,z0,side,strategy,nx,ny,nz)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    integer, intent(in) :: level, side, nx, ny, nz
    integer, intent(in) :: callSequence(3*level+3)
    real(kind=RealKind), intent(in) :: x0, z0
    integer, intent(inout) :: strategy
    integer :: i, j, k
    real(kind=RealKind) :: xnew, znew
    type(zoneType), pointer :: neighbourContainerCell
    integer :: newCallSequence(3*level+3)

    if (level.gt.0) then
       if ((side.eq.0 .and. callSequence(3*level+2).eq.1) .or. &
            (side.eq.1 .and. callSequence(3*level+2).eq.2)) then
          if (callSequence(3*level+1).eq.1) then
             xnew = 0.5*x0
          else
             xnew = 0.5*x0 + 0.5
          endif
          if (callSequence(3*level+3).eq.1) then
             znew = 0.5*z0
          else
             znew = 0.5*z0 + 0.5
          endif
          call findXZNeighbour(currentCell%parent,level-1,callSequence(1:3*level),xnew,znew,side,strategy,nx,ny,nz)
       else
          newCallSequence = callSequence
          if (side.eq.0) then
             newCallSequence(3*level+2) = 1
          else
             newCallSequence(3*level+2) = 2
          endif
          neighbourContainerCell => currentCell%parent%cell(newCallSequence(3*level+1), &
               newCallSequence(3*level+2),newCallSequence(3*level+3))
          call zoomXZNeighbour(neighbourContainerCell,level,newCallSequence,x0,z0,side)
       endif
    else
       if ((side.eq.0 .and. callSequence(2).eq.1) .or. &
            (side.eq.1 .and. callSequence(2).eq.ny)) then
          !        write(*,*) 'reached boundary'
          strategy = boundary
       else
          newCallSequence = callSequence
          if (side.eq.0) then
             newCallSequence(2) = callSequence(2)-1
          else
             newCallSequence(2) = callSequence(2)+1
          endif
          neighbourContainerCell => currentCell%parent%cell(newCallSequence(1), &
               newCallSequence(2),newCallSequence(3))
          call zoomXZNeighbour(neighbourContainerCell,level,newCallSequence,x0,z0,side)
       endif
    endif

  end subroutine findXZNeighbour

  recursive subroutine zoomXYNeighbour(neighbourContainerCell,level,callSequence,x0,y0,side)

    use definitions

    implicit none
    type(zoneType), target :: neighbourContainerCell
    integer, intent(in) :: level, side
    integer, intent(in) :: callSequence(3*level+3)
    real(kind=RealKind), intent(in) :: x0, y0
    integer :: i, j, k
    real(kind=RealKind) :: xnew, ynew
    integer :: newCallSequence(3*level+6)

    if (neighbourContainerCell%refined) then
       if (x0.lt.0.5) then
          xnew = 2.*x0
          i = 1
       else
          xnew = 2.*x0 - 1.
          i = 2
       endif
       if (y0.lt.0.5) then
          ynew = 2.*y0
          j = 1
       else
          ynew = 2.*y0 - 1.
          j = 2
       endif
       if (side.eq.0) then
          k = 2
       else
          k = 1
       endif
       newCallSequence(1:3*level+3) = callSequence
       newCallSequence(3*level+4:3*level+6) = (/i,j,k/)
       call zoomXYNeighbour(neighbourContainerCell%cell(i,j,k),level+1,newCallSequence,xnew,ynew,side)
    else
       neighbourCell => neighbourContainerCell
       xneighbour = x0
       yneighbour = y0
       neighbourCallSequence(1:3*neighbourCell%level+3) = callSequence
    endif

  end subroutine zoomXYNeighbour

  recursive subroutine zoomYZNeighbour(neighbourContainerCell,level,callSequence,y0,z0,side)

    use definitions

    implicit none
    type(zoneType), target :: neighbourContainerCell
    integer, intent(in) :: level, side
    integer, intent(in) :: callSequence(3*level+3)
    real(kind=RealKind), intent(in) :: y0, z0
    integer :: i, j, k
    real(kind=RealKind) :: ynew, znew
    integer :: newCallSequence(3*level+6)

    if (neighbourContainerCell%refined) then
       if (side.eq.0) then
          i = 2
       else
          i = 1
       endif
       if (y0.lt.0.5) then
          ynew = 2.*y0
          j = 1
       else
          ynew = 2.*y0 - 1.
          j = 2
       endif
       if (z0.lt.0.5) then
          znew = 2.*z0
          k = 1
       else
          znew = 2.*z0 - 1.
          k = 2
       endif
       newCallSequence(1:3*level+3) = callSequence
       newCallSequence(3*level+4:3*level+6) = (/i,j,k/)
       call zoomYZNeighbour(neighbourContainerCell%cell(i,j,k),level+1,newCallSequence,ynew,znew,side)
    else
       neighbourCell => neighbourContainerCell
       yneighbour = y0
       zneighbour = z0
       neighbourCallSequence(1:3*neighbourCell%level+3) = callSequence
    endif

  end subroutine zoomYZNeighbour

  recursive subroutine zoomXZNeighbour(neighbourContainerCell,level,callSequence,x0,z0,side)

    use definitions

    implicit none
    type(zoneType), target :: neighbourContainerCell
    integer, intent(in) :: level, side
    integer, intent(in) :: callSequence(3*level+3)
    real(kind=RealKind), intent(in) :: x0, z0
    integer :: i, j, k
    real(kind=RealKind) :: xnew, znew
    integer :: newCallSequence(3*level+6)

    if (neighbourContainerCell%refined) then
       if (x0.lt.0.5) then
          xnew = 2.*x0
          i = 1
       else
          xnew = 2.*x0 - 1.
          i = 2
       endif
       if (side.eq.0) then
          j = 2
       else
          j = 1
       endif
       if (z0.lt.0.5) then
          znew = 2.*z0
          k = 1
       else
          znew = 2.*z0 - 1.
          k = 2
       endif
       newCallSequence(1:3*level+3) = callSequence
       newCallSequence(3*level+4:3*level+6) = (/i,j,k/)
       call zoomXZNeighbour(neighbourContainerCell%cell(i,j,k),level+1,newCallSequence,xnew,znew,side)
    else
       neighbourCell => neighbourContainerCell
       xneighbour = x0
       zneighbour = z0
       neighbourCallSequence(1:3*neighbourCell%level+3) = callSequence
    endif

  end subroutine zoomXZNeighbour

  subroutine checkPoint(point)

    use definitions

    implicit none
    type(pointType) :: point

    if (point%x.lt.0. .or. point%x.gt.1. .or. &
         point%y.lt.0. .or. point%y.gt.1. .or. &
         point%z.lt.0. .or. point%z.gt.1.) then
       write(*,*) 'error in coordinates', point%x, point%y, point%z
       stop
    endif

  end subroutine checkPoint

  subroutine printPixel(pixel,iray)

    use definitions

    implicit none
    type(pixelType) :: pixel
    integer, intent(in) :: iray

    select case (pixel%level)
    case(1)
       write(*,'("pixel =", i5)') iray
    case(2)
       write(*,'("pixel =    ", i5)') iray
    case(3)
       write(*,'("pixel =        ", i5)') iray
    case(4)
       write(*,'("pixel =            ", i5)') iray
    case(5)
       write(*,'("pixel =                ", i5)') iray
    case(6)
       write(*,'("pixel =                    ", i5)') iray
    case(7)
       write(*,'("pixel =                        ", i5)') iray
    case(8)
       write(*,'("pixel =                            ", i5)') iray
    case(9)
       write(*,'("pixel =                                ", i5)') iray
    case(10)
       write(*,'("pixel =                                    ", i5)') iray
    end select

  end subroutine printPixel

  recursive subroutine absoluteCoordinates(level,callSequence,point,nx,ny,nz)

    use definitions

    implicit none
    integer, intent(in) :: level, nx, ny, nz
    integer, intent(in), dimension(3*level+3) :: callSequence
    type(pointType), intent(in) :: point
    type(pointType) :: newPoint

    if (level.gt.0) then
       select case (callSequence(3*level+1))
       case(1)
          newPoint%x = 0.5*point%x
       case(2)
          newPoint%x = 0.5*point%x + 0.5
       end select
       select case (callSequence(3*level+2))
       case(1)
          newPoint%y = 0.5*point%y
       case(2)
          newPoint%y = 0.5*point%y + 0.5
       end select
       select case (callSequence(3*level+3))
       case(1)
          newPoint%z = 0.5*point%z
       case(2)
          newPoint%z = 0.5*point%z + 0.5
       end select
       call absoluteCoordinates(level-1,callSequence(1:3*level),newPoint,nx,ny,nz)
    else
       xbase = (float(callSequence(1)-1)+point%x)/float(nx)
       ybase = (float(callSequence(2)-1)+point%y)/float(ny)
       zbase = (float(callSequence(3)-1)+point%z)/float(nz)
    endif

  end subroutine absoluteCoordinates

  recursive subroutine localizeSplitContinuationCell(x,y,z,level,nx,ny,nz)

    use definitions

    implicit none
    real(kind=RealKind), intent(in) :: x, y, z
    integer, intent(in) :: level, nx, ny, nz
    integer :: i, j, k
    real(kind=RealKind) :: xnew, ynew, znew

    if (level.eq.0) then
       i = int(x*nx) + 1
       j = int(y*ny) + 1
       k = int(z*nz) + 1
       splitContinuationCell => baseGrid
    else
       if (x.lt.0.5) then
          i = 1
       else
          i = 2
       endif
       if (y.lt.0.5) then
          j = 1
       else
          j = 2
       endif
       if (z.lt.0.5) then
          k = 1
       else
          k = 2
       endif
    endif

    splitContinuationCallSequence(3*level+1:3*level+3) = (/i,j,k/)
    splitContinuationCell => splitContinuationCell%cell(i,j,k)

    if (level.eq.0) then
       xnew = x*float(nx) - float(i-1)
       ynew = y*float(ny) - float(j-1)
       znew = z*float(nz) - float(k-1)
    else
       select case (i)
       case(1)
          xnew = 2.*x
       case(2)
          xnew = 2.*x - 1.
       end select
       select case (j)
       case(1)
          ynew = 2.*y
       case(2)
          ynew = 2.*y - 1.
       end select
       select case (k)
       case(1)
          znew = 2.*z
       case(2)
          znew = 2.*z - 1.
       end select
    endif

    if (splitContinuationCell%refined) then
       call localizeSplitContinuationCell(xnew,ynew,znew,level+1,nx,ny,nz)
    else
       splitContinuation%x = xnew
       splitContinuation%y = ynew
       splitContinuation%z = znew
    endif

  end subroutine localizeSplitContinuationCell

  recursive subroutine startNewLongRay(startingCell,startingPoint,currentPixel,irayStarting,level, &
       startingSequence,startingRadius,ndot1,depthStart1,depthStart2,depthStart3,depthStartDust,nx,ny,nz)

    use definitions
    use localDefinitions

    implicit none
    type(zoneType), target :: startingCell
    type(zoneType), pointer :: currentCell
    type(pointType), intent(in) :: startingPoint
    type(pointType) :: currentPoint
    type(pixelType), target :: currentPixel
    type(pixelType), pointer :: leafPixel
    integer*8, intent(in) :: irayStarting
    integer, intent(in) :: level, nx, ny, nz
    integer, intent(in) :: startingSequence(3*level+3)
    integer, dimension(33) :: callSequence ! maximum 10 levels of cell refinement
    real(kind=RealKind), intent(in) :: startingRadius ! always measured in units of the base grid cell size
    real(kind=RealKind) :: radius, oldRadius ! always measured in units of the base grid cell size
    real(kind=RealKind), intent(in) :: ndot1, depthStart1, depthStart2, depthStart3, depthStartDust
    integer :: i, j, k, segmentDirection, side, strategy, count, newLevel, nside, iradius, ienergy
    real(kind=RealKind) :: tmp1, tmp2, tmp3, tmp, prox, proy, proz, &
         tau1, tau2, tau3, tauDust, physicalCellSize, physicalLength, krate1, krate2, krate3, &
         depth1, depth2, depth3, depthDust, etmp1, etmp2, freq, &
         outputDepth1, outputDepth2, outputDepth3, outputDepthDust, &
         outputDepthThreshold1, outputDepthThreshold2, outputDepthThreshold3, outputDepthThresholdDust, &
         ratio, opacityCoefficient
    type(segmentType), pointer :: newSegment
    type(segmentType) :: currentSegment
    integer*8 :: longIntegerArgument, iray

    if (startingCell%level.ne.level) then
       print*, 'error in level'
       stop
    endif

    currentCell => startingCell
    currentPoint = startingPoint
    radius = startingRadius
    depth1 = depthStart1
    depth2 = depthStart2
    depth3 = depthStart3
    depthDust = depthStartDust
    newLevel = currentCell%level
    count = 3*newLevel + 3
    callSequence(1:count) = startingSequence
    strategy = proceed

    do while (strategy.eq.proceed)

       oldRadius = radius
       call drawSegment(currentCell,currentPoint,currentPixel,newLevel, &
            callSequence(1:count),radius,strategy,currentSegment%len,nx,ny,nz)

       !     print*, 'radius =', radius*physicalBoxSize/(float(nx)*kpc)
       !     print*, 'radius =', radius
       physicalCellSize = physicalBoxSize / (float(2**newLevel)*float(nx))
       physicalLength = physicalCellSize * currentSegment%len

       ! optical depths at each reaction's threshold
       tau1 = physicalLength * currentCell%HI * 6.3e-18
       tau2 = physicalLength * currentCell%HeI * 7.42e-18
       tau3 = physicalLength * currentCell%HeII * 1.58e-18

       ! dust optical depth at the Ly-limit (hydrogenIonization)
       select case (dustApproximation)
       case (noDust)
          tauDust = 0.
       case (completeSublimation)
          tauDust = physicalLength * currentCell%HI * 5.4116737e-22 * currentCell%abun2/0.2
          !         if (tauDust.gt.maxDepth) then
          !            print*, '>>>', maxDepth, currentCell%abun2
          !            maxDepth = tauDust
          !         endif
       case (noSublimation)
          tauDust = physicalLength * psi*currentCell%rho/mh * 5.4116737e-22 * currentCell%abun2/0.2
       end select

       do iradius = 1, nradius
          tmp = outputRadius(iradius)*kpc
          tmp1 = oldRadius*physicalBoxSize/float(nx)
          tmp2 = radius*physicalBoxSize/float(nx)
          if (tmp.ge.tmp1 .and. tmp.le.tmp2) then
             ratio = (tmp-tmp1)/(tmp2-tmp1)
             ndotRemaining(iradius) = ndotRemaining(iradius) + &
                  ndot1*exp(-(ratio*(tau1+tauDust)+depth1+depthDust))
             if (iradius.eq.nradius) then
                ! optical depths at each reaction's threshold
                outputDepthThreshold1 = ratio*tau1 + depth1
                outputDepthThreshold2 = ratio*tau2 + depth2
                outputDepthThreshold3 = ratio*tau3 + depth3
                ! dust optical depth at the Ly-limit (hydrogenIonization)
                outputDepthThresholdDust = ratio*tauDust + depthDust
                ndotDust = ndotDust + ndot1*exp(-outputDepthThresholdDust)
                do ienergy = 1, nenergy
                   ! convert to frequency-dependent optical depths
                   outputDepth1 = outputSigma24(ienergy) / 6.30e-18 * outputDepthThreshold1
                   outputDepth2 = outputSigma26(ienergy) / 7.42e-18 * outputDepthThreshold2
                   outputDepth3 = outputSigma25(ienergy) / 1.58e-18 * outputDepthThreshold3
                   outputDepthDust = outputSigmaDust(ienergy) / 5.4116737e-22 * outputDepthThresholdDust
                   ! in units of stellarPopulation(iStar,iSpectrum,coefSpectrum,nu1) [erg/s/Hz]
                   ndotSpectrum(ienergy) = ndotSpectrum(ienergy) + &
                        ndot1*exp(-(outputDepth1+outputDepth2+outputDepth3+outputDepthDust))
                enddo
             endif
          endif
       enddo

       if (strategy.eq.boundary) then
          tmp = radius*physicalBoxSize/(float(nx)*kpc)
          do iradius = 1, nradius
             if (outputRadius(iradius).gt.tmp) ndotBoundary(iradius) = ndotBoundary(iradius) + ndot1
          enddo
       endif
       ! print*, 'tau =', tau1, tau2, tau3, tauDust
       ! print*, 'kappa =', currentCell%kappa1, currentCell%kappa2, currentCell%kappa3
       ! print*, currentCell%HI, currentCell%HeI, currentCell%HeII
       ! print*, physicalLength, currentCell%kappa1
       ! print*, 'xneu =', currentCell%HI * mh / (psi * currentCell%rho)
       ! print*, 'rho =', currentCell%rho

       if (min(depth1+tau1,depth2+tau2,depth3+tau3,depthDust+tauDust).gt.100.) strategy = boundary

       ! reaction [1/s] and heating [erg/s] rates for the entire cell

       ! might need to refine optical depths here

       call getRatesHydrogenHelium(1,depth1,depth2,depth3,depthDust,tmp1,etmp1)
       call getRatesHydrogenHelium(1,depth1+tau1,depth2,depth3,depthDust,tmp2,etmp2)
       currentCell%krate24 = currentCell%krate24 + ndot1 * (tmp1-tmp2) ! [1/s]
       currentCell%crate24 = currentCell%crate24 + ndot1 * (etmp1-etmp2) ! [erg/s]

       call getRatesHydrogenHelium(2,depth1,depth2,depth3,depthDust,tmp1,etmp1)
       call getRatesHydrogenHelium(2,depth1,depth2+tau2,depth3,depthDust,tmp2,etmp2)
       currentCell%krate26 = currentCell%krate26 + ndot1 * (tmp1-tmp2) ! [1/s]
       currentCell%crate26 = currentCell%crate26 + ndot1 * (etmp1-etmp2) ! [erg/s]

       call getRatesHydrogenHelium(3,depth1,depth2,depth3,depthDust,tmp1,etmp1)
       call getRatesHydrogenHelium(3,depth1,depth2,depth3+tau3,depthDust,tmp2,etmp2)
       currentCell%krate25 = currentCell%krate25 + ndot1 * (tmp1-tmp2) ! [1/s]
       currentCell%crate25 = currentCell%crate25 + ndot1 * (etmp1-etmp2) ! [erg/s]

       ! print*, currentCell%krate24, currentCell%krate26, currentCell%krate25
       ! print*, currentCell%crate24, currentCell%crate26, currentCell%crate25

       ! optical depths at each reaction's threshold
       depth1 = depth1 + tau1
       depth2 = depth2 + tau2
       depth3 = depth3 + tau3
       depthDust = depthDust + tauDust

       if (strategy.eq.proceed) then
          currentCell => neighbourCell
          newLevel = currentCell%level
          count = 3*newLevel+3
          callSequence(1:count) = neighbourCallSequence(1:count)
       endif

    enddo

    if (strategy.eq.split) then

       ! print*, 'split in cell:', currentCell%rho, currentCell%level
       ! print*, callSequence(1:count)
       ! print*, currentPoint
       ! print*, 'radius =', radius*float(2**newLevel)

       ! if (currentPixel%level.ge.nAngularLevel) then
       !    print*, 'increase nAngularLevel', nAngularLevel
       !    print*, 'radius =', radius
       !    print*, level
       !    print*, startingSequence
       ! endif

       if (.not.currentPixel%refined) then
          allocate(currentPixel%pixel(4))
          currentPixel%refined = .true.
          currentPixel%pixel(1:4)%level = -99
       endif

       ! allocate(currentPixel%pixel(4))
       ! currentPixel%refined = .true.

       do iray = 1, 4

          leafPixel => currentPixel%pixel(iray)
          if (leafPixel%level.ne.currentPixel%level+1) then
             leafPixel%level = currentPixel%level + 1
             leafPixel%refined = .false.
             leafPixel%parent => currentPixel
             nside = 2**(leafPixel%level-1)
             longIntegerArgument = 4*irayStarting+iray-5
             ! print*, '>>>', longIntegerArgument, irayStarting, iray
             call pix2ang_nest(nside,longIntegerArgument,leafPixel%phi,leafPixel%theta)
          endif

          if (currentPixel%level+1.gt.highestPixelLevel) highestPixelLevel = currentPixel%level + 1

          ! print*, leafPixel%level, iray, leafPixel%phi, leafPixel%theta
          ! if (leafPixel%level.eq.1 .and. iray.eq.4) stop

          call absoluteCoordinates(newLevel,callSequence(1:count),currentPoint,nx,ny,nz)

          ! print*, 'base =', xbase, ybase, zbase, currentPixel%level

          xbase = xbase + radius/float(nx)* &
               (cos(leafPixel%phi)*cos(leafPixel%theta)- &
               cos(currentPixel%phi)*cos(currentPixel%theta))
          ybase = ybase + radius/float(ny)* &
               (sin(leafPixel%phi)*cos(leafPixel%theta)- &
               sin(currentPixel%phi)*cos(currentPixel%theta))
          zbase = zbase + radius/float(nz)* &
               (sin(leafPixel%theta)-sin(currentPixel%theta))

          ! print*, 'base =', xbase, ybase, zbase

          if (xbase.lt.0. .or. xbase.gt.1. .or. ybase.lt.0. .or. ybase.gt.1. .or. &
               zbase.lt.0. .or. zbase.gt.1.) then
             strategy = boundary
             tmp = radius*physicalBoxSize/(float(nx)*kpc)
             do iradius = 1, nradius
                if (outputRadius(iradius).gt.tmp) ndotBoundary(iradius) = ndotBoundary(iradius) + ndot1/4.
             enddo
          endif

          if (strategy.ne.boundary) then

             call localizeSplitContinuationCell(xbase,ybase,zbase,0,nx,ny,nz)
             call checkPoint(splitContinuation)

             ! print*, splitContinuation
             ! stop

             newLevel = splitContinuationCell%level
             count = 3*newLevel + 3

             ! print*, 'continued in cell:', splitContinuationCell%rho, newLevel
             ! print*, splitContinuationCallSequence(1:count)
             ! print*, splitContinuation
             ! write(*,*) '1) r =', radius*float(2**newLevel), newLevel

             ! call printPixel(leafPixel,iray)
             call startNewLongRay(splitContinuationCell,splitContinuation, &
                  leafPixel,4*irayStarting+iray-4, &
                  newLevel,splitContinuationCallSequence(1:count),radius, &
                  ndot1/4.d0,depth1,depth2,depth3,depthDust,nx,ny,nz)

             newLevel = currentCell%level
             count = 3*newLevel + 3

             ! write(*,*) '2) r =', radius*float(2**newLevel), newLevel
             ! print*, 'back in cell:', currentCell%rho, newLevel
             ! print*, callSequence(1:count)
             ! print*, currentPoint
             ! print*, 'radius =', radius*float(2**newLevel)

          endif

       enddo

       ! deallocate(currentPixel%pixel)
       ! currentPixel%refined = .false.

    endif

  end subroutine startNewLongRay

  recursive subroutine analyzeCell(currentCell,nx)

    use definitions

    implicit none
    type(zoneType), target :: currentCell
    integer, intent(in) :: nx
    integer :: i, j, k
    type(segmentType), pointer :: thisSegment

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call analyzeCell(currentCell%cell(i,j,k),nx)
             enddo
          enddo
       enddo
    else
       if (associated(currentCell%segment)) then
          thisSegment => currentCell%segment
          print*, thisSegment%ndot1, thisSegment%len/(float(2**currentCell%level)*float(nx))
          do while (associated(thisSegment%segment))
             thisSegment => thisSegment%segment
             print*, thisSegment%ndot1, thisSegment%len/(float(2**currentCell%level)*float(nx))
          enddo
       endif
    endif

  end subroutine analyzeCell

  recursive subroutine scanCube(currentCell,level,callSequence,nx,ny,nz)

    use definitions

    implicit none
    type(zoneType), target :: currentCell
    integer, intent(in) :: level, nx, ny, nz
    integer, intent(in) :: callSequence(3*level+3)
    integer :: i, j, k, callLength
    integer :: newCallSequence(3*level+6)
    type(pointType) :: point

    callLength = 3*level + 3

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                newCallSequence(1:callLength) = callSequence
                newCallSequence(callLength+1:callLength+3) = (/i,j,k/)
                call scanCube(currentCell%cell(i,j,k),level+1,newCallSequence,nx,ny,nz)
             enddo
          enddo
       enddo
    else
       if (.not.associated(currentCell%segment)) then
          print*, 'no rays inside'
          print*, level
          print*, callSequence
          stop
       endif
       point%x = 0.5
       point%y = 0.5
       point%z = 0.5
       call absoluteCoordinates(level,callSequence,point,nx,ny,nz)
       print*, sqrt((xbase-xsrc)**2+(ybase-ysrc)**2+(zbase-zsrc)**2), &
            currentCell%krate24 * float(2**currentCell%level)**3
    endif

  end subroutine scanCube

  recursive subroutine solveRateEquations(currentCell,nx,runUVBTransfer)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    integer :: i, j, k, indixe, iRedshift, iDensity, iSpecificEnergy
    integer, intent(in) :: nx
    logical, intent(in) :: runUVBTransfer
    real(kind=RealKind) :: HI, HII, HeI, HeII, HeIII, tgas, de, de1, de2, res, res1, res2, HeIprev
    real(kind=RealKind) :: de10, de20, res10, res20
    real(kind=RealKind) :: HIp, HIIp, HeIp, HeIIp, HeIIIp, tgasp, dep
    real(kind=RealKind) :: k1, k2, k3, k4, k5, k6, &
         logtem, t1, t2, tdef, nden, en
    real(kind=RealKind) :: ceHI, ceHeI, ceHeII, &
         ciHI, ciHeI, ciHeIS, ciHeII, reHII, reHeII1, reHeII2, reHeIII, brem, lineHI, &
         krate24, krate25, krate26
    real(kind=RealKind) :: scoef, acoef, ttmin, nh, nhe, x1, x2, x3, delta, &
         comp1, comp2, tmp, tmp1, tmp2, tmp3, physicalCellSize, &
         meanFreePathLymanLimit

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call solveRateEquations(currentCell%cell(i,j,k),nx,runUVBTransfer)
             enddo
          enddo
       enddo
    else

       icosmic = icosmic + 1

       if (100000*(icosmic/100000).eq.icosmic) then
          write(*,fmt='(f7.2,"% done")') real(icosmic)/real(ncosmic)*100.
       endif

       nh = psi * currentCell%rho / mh
       nhe = (1.-psi) * currentCell%rho / mhe

       HI = min(currentCell%HI,nh)
       HII = nh - currentCell%HI
       HeI = currentCell%HeI
       HeII = currentCell%HeII
       HeIII = nhe - currentCell%HeI - currentCell%HeII

       if (HeIII.lt.0.) then
          currentCell%HeII = nhe - currentCell%HeI
          HeIII = 0.
          if (HeII.lt.0.) then
             currentCell%HeI = nhe
             currentCell%HeII = 0.
             HeII = 0.
             HeIII = 0.
          endif
       endif

       tgas = currentCell%tgas
       de = HII + HeII + 2.*HeIII ! electron # density in [1/cm^3]

       physicalCellSize = physicalBoxSize / (float(2**currentCell%level)*float(nx))

       ! convert rates from [1/s] per cell to [1/s] per absorber particle

       if (HI.gt.0.d0) then
          krate24 = currentCell%krate24 / (physicalCellSize**3 * HI)
       else
          krate24 = 0.
       endif

       if (HeII.gt.0.d0) then
          krate25 = currentCell%krate25 / (physicalCellSize**3 * HeII)
       else
          krate25 = 0.
       endif

       if (HeI.gt.0.d0) then
          krate26 = currentCell%krate26 / (physicalCellSize**3 * HeI)
       else
          krate26 = 0.
       endif

       krate24 = max(krate24, 0.d0)
       krate25 = max(krate25, 0.d0)
       krate26 = max(krate26, 0.d0)


       if (runUVBTransfer) then
          ! add non-uniform UVB
          tmp1 = 4.*pi*currentCell%Jmean1 ! [erg/(cm^2 s Hz)]
          tmp2 = 4.*pi*currentCell%Jmean2
          tmp3 = 4.*pi*currentCell%Jmean3
          krate24 = krate24 + tmp1 * group1%ksi24 + tmp2 * group2%ksi24 + tmp3 * group3%ksi24 ! [1/s]
          krate25 = krate25 + tmp3 * group3%ksi25
          krate26 = krate26 + tmp2 * group2%ksi26 + tmp3 * group3%ksi26
       else
          ! add uniform UVB
          meanFreePathLymanLimit = 1./(HI*6.3e-18 + HeI*7.42e-18 + HeII*1.58e-18)
          if (meanFreePathLymanLimit.ge.selfShieldingThreshold) then
             krate24 = krate24 + 4.*pi * (uniformQuasar*quasar%ksi24 + uniformStellar*stellar%ksi24)
             krate25 = krate25 + 4.*pi * (uniformQuasar*quasar%ksi25 + uniformStellar*stellar%ksi25)
             krate26 = krate26 + 4.*pi * (uniformQuasar*quasar%ksi26 + uniformStellar*stellar%ksi26)
          endif
       endif

       de = HII + HeII + 2.*HeIII ! electron # density in [1/cm^3]

       ! Compute temp-centered temperature (and log)

       logtem = log(tgas)
       logtem = max(logtem, logtem0)
       logtem = min(logtem, logtem9)

       ! Find index into table and precompute interpolation values

       indixe = min(nratec-1,max(1,int((logtem-logtem0)/dlogtem)+1))
       t1 = (logtem0 + (indixe - 1)*dlogtem)
       t2 = (logtem0 + (indixe    )*dlogtem)
       tdef = t2 - t1

       ! Do linear table lookup (in log temperature)

       k1 = k1a(indixe) + (logtem - t1)*(k1a(indixe+1)-k1a(indixe))/tdef
       k2 = k2a(indixe) + (logtem - t1)*(k2a(indixe+1)-k2a(indixe))/tdef
       k3 = k3a(indixe) + (logtem - t1)*(k3a(indixe+1)-k3a(indixe))/tdef
       k4 = k4a(indixe) + (logtem - t1)*(k4a(indixe+1)-k4a(indixe))/tdef
       k5 = k5a(indixe) + (logtem - t1)*(k5a(indixe+1)-k5a(indixe))/tdef
       k6 = k6a(indixe) + (logtem - t1)*(k6a(indixe+1)-k6a(indixe))/tdef

       ! ionization equilibrium

       de1 = 1.e-30
       de = de1
       HeI = (de - nh / (1. + k2*de / (k1*de + krate24)) - 2.*nhe) &
            / ((k3*de + krate26)/(k4*de) - 2. - 2.*(k3*de + krate26)/(k4*de))
       res1 = k3*HeI*de + k6*(nhe - HeI - HeI*(k3*de + krate26)/(k4*de))*de + &
            krate26*HeI - HeI*(k3*de + krate26)/(k4*de) * (k4*de + k5*de + krate25)

       de2 = nh + 2.*nhe
       de = de2
       HeI = (de - nh / (1. + k2*de / (k1*de + krate24)) - 2.*nhe) &
            / ((k3*de + krate26)/(k4*de) - 2. - 2.*(k3*de + krate26)/(k4*de))
       res2 = k3*HeI*de + k6*(nhe - HeI - HeI*(k3*de + krate26)/(k4*de))*de + &
            krate26*HeI - HeI*(k3*de + krate26)/(k4*de) * (k4*de + k5*de + krate25)

       HeIprev = - 1.

       de10 = de1
       de20 = de2
       res10 = res1
       res20 = res2

       do while (abs(HeI-HeIprev)/nhe.gt.1.d-10)
          HeIprev = HeI
          !           de = (de2*abs(res1)+de1*abs(res2))/(abs(res2)+abs(res1))
          de = 0.5 * (de1+de2)
          HeI = (de - nh / (1. + k2*de / (k1*de + krate24)) - 2.*nhe) &
               / ((k3*de + krate26)/(k4*de) - 2. - 2.*(k3*de + krate26)/(k4*de))
          res = k3*HeI*de + k6*(nhe - HeI - HeI*(k3*de + krate26)/(k4*de))*de + &
               krate26*HeI - HeI*(k3*de + krate26)/(k4*de) * (k4*de + k5*de + krate25)
          if (opposite(res,res1)) then
             de2 = de
             res2 = res
          else
             de1 = de
             res1 = res
          endif
          !           print*, de, res, HeI/nhe
       enddo

       HeII = HeI*(k3*de + krate26)/(k4*de)
       HeIII = nhe - HeI - HeII
       HII  = nh / (1. + k2*de / (k1*de + krate24))
       HI = k2*HII*de / (k1*de + krate24)

       if (HI/nh.ge.0. .and. HI/nh.le.1.) then
          continue
       else
          print*, 'H:', HI, HII
          print*, 'He:', HeI, HeII, HeIII
          stop
       endif
       if (HeI/nhe.ge.0. .and. HeI/nhe.le.1.) then
          continue
       else
          print*, 'H:', HI, HII
          print*, 'He:', HeI, HeII, HeIII
          print*, de10, de20
          print*, res10, res20
          print*, '----------------'
          print*, nhe, currentCell%HeI, currentCell%HeII
          stop
       endif

       ! balance check
       ! tmp = (HI+HII)/(HIp+HIIp)
       ! HI = HIp * tmp
       ! HII = HIIp * tmp
       ! tmp = (HeI+HeII+HeIII)/(HeIp+HeIIp+HeIIIp)
       ! HeI = HeIp * tmp
       ! HeII = HeIIp * tmp
       ! HeIII = HeIIIp * tmp
       ! nhe = (1.-psi) * currentCell%rho / mhe
       ! HeIII = (1.-psi) * currentCell%rho / mhe - currentCell%HeI - currentCell%HeII
       ! xneu = currentCell%HI * mh / (psi * currentCell%rho)
       ! HII = psi * currentCell%rho / mh - currentCell%HI

       tmp1 = abs(HI-currentCell%HI) * mh / (psi * currentCell%rho)
       tmp2 = abs(HeI-currentCell%HeI) * mhe / ((1.-psi) * currentCell%rho)
       tmp3 = abs(HeII-currentCell%HeII) * mhe / ((1.-psi) * currentCell%rho)
       tmp = max(tmp1,tmp2,tmp3)

       currentCell%HI = HI
       currentCell%HeI = HeI
       currentCell%HeII = HeII

    endif

  end subroutine solveRateEquations

  recursive subroutine initialIonizationEquilibrium(currentCell)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    integer :: i, j, k, indixe, iRedshift, iDensity, iSpecificEnergy
    real(kind=RealKind) :: HI, HII, HeI, HeII, HeIII, tgas, de, de1, de2, res, res1, res2, &
         HIprev, HeIprev, HeIIprev

    real(kind=RealKind) :: k1, k2, k3, k4, k5, k6, &
         logtem, t1, t2, tdef, nden, en, &
         krate24, krate25, krate26, &
         nh, nhe, tmp, tmp1, tmp2, tmp3, meanFreePathLymanLimit

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call initialIonizationEquilibrium(currentCell%cell(i,j,k))
             enddo
          enddo
       enddo
    else

       icosmic = icosmic + 1

       if (100000*(icosmic/100000).eq.icosmic) then
          write(*,fmt='(f7.2,"% done")') real(icosmic)/real(ncosmic)*100.
       endif

       HI = min(currentCell%HI,psi*currentCell%rho/mh)
       HII = psi * currentCell%rho / mh - currentCell%HI
       HeI = currentCell%HeI
       HeII = currentCell%HeII
       HeIII = (1.-psi) * currentCell%rho / mhe - currentCell%HeI - currentCell%HeII

       ! print*, HI/(HI+HII), HeI/(HeI+HeII+HeIII), HeII/(HeI+HeII+HeIII)

       if (HeIII.lt.0.) then
          currentCell%HeII = (1.-psi) * currentCell%rho / mhe - currentCell%HeI
          HeIII = 0.
          if (currentCell%HeII.lt.0.) then
             currentCell%HeI = (1.-psi) * currentCell%rho / mhe
             currentCell%HeII = 0.
             HeII = 0.
             HeIII = 0.
          endif
       endif

       tgas = currentCell%tgas
       de = HII + HeII + 2.*HeIII ! electron # density in [1/cm^3]

       ! add uniform UVB

       meanFreePathLymanLimit = 1./(HI*6.3e-18 + HeI*7.42e-18 + HeII*1.58e-18)
       if (meanFreePathLymanLimit.ge.selfShieldingThreshold) then
          krate24 = 4.*pi * (uniformQuasar*quasar%ksi24 + uniformStellar*stellar%ksi24)
          krate25 = 4.*pi * (uniformQuasar*quasar%ksi25 + uniformStellar*stellar%ksi25)
          krate26 = 4.*pi * (uniformQuasar*quasar%ksi26 + uniformStellar*stellar%ksi26)
       else
          krate24 = 0.
          krate25 = 0.
          krate26 = 0.
       endif

       de = HII + HeII + 2.*HeIII ! electron # density in [1/cm^3]

       ! Compute temp-centered temperature (and log)

       logtem = log(tgas)
       logtem = max(logtem, logtem0)
       logtem = min(logtem, logtem9)

       ! Find index into table and precompute interpolation values

       indixe = min(nratec-1,max(1,int((logtem-logtem0)/dlogtem)+1))
       t1 = (logtem0 + (indixe - 1)*dlogtem)
       t2 = (logtem0 + (indixe    )*dlogtem)
       tdef = t2 - t1

       ! Do linear table lookup (in log temperature)

       k1 = k1a(indixe) + (logtem - t1)*(k1a(indixe+1)-k1a(indixe))/tdef
       k2 = k2a(indixe) + (logtem - t1)*(k2a(indixe+1)-k2a(indixe))/tdef
       k3 = k3a(indixe) + (logtem - t1)*(k3a(indixe+1)-k3a(indixe))/tdef
       k4 = k4a(indixe) + (logtem - t1)*(k4a(indixe+1)-k4a(indixe))/tdef
       k5 = k5a(indixe) + (logtem - t1)*(k5a(indixe+1)-k5a(indixe))/tdef
       k6 = k6a(indixe) + (logtem - t1)*(k6a(indixe+1)-k6a(indixe))/tdef

       nh = psi * currentCell%rho / mh
       nhe = (1.-psi) * currentCell%rho / mhe

       de1 = 1.e-20
       de = de1
       HeI = (de - nh / (1. + k2*de / (k1*de + krate24)) - 2.*nhe) &
            / ((k3*de + krate26)/(k4*de) - 2. - 2.*(k3*de + krate26)/(k4*de))
       res1 = k3*HeI*de + k6*(nhe - HeI - HeI*(k3*de + krate26)/(k4*de))*de + &
            krate26*HeI - HeI*(k3*de + krate26)/(k4*de) * (k4*de + k5*de + krate25)

       de2 = nh + 2.*nhe
       de = de2
       HeI = (de - nh / (1. + k2*de / (k1*de + krate24)) - 2.*nhe) &
            / ((k3*de + krate26)/(k4*de) - 2. - 2.*(k3*de + krate26)/(k4*de))
       res2 = k3*HeI*de + k6*(nhe - HeI - HeI*(k3*de + krate26)/(k4*de))*de + &
            krate26*HeI - HeI*(k3*de + krate26)/(k4*de) * (k4*de + k5*de + krate25)

       HIprev = - 1.
       HeIprev = - 1.
       HeIIprev = - 1.

       ! do while (abs(HeI-HeIprev)/nhe.gt.1.d-10)
       do while (HeI.ne.HeIprev)
          HIprev = HI
          HeIprev = HeI
          HeIIprev = HeII
          ! de = (de2*abs(res1)+de1*abs(res2))/(abs(res2)+abs(res1))
          de = 0.5 * (de1+de2)
          HeI = (de - nh / (1. + k2*de / (k1*de + krate24)) - 2.*nhe) &
               / ((k3*de + krate26)/(k4*de) - 2. - 2.*(k3*de + krate26)/(k4*de))
          res = k3*HeI*de + k6*(nhe - HeI - HeI*(k3*de + krate26)/(k4*de))*de + &
               krate26*HeI - HeI*(k3*de + krate26)/(k4*de) * (k4*de + k5*de + krate25)
          if (opposite(res,res1)) then
             de2 = de
          else
             de1 = de
          endif
          ! print*, de, res
          ! print*, real(HI/nh), real(HeI/nhe)
       enddo
       if (HIprev.ne.HI) then
          print*, 'H not converged to equilibrium ...'
          print*, HIprev, HI
          stop
       endif
       if (HeIIprev.ne.HeII) then
          print*, 'He not converged to equilibrium ...'
          print*, HeIIprev, HeII
          stop
       endif

       HeII = HeI*(k3*de + krate26)/(k4*de)
       HeIII = nhe - HeI - HeII
       HII  = nh / (1. + k2*de / (k1*de + krate24))
       HI = k2*HII*de / (k1*de + krate24)

       ! de = HII + HeII + 2.*HeIII
       ! HI = k2*HII*de / (k1*de + krate24)
       ! HII = psi * currentCell%rho / mh - HI
       ! HeI = k4*HeII*de/(k3*de + krate26)
       ! HeII = (k3*HeI*de + k6*HeIII*de + krate26*HeI)/(k4*de + k5*de + krate25)
       ! HeIII = (1.-psi) * currentCell%rho / mhe - HeI - HeII

       if (HI/nh.ge.0. .and. HI/nh.le.1.) then
          continue
       else
          print*, '>>>', HI
          stop
       endif
       if (HeI/nhe.ge.0. .and. HeI/nhe.le.1.) then
          continue
       else
          print*, '>>>', HeI, HeII
          stop
       endif
       ! else
       !    HI = k2*nh/(k1+k2)
       !    HII = nh - HI
       !    HeII = k6*nhe/(k5 + k6*k4/k3 + k6)
       !    HeI = k4*HeII/k3
       !    HeIII = nhe - k4*HeII/k3 - HeII
       !    de = HII + HeII + 2.*HeIII
       ! endif

       ! HI = k2*(nh-HI)*((nh-HI)+HeII+2.*(nhe-HeI-HeII)) / (k1*((nh-HI)+HeII+2.*(nhe-HeI-HeII)) + krate24)
       ! HeI = k4*HeII*((nh-HI)+HeII+2.*(nhe-HeI-HeII))/(k3*((nh-HI)+HeII+2.*(nhe-HeI-HeII)) + krate26)
       ! HeII = (k3*HeI*((nh-HI)+HeII+2.*(nhe-HeI-HeII)) + k6*(nhe-HeI-HeII)*((nh-HI)+HeII+2.*(nhe-HeI-HeII)) + krate26*HeI)/((k4+k5)*((nh-HI)+HeII+2.*(nhe-HeI-HeII)) + krate25)

       ! HI = k2*(nh-HI)*((nh-HI)+HeII+2.*(nhe-HeI-HeII)) / (k1*((nh-HI)+HeII+2.*(nhe-HeI-HeII)))
       ! HeI = k4*HeII*((nh-HI)+HeII+2.*(nhe-HeI-HeII)) / (k3*((nh-HI)+HeII+2.*(nhe-HeI-HeII)))
       ! HeII = (k3*HeI*((nh-HI)+HeII+2.*(nhe-HeI-HeII)) + k6*(nhe-HeI-HeII)*((nh-HI)+HeII+2.*(nhe-HeI-HeII)))/((k4+k5) * ((nh-HI)+HeII+2.*(nhe-HeI-HeII))

       currentCell%HI = HI
       currentCell%HeI = HeI
       currentCell%HeII = HeII
       currentCell%tgas = tgas

    endif

  end subroutine initialIonizationEquilibrium

  recursive subroutine thermalEquilibrium(currentCell)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    integer :: i, j, k, indixe, iRedshift, iDensity, iSpecificEnergy
    real(kind=RealKind) :: HI, HII, HeI, HeII, HeIII, tgas, de
    real(kind=RealKind) :: HIp, HIIp, HeIp, HeIIp, HeIIIp, tgasp, dep

    real(kind=RealKind) :: k1, k2, k3, k4, k5, k6, &
         dtit, edot, logtem, t1, t2, tdef, nden, en

    real(kind=RealKind) :: ceHI, ceHeI, ceHeII, &
         ciHI, ciHeI, ciHeIS, ciHeII, reHII, reHeII1, reHeII2, reHeIII, brem, lineHI, &
         crate24, crate25, crate26, crate

    real(kind=RealKind) :: scoef, acoef, ttmin, dt, nh, nhe, x1, x2, x3, delta, &
         comp1, comp2, tmp, tmp1, tmp2, tmp3, physicalCellSize, &
         meanFreePathLymanLimit

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call thermalEquilibrium(currentCell%cell(i,j,k))
             enddo
          enddo
       enddo
    else

       icosmic = icosmic + 1

       if (100000*(icosmic/100000).eq.icosmic) then
          write(*,fmt='(f7.2,"% done")') real(icosmic)/real(ncosmic)*100.
       endif

       HI = min(currentCell%HI,psi*currentCell%rho/mh)
       HII = psi * currentCell%rho / mh - currentCell%HI
       HeI = currentCell%HeI
       HeII = currentCell%HeII
       HeIII = (1.-psi) * currentCell%rho / mhe - currentCell%HeI - currentCell%HeII

       if (HeIII.lt.0.) then
          currentCell%HeII = (1.-psi) * currentCell%rho / mhe - currentCell%HeI
          HeIII = 0.
          if (currentCell%HeII.lt.0.) then
             currentCell%HeI = (1.-psi) * currentCell%rho / mhe
             currentCell%HeII = 0.
             HeII = 0.
             HeIII = 0.
          endif
       endif

       tgas = currentCell%tgas
       de = HII + HeII + 2.*HeIII ! electron # density in [1/cm^3]

       ! add uniform UVB

       meanFreePathLymanLimit = 1./(HI*6.3e-18 + HeI*7.42e-18 + HeII*1.58e-18)
       if (meanFreePathLymanLimit.ge.selfShieldingThreshold) then
          crate24 = 4.*pi * (uniformQuasar*quasar%gammaHI + uniformStellar*stellar%gammaHI)
          crate25 = 4.*pi * (uniformQuasar*quasar%gammaHeII + uniformStellar*stellar%gammaHeII)
          crate26 = 4.*pi * (uniformQuasar*quasar%gammaHeI + uniformStellar*stellar%gammaHeI)
       else
          crate24 = 0.
          crate25 = 0.
          crate26 = 0.
       endif

       crate = crate24*HI + crate25*HeII + crate26*HeI ! volumetric heating rate [erg/(s cm^3)]

       de = HII + HeII + 2.*HeIII ! electron # density in [1/cm^3]

       ! Compute temp-centered temperature (and log)

       logtem = log(tgas)
       logtem = max(logtem, logtem0)
       logtem = min(logtem, logtem9)

       ! Find index into table and precompute interpolation values

       indixe = min(nratec-1,max(1,int((logtem-logtem0)/dlogtem)+1))
       t1 = (logtem0 + (indixe - 1)*dlogtem)
       t2 = (logtem0 + (indixe    )*dlogtem)
       tdef = t2 - t1

       ! Do linear table lookup (in log temperature)

       k1 = k1a(indixe) + (logtem - t1)*(k1a(indixe+1)-k1a(indixe))/tdef
       k2 = k2a(indixe) + (logtem - t1)*(k2a(indixe+1)-k2a(indixe))/tdef
       k3 = k3a(indixe) + (logtem - t1)*(k3a(indixe+1)-k3a(indixe))/tdef
       k4 = k4a(indixe) + (logtem - t1)*(k4a(indixe+1)-k4a(indixe))/tdef
       k5 = k5a(indixe) + (logtem - t1)*(k5a(indixe+1)-k5a(indixe))/tdef
       k6 = k6a(indixe) + (logtem - t1)*(k6a(indixe+1)-k6a(indixe))/tdef

       !---------------------- compute cooling rates

       ! Set compton cooling coefficients (and temperature)

       comp1 = compa * (1. + currentRedshift)**4
       comp2 = 2.73  * (1. + currentRedshift)

       ! Lookup cooling values and do a linear temperature in log(T)

       ceHI = ceHIa(indixe) + (logtem-t1) * (ceHIa(indixe+1)-ceHIa(indixe)) / tdef
       ceHeI = ceHeIa(indixe) + (logtem-t1) * (ceHeIa(indixe+1)-ceHeIa(indixe)) / tdef
       ceHeII = ceHeIIa(indixe) + (logtem-t1) * (ceHeIIa(indixe+1)-ceHeIIa(indixe)) / tdef
       ciHI = ciHIa(indixe) + (logtem-t1) * (ciHIa(indixe+1)-ciHIa(indixe)) / tdef
       ciHeI = ciHeIa(indixe) + (logtem-t1) * (ciHeIa(indixe+1)-ciHeIa(indixe)) / tdef
       ciHeIS = ciHeISa(indixe) + (logtem-t1) * (ciHeISa(indixe+1)-ciHeISa(indixe)) / tdef
       ciHeII = ciHeIIa(indixe) + (logtem-t1) * (ciHeIIa(indixe+1)-ciHeIIa(indixe)) / tdef
       reHII = reHIIa(indixe) + (logtem-t1) * (reHIIa(indixe+1)-reHIIa(indixe)) / tdef
       reHeII1=reHeII1a(indixe) + (logtem-t1) * (reHeII1a(indixe+1)-reHeII1a(indixe)) / tdef
       reHeII2=reHeII2a(indixe) + (logtem-t1) * (reHeII2a(indixe+1)-reHeII2a(indixe)) / tdef
       reHeIII=reHeIIIa(indixe) + (logtem-t1) * (reHeIIIa(indixe+1)-reHeIIIa(indixe)) / tdef
       brem = brema(indixe) + (logtem-t1) * (brema(indixe+1)-brema(indixe)) / tdef
       lineHI = lineHIa(indixe) + (logtem-t1) * (lineHIa(indixe+1)-lineHIa(indixe)) / tdef

       ! Compute the cooling function

       edot = ( &           ! [erg/cm^3/s]
            !
            !                    Collisional excitations
            !
            - ceHI  *HI  *de &                ! ce of HI
            - ceHeI *HeI*de**2 &              ! ce of HeI
            - ceHeII*HeII*de &                ! ce of HeII
            !
            !                    Collisional ionizations
            !
            - ciHI  *HI  *de &        ! ci of HI
            - ciHeI *HeI *de &        ! ci of HeI
            - ciHeII*HeII*de &        ! ci of HeII
            - ciHeIS*HeII*de**2 & ! ci of HeIS
            !
            !                    Recombinations
            !
            - reHII  *HII  *de &     ! re of HII
            - reHeII1*HeII *de &     ! re of HeII
            - reHeII2*HeII *de &     ! re of HeII
            - reHeIII*HeIII*de &     ! re of HeIII
            !
            !                    Compton cooling or heating
            !
            - comp1*(tgas-comp2)*de &
            !
            !                    X-ray compton heating
            !
            - comp_xraya * (tgas-comp_temp)*de &
            !
            !                    Bremsstrahlung
            !
            - brem*(HII+HeII+4.*HeIII)*de) &
            !
            !                    HI line excitation cooling
            !
            - lineHI * HI * de * 0.

       edot = edot + crate ! [erg / (s cm^3)]
       currentCell%hydroHeating = - edot

       ! hydroHeating should be positive if there is local
       ! dissipation/feedback, and be almost zero (on the positive side) with
       ! little extra heating (cooling balances heating by the UVB). If it's
       ! negative then UVB heating dominates cooling, and there is no need to
       ! invoke local heating terms -- just set it to zero.

       if (currentCell%hydroHeating.lt.0.) currentCell%hydroHeating = 0.

    endif

  end subroutine thermalEquilibrium

  recursive subroutine writeCell(currentCell,level)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    integer, intent(in) :: level
    integer :: i, j, k

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call writeCell(currentCell%cell(i,j,k),level+1)
             enddo
          enddo
       enddo
    else
       icosmic = icosmic + 1
       cellArrayLevel(icosmic) = level
       cellArrayHI(icosmic) = sngl(currentCell%HI)
       cellArrayHeI(icosmic) = sngl(currentCell%HeI)
       cellArrayHeII(icosmic) = sngl(currentCell%HeII)
       cellArrayTemp(icosmic) = sngl(currentCell%tgas)
       cellArrayDensity(icosmic) = sngl(currentCell%rho)
       if (readKinematics) then
          cellArrayVelx(icosmic) = sngl(currentCell%velx)
          cellArrayVely(icosmic) = sngl(currentCell%vely)
          cellArrayVelz(icosmic) = sngl(currentCell%velz)
       endif
       if (readMetals) then
          cellArrayAbun2(icosmic) = sngl(currentCell%abun2)
       endif
    endif

  end subroutine writeCell

  recursive subroutine cleanUpSegments(currentCell)

    use definitions

    implicit none
    type(zoneType), target :: currentCell
    integer :: i, j, k

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call cleanUpSegments(currentCell%cell(i,j,k))
             enddo
          enddo
       enddo
    else
       if (associated(currentCell%lastSegment)) then
          nullify(currentCell%lastSegment)
          call deleteSegments(currentCell)
          stop
       endif
    endif

    ! if (associated(currentCell%lastSegment)) then
    !    allocate(currentCell%lastSegment%segment)
    !    currentCell%lastSegment => currentCell%lastSegment%segment
    ! else
    !    allocate(currentCell%segment)
    !    currentCell%lastSegment => currentCell%segment
    ! endif
    ! newSegment => currentCell%lastSegment

  end subroutine cleanUpSegments

  recursive subroutine deleteSegments(currentCell)

    use definitions

    implicit none
    type(zoneType), target :: currentCell

    print*, 'deleteSegments not implemented yet'
    stop

  end subroutine deleteSegments

  recursive subroutine setZeroRates(currentCell)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    integer :: i, j, k

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call setZeroRates(currentCell%cell(i,j,k))
             enddo
          enddo
       enddo
    else

       currentCell%krate24 = 0.d0
       currentCell%krate25 = 0.d0
       currentCell%krate26 = 0.d0
       currentCell%crate24 = 0.d0
       currentCell%crate25 = 0.d0
       currentCell%crate26 = 0.d0

    endif

  end subroutine setZeroRates

  subroutine getRatesHydrogenHelium(reaction,tau1,tau2,tau3,tauDust,numberRate,heatingRate)

    use definitions

    implicit none
    integer, intent(in) :: reaction
    real(kind=RealKind), intent(in) :: tau1, tau2, tau3, tauDust
    real(kind=RealKind), intent(out) :: numberRate, heatingRate
    integer :: idepth1, idepth2, idepth3, idepthDust
    real(kind=RealKind) :: c1, c2, c3, cDust, nr1, nr2, hr1, hr2

    if (tau1.gt.maxOpticalDepth1 .or. tau2.gt.maxOpticalDepth2 .or. &
         tau3.gt.maxOpticalDepth3 .or. tauDust.gt.maxOpticalDepthDust) then
       numberRate = 0.
       heatingRate = 0.
       return
    endif

    idepth1 = int(tau1/maxOpticalDepth1*float(ndepth1))
    idepth2 = int(tau2/maxOpticalDepth2*float(ndepth2))
    idepth3 = int(tau3/maxOpticalDepth3*float(ndepth3))

    c1 = tau1*float(ndepth1)/maxOpticalDepth1 - float(idepth1)
    c2 = tau2*float(ndepth2)/maxOpticalDepth2 - float(idepth2)
    c3 = tau3*float(ndepth3)/maxOpticalDepth3 - float(idepth3)

    if (dustApproximation.eq.noDust) then
       idepthDust = 0
       cDust = 0.
    else
       idepthDust = int(tauDust/maxOpticalDepthDust*float(ndepthDust))
       cDust = tauDust*float(ndepthDust)/maxOpticalDepthDust - float(idepthDust)
    endif

    ! if (tau1+tau2+tau3+tauDust.gt.0.) then
    !    print*, tau1, tau2, tau3, tauDust
    !    print*, idepth1,idepth2,idepth3,idepthDust
    !    stop
    ! endif

    if (min(idepth1,idepth2,idepth3,idepthDust).lt.0) then
       print*, 'error in idepth123'
       print*, tau1, tau2, tau3, idepthDust
       stop
    endif

    select case (reaction)
    case(1)
       nr1 = c1 * ((1.-c3)*(1.-c2)*log(reactionRate1(idepth1+1,idepth2,idepth3,idepthDust))+ &
            c3*(1.-c2)*log(reactionRate1(idepth1+1,idepth2,idepth3+1,idepthDust)) + &
            c2*(1.-c3)*log(reactionRate1(idepth1+1,idepth2+1,idepth3,idepthDust)) + &
            c3*c2*log(reactionRate1(idepth1+1,idepth2+1,idepth3+1,idepthDust))) &
            + (1.-c1) * ((1.-c3)*(1.-c2)*log(reactionRate1(idepth1,idepth2,idepth3,idepthDust)) + &
            c3*(1.-c2)*log(reactionRate1(idepth1,idepth2,idepth3+1,idepthDust)) + &
            c2*(1.-c3)*log(reactionRate1(idepth1,idepth2+1,idepth3,idepthDust)) + &
            c3*c2*log(reactionRate1(idepth1,idepth2+1,idepth3+1,idepthDust)))
       nr2 = c1 * ((1.-c3)*(1.-c2)*log(reactionRate1(idepth1+1,idepth2,idepth3,idepthDust+1))+ &
            c3*(1.-c2)*log(reactionRate1(idepth1+1,idepth2,idepth3+1,idepthDust+1)) + &
            c2*(1.-c3)*log(reactionRate1(idepth1+1,idepth2+1,idepth3,idepthDust+1)) + &
            c3*c2*log(reactionRate1(idepth1+1,idepth2+1,idepth3+1,idepthDust+1))) &
            + (1.-c1) * ((1.-c3)*(1.-c2)*log(reactionRate1(idepth1,idepth2,idepth3,idepthDust+1)) + &
            c3*(1.-c2)*log(reactionRate1(idepth1,idepth2,idepth3+1,idepthDust+1)) + &
            c2*(1.-c3)*log(reactionRate1(idepth1,idepth2+1,idepth3,idepthDust+1)) + &
            c3*c2*log(reactionRate1(idepth1,idepth2+1,idepth3+1,idepthDust+1)))
       numberRate = exp((1.-cDust)*nr1 + cDust*nr2) ! [1/s]
       hr1 = c1 * ((1.-c3)*(1.-c2)*log(energyRate1(idepth1+1,idepth2,idepth3,idepthDust))+ &
            c3*(1.-c2)*log(energyRate1(idepth1+1,idepth2,idepth3+1,idepthDust)) + &
            c2*(1.-c3)*log(energyRate1(idepth1+1,idepth2+1,idepth3,idepthDust)) + &
            c3*c2*log(energyRate1(idepth1+1,idepth2+1,idepth3+1,idepthDust))) &
            + (1.-c1) * ((1.-c3)*(1.-c2)*log(energyRate1(idepth1,idepth2,idepth3,idepthDust)) + &
            c3*(1.-c2)*log(energyRate1(idepth1,idepth2,idepth3+1,idepthDust)) + &
            c2*(1.-c3)*log(energyRate1(idepth1,idepth2+1,idepth3,idepthDust)) + &
            c3*c2*log(energyRate1(idepth1,idepth2+1,idepth3+1,idepthDust)))
       hr2 = c1 * ((1.-c3)*(1.-c2)*log(energyRate1(idepth1+1,idepth2,idepth3,idepthDust+1))+ &
            c3*(1.-c2)*log(energyRate1(idepth1+1,idepth2,idepth3+1,idepthDust+1)) + &
            c2*(1.-c3)*log(energyRate1(idepth1+1,idepth2+1,idepth3,idepthDust+1)) + &
            c3*c2*log(energyRate1(idepth1+1,idepth2+1,idepth3+1,idepthDust+1))) &
            + (1.-c1) * ((1.-c3)*(1.-c2)*log(energyRate1(idepth1,idepth2,idepth3,idepthDust+1)) + &
            c3*(1.-c2)*log(energyRate1(idepth1,idepth2,idepth3+1,idepthDust+1)) + &
            c2*(1.-c3)*log(energyRate1(idepth1,idepth2+1,idepth3,idepthDust+1)) + &
            c3*c2*log(energyRate1(idepth1,idepth2+1,idepth3+1,idepthDust+1)))
       heatingRate = exp((1.-cDust)*hr1 + cDust*hr2) ! [erg/s]
    case(2)
       nr1 = c1 * ((1.-c3)*(1.-c2)*log(reactionRate2(idepth1+1,idepth2,idepth3,idepthDust))+ &
            c3*(1.-c2)*log(reactionRate2(idepth1+1,idepth2,idepth3+1,idepthDust)) + &
            c2*(1.-c3)*log(reactionRate2(idepth1+1,idepth2+1,idepth3,idepthDust)) + &
            c3*c2*log(reactionRate2(idepth1+1,idepth2+1,idepth3+1,idepthDust))) &
            + (1.-c1) * ((1.-c3)*(1.-c2)*log(reactionRate2(idepth1,idepth2,idepth3,idepthDust)) + &
            c3*(1.-c2)*log(reactionRate2(idepth1,idepth2,idepth3+1,idepthDust)) + &
            c2*(1.-c3)*log(reactionRate2(idepth1,idepth2+1,idepth3,idepthDust)) + &
            c3*c2*log(reactionRate2(idepth1,idepth2+1,idepth3+1,idepthDust)))
       nr2 = c1 * ((1.-c3)*(1.-c2)*log(reactionRate2(idepth1+1,idepth2,idepth3,idepthDust+1))+ &
            c3*(1.-c2)*log(reactionRate2(idepth1+1,idepth2,idepth3+1,idepthDust+1)) + &
            c2*(1.-c3)*log(reactionRate2(idepth1+1,idepth2+1,idepth3,idepthDust+1)) + &
            c3*c2*log(reactionRate2(idepth1+1,idepth2+1,idepth3+1,idepthDust+1))) &
            + (1.-c1) * ((1.-c3)*(1.-c2)*log(reactionRate2(idepth1,idepth2,idepth3,idepthDust+1)) + &
            c3*(1.-c2)*log(reactionRate2(idepth1,idepth2,idepth3+1,idepthDust+1)) + &
            c2*(1.-c3)*log(reactionRate2(idepth1,idepth2+1,idepth3,idepthDust+1)) + &
            c3*c2*log(reactionRate2(idepth1,idepth2+1,idepth3+1,idepthDust+1)))
       numberRate = exp((1.-cDust)*nr1 + cDust*nr2) ! [1/s]
       hr1 = c1 * ((1.-c3)*(1.-c2)*log(energyRate2(idepth1+1,idepth2,idepth3,idepthDust))+ &
            c3*(1.-c2)*log(energyRate2(idepth1+1,idepth2,idepth3+1,idepthDust)) + &
            c2*(1.-c3)*log(energyRate2(idepth1+1,idepth2+1,idepth3,idepthDust)) + &
            c3*c2*log(energyRate2(idepth1+1,idepth2+1,idepth3+1,idepthDust))) &
            + (1.-c1) * ((1.-c3)*(1.-c2)*log(energyRate2(idepth1,idepth2,idepth3,idepthDust)) + &
            c3*(1.-c2)*log(energyRate2(idepth1,idepth2,idepth3+1,idepthDust)) + &
            c2*(1.-c3)*log(energyRate2(idepth1,idepth2+1,idepth3,idepthDust)) + &
            c3*c2*log(energyRate2(idepth1,idepth2+1,idepth3+1,idepthDust)))
       hr2 = c1 * ((1.-c3)*(1.-c2)*log(energyRate2(idepth1+1,idepth2,idepth3,idepthDust+1))+ &
            c3*(1.-c2)*log(energyRate2(idepth1+1,idepth2,idepth3+1,idepthDust+1)) + &
            c2*(1.-c3)*log(energyRate2(idepth1+1,idepth2+1,idepth3,idepthDust+1)) + &
            c3*c2*log(energyRate2(idepth1+1,idepth2+1,idepth3+1,idepthDust+1))) &
            + (1.-c1) * ((1.-c3)*(1.-c2)*log(energyRate2(idepth1,idepth2,idepth3,idepthDust+1)) + &
            c3*(1.-c2)*log(energyRate2(idepth1,idepth2,idepth3+1,idepthDust+1)) + &
            c2*(1.-c3)*log(energyRate2(idepth1,idepth2+1,idepth3,idepthDust+1)) + &
            c3*c2*log(energyRate2(idepth1,idepth2+1,idepth3+1,idepthDust+1)))
       heatingRate = exp((1.-cDust)*hr1 + cDust*hr2) ! [erg/s]
    case(3)
       nr1 = c1 * ((1.-c3)*(1.-c2)*log(reactionRate3(idepth1+1,idepth2,idepth3,idepthDust))+ &
            c3*(1.-c2)*log(reactionRate3(idepth1+1,idepth2,idepth3+1,idepthDust)) + &
            c2*(1.-c3)*log(reactionRate3(idepth1+1,idepth2+1,idepth3,idepthDust)) + &
            c3*c2*log(reactionRate3(idepth1+1,idepth2+1,idepth3+1,idepthDust))) &
            + (1.-c1) * ((1.-c3)*(1.-c2)*log(reactionRate3(idepth1,idepth2,idepth3,idepthDust)) + &
            c3*(1.-c2)*log(reactionRate3(idepth1,idepth2,idepth3+1,idepthDust)) + &
            c2*(1.-c3)*log(reactionRate3(idepth1,idepth2+1,idepth3,idepthDust)) + &
            c3*c2*log(reactionRate3(idepth1,idepth2+1,idepth3+1,idepthDust)))
       nr2 = c1 * ((1.-c3)*(1.-c2)*log(reactionRate3(idepth1+1,idepth2,idepth3,idepthDust+1))+ &
            c3*(1.-c2)*log(reactionRate3(idepth1+1,idepth2,idepth3+1,idepthDust+1)) + &
            c2*(1.-c3)*log(reactionRate3(idepth1+1,idepth2+1,idepth3,idepthDust+1)) + &
            c3*c2*log(reactionRate3(idepth1+1,idepth2+1,idepth3+1,idepthDust+1))) &
            + (1.-c1) * ((1.-c3)*(1.-c2)*log(reactionRate3(idepth1,idepth2,idepth3,idepthDust+1)) + &
            c3*(1.-c2)*log(reactionRate3(idepth1,idepth2,idepth3+1,idepthDust+1)) + &
            c2*(1.-c3)*log(reactionRate3(idepth1,idepth2+1,idepth3,idepthDust+1)) + &
            c3*c2*log(reactionRate3(idepth1,idepth2+1,idepth3+1,idepthDust+1)))
       numberRate = exp((1.-cDust)*nr1 + cDust*nr2) ! [1/s]
       hr1 = c1 * ((1.-c3)*(1.-c2)*log(energyRate3(idepth1+1,idepth2,idepth3,idepthDust))+ &
            c3*(1.-c2)*log(energyRate3(idepth1+1,idepth2,idepth3+1,idepthDust)) + &
            c2*(1.-c3)*log(energyRate3(idepth1+1,idepth2+1,idepth3,idepthDust)) + &
            c3*c2*log(energyRate3(idepth1+1,idepth2+1,idepth3+1,idepthDust))) &
            + (1.-c1) * ((1.-c3)*(1.-c2)*log(energyRate3(idepth1,idepth2,idepth3,idepthDust)) + &
            c3*(1.-c2)*log(energyRate3(idepth1,idepth2,idepth3+1,idepthDust)) + &
            c2*(1.-c3)*log(energyRate3(idepth1,idepth2+1,idepth3,idepthDust)) + &
            c3*c2*log(energyRate3(idepth1,idepth2+1,idepth3+1,idepthDust)))
       hr2 = c1 * ((1.-c3)*(1.-c2)*log(energyRate3(idepth1+1,idepth2,idepth3,idepthDust+1))+ &
            c3*(1.-c2)*log(energyRate3(idepth1+1,idepth2,idepth3+1,idepthDust+1)) + &
            c2*(1.-c3)*log(energyRate3(idepth1+1,idepth2+1,idepth3,idepthDust+1)) + &
            c3*c2*log(energyRate3(idepth1+1,idepth2+1,idepth3+1,idepthDust+1))) &
            + (1.-c1) * ((1.-c3)*(1.-c2)*log(energyRate3(idepth1,idepth2,idepth3,idepthDust+1)) + &
            c3*(1.-c2)*log(energyRate3(idepth1,idepth2,idepth3+1,idepthDust+1)) + &
            c2*(1.-c3)*log(energyRate3(idepth1,idepth2+1,idepth3,idepthDust+1)) + &
            c3*c2*log(energyRate3(idepth1,idepth2+1,idepth3+1,idepthDust+1)))
       heatingRate = exp((1.-cDust)*hr1 + cDust*hr2) ! [erg/s]
    endselect

  end subroutine getRatesHydrogenHelium

  function enclosedSphericalSurface(x0,y0,z0,r) result(surface)

    use definitions

    implicit none
    real(kind=RealKind), intent(in) :: x0, y0, z0, r
    real(kind=RealKind) :: surface, theta, phi, x, y, z
    integer :: level, i, nside
    integer*8 :: longIntegerArgument
    integer, parameter :: ntheta = 30, nphi = 30

    surface = 0.

    level = 6
    nside = 2**(level-1)
    do i = 1, 12*4**(level-1)
       longIntegerArgument = i-1
       !     print*, '>>>', longIntegerArgument
       call pix2ang_nest(nside,longIntegerArgument,phi,theta)
       x = x0 + r*cos(theta)*cos(phi)
       y = y0 + r*cos(theta)*sin(phi)
       z = z0 + r*sin(theta)
       if (x.ge.0. .and. x.le.1. .and. y.ge.0. .and. y.le.1. .and. z.ge.0. .and. z.le.1.) &
            surface = surface + cos(theta)
    enddo
    surface = surface / float(12*4**(level-1)) / quarterPi

    if (surface.gt.1.) surface = 1.

  end function enclosedSphericalSurface

  recursive subroutine saveRestoreOriginalFields(currentCell)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    integer :: i, j, k

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call saveRestoreOriginalFields(currentCell%cell(i,j,k))
             enddo
          enddo
       enddo
    else
       !     currentCell%tgasprev = currentCell%tgas
       currentCell%HIprev = currentCell%HI
       currentCell%HeIprev = currentCell%HeI
       currentCell%HeIIprev = currentCell%HeII
    endif

  end subroutine saveRestoreOriginalFields

  recursive subroutine computeMass(currentCell,nx)

    use definitions

    implicit none
    type(zoneType), target :: currentCell
    integer, intent(in) :: nx
    integer :: i, j, k
    real(kind=RealKind) :: physicalCellSizeCube

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call computeMass(currentCell%cell(i,j,k),nx)
             enddo
          enddo
       enddo
    else
       physicalCellSizeCube = (physicalBoxSize / (float(2**currentCell%level)*float(nx)))**3
       neutralHydrogenMass = neutralHydrogenMass + currentCell%HI * mh * physicalCellSizeCube / msun
       totalHydrogenMass = totalHydrogenMass + psi * currentCell%rho * physicalCellSizeCube / msun
    endif

  end subroutine computeMass

  recursive subroutine computeExpansionParameters(nh)

    use definitions

    implicit none
    real(kind=RealKind), intent(in) :: nh
    integer :: i
    real(kind=RealKind) :: lognh, tmp
    integer, parameter :: nExpansionTable = 10
    real(kind=RealKind), dimension(nExpansionTable) :: logInitialDensity, logFinalDensity, logFinalRadius

    logInitialDensity = (/ 0.00000, 0.333333, 0.666667,  1.00000,  1.33333,  1.66667,  2.00000,  2.33333,  2.66667,  3.00000 /)
    logFinalRadius = (/ 2.99506, 2.77808, 2.57210, 2.37683, 2.19731, 2.02898, 1.87315, 1.73656, 1.61294, 1.50202 /)
    logFinalDensity = (/ -0.0222764, 0.295050, 0.579490, 0.831870,  1.03717,  1.20892,  1.34321,  1.41970,  1.45725,  1.45667 /)

    sourceTotalHydrogenDensity = nh

    lognh = log10(nh)
    i = 1
    do while (lognh.gt.logInitialDensity(i) .and. i.lt.nExpansionTable)
       i = i + 1
    enddo
    i = max(i,2)
    tmp = (lognh-logInitialDensity(i-1))/(logInitialDensity(i)-logInitialDensity(i-1))
    finalRadius = 10.**(tmp*(logFinalRadius(i)-logFinalRadius(i-1))+logFinalRadius(i-1)) * pc
    densityCoefficient = 10.**(tmp*(logFinalDensity(i)-logFinalDensity(i-1))+logFinalDensity(i-1)) / nh

    if (lognh.lt.logInitialDensity(1)) then
       tmp = (lognh+6.)/(logInitialDensity(1)+6.)
       densityCoefficient = 10.**(tmp*(logFinalDensity(1)+6.)-6.) / nh
    endif

    !  print*, lognh, finalRadius/kpc, densityCoefficient

  end subroutine computeExpansionParameters

  recursive subroutine findExpansion(currentCell,xcell,ycell,zcell,nx)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    real(kind=RealKind), intent(in) :: xcell, ycell, zcell
    integer, intent(in) :: nx
    integer :: i, j, k
    real(kind=RealKind) :: xnew, ynew, znew, shift, dist

    if (currentCell%refined) then
       shift = 0.25 / (float(2**currentCell%level)*float(nx))
       do i = 1, 2
          if (i.eq.1) then
             xnew = xcell - shift
          else
             xnew = xcell + shift
          endif
          do j = 1, 2
             if (j.eq.1) then
                ynew = ycell - shift
             else
                ynew = ycell + shift
             endif
             do k = 1, 2
                if (k.eq.1) then
                   znew = zcell - shift
                else
                   znew = zcell + shift
                endif
                call findExpansion(currentCell%cell(i,j,k),xnew,ynew,znew,nx)
             enddo
          enddo
       enddo
    else
       dist = physicalBoxSize * sqrt((xbase-xcell)**2+(ybase-ycell)**2+(zbase-zcell)**2)
       !     if (dist.lt.finalRadius .and. psi*currentCell%rho/mh.lt.sourceTotalHydrogenDensity) then
       if (dist.lt.finalRadius .and. psi*currentCell%rho/mh.le.1.0001*sourceTotalHydrogenDensity) then
          currentCell%rhoCoef = min(currentCell%rhoCoef,densityCoefficient)
       endif
    endif

  end subroutine findExpansion

  recursive subroutine applyExpansion(currentCell)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    integer :: i, j, k

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call applyExpansion(currentCell%cell(i,j,k))
             enddo
          enddo
       enddo
    else
       ! icosmic = icosmic + 1
       ! ddd = ddd + currentCell%rhoCoef
       if (currentCell%rhoCoef.lt.1.) then
          currentCell%rho  = currentCell%rho  * currentCell%rhoCoef
          currentCell%HI   = currentCell%HI   * currentCell%rhoCoef
          currentCell%HeI  = currentCell%HeI  * currentCell%rhoCoef
          currentCell%HeII = currentCell%HeII * currentCell%rhoCoef
       endif
    endif

  end subroutine applyExpansion

  recursive subroutine addStarToTree(currentNode,x0,y0,z0,level,starWeight)

    use definitions

    implicit none
    type(nodeType), target :: currentNode
    real(kind=RealKind), intent(in) :: x0, y0, z0
    integer, intent(in) :: level, starWeight
    integer :: i, j, k, oldWeight
    real(kind=RealKind) :: xnew, ynew, znew

    !  print*, 'entered node of level', level, 'with', currentNode%n, 'particles'

    oldWeight = currentNode%weight

    currentNode%n = currentNode%n + 1
    currentNode%weight = currentNode%weight + starWeight

    if (currentNode%n.eq.1) then
       currentNode%xpos = x0
       currentNode%ypos = y0
       currentNode%zpos = z0
    endif

    if (currentNode%n.eq.2) then ! add the first (already existing) star to the refined node
       allocate(currentNode%node(2,2,2))
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call nullifyNode(currentNode%node(i,j,k))
                currentNode%node(i,j,k)%size = 0.5**(level+1)
                currentNode%node(i,j,k)%parent => currentNode
             enddo
          enddo
       enddo
       call selectChild(currentNode%xpos,currentNode%ypos,currentNode%zpos,i,j,k,xnew,ynew,znew)
       !     write(*,1031) level+1, i, j, k
       call addStarToTree(currentNode%node(i,j,k),xnew,ynew,znew,level+1,oldWeight)
    endif

    if (currentNode%n.ge.2) then ! add the current star to the refined node
       call selectChild(x0,y0,z0,i,j,k,xnew,ynew,znew)
       !     write(*,1032) level+1, i, j, k
       call addStarToTree(currentNode%node(i,j,k),xnew,ynew,znew,level+1,starWeight)
    endif

    if (currentNode%n.ge.2) then ! add the coordinates of the current star to the node
       currentNode%xpos = (currentNode%xpos*oldWeight+x0*starWeight)/currentNode%weight
       currentNode%ypos = (currentNode%ypos*oldWeight+y0*starWeight)/currentNode%weight
       currentNode%zpos = (currentNode%zpos*oldWeight+z0*starWeight)/currentNode%weight
    endif

1031 format('existing->refine',i5,'     ',3i5)
1032 format('new->refine',i5,'     ',3i5)

  end subroutine addStarToTree

  recursive subroutine countNodes(currentNode,level)

    use definitions

    implicit none
    type(nodeType) :: currentNode
    integer, intent(in) :: level
    integer :: i, j, k

    if (level.gt.highestNodeLevel) highestNodeLevel = level

    if (currentNode%n.gt.0) nNodes = nNodes + 1

    if (currentNode%n.ge.2) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call countNodes(currentNode%node(i,j,k),level+1)
             enddo
          enddo
       enddo
    endif

  end subroutine countNodes

  subroutine selectChild(x0,y0,z0,i,j,k,xnew,ynew,znew)

    use definitions

    implicit none
    real(kind=RealKind), intent(in) :: x0, y0, z0
    integer, intent(out) :: i, j, k
    real(kind=RealKind), intent(out) :: xnew, ynew, znew

    if (x0.le.0.5) then
       i = 1
       xnew = 2.*x0
    else
       i = 2
       xnew = 2.*x0 - 1.
    endif
    if (y0.le.0.5) then
       j = 1
       ynew = 2.*y0
    else
       j = 2
       ynew = 2.*y0 - 1.
    endif
    if (z0.le.0.5) then
       k = 1
       znew = 2.*z0
    else
       k = 2
       znew = 2.*z0 - 1.
    endif

  end subroutine selectChild

  subroutine nullifyNode(currentNode)

    use definitions

    implicit none
    type(nodeType) :: currentNode

    currentNode%n = 0
    currentNode%weight = 0
    currentNode%xpos = 0.
    currentNode%ypos = 0.
    currentNode%zpos = 0.
    currentNode%size = 0.

  end subroutine nullifyNode

  recursive subroutine convertPosToAbsoluteCoordinates(currentNode,level,x1,y1,z1)

    use definitions

    implicit none
    type(nodeType) :: currentNode
    integer, intent(in) :: level
    real(kind=RealKind), intent(in) :: x1, y1, z1
    integer :: i, j, k
    real(kind=RealKind) :: oneOverZoom, xnew, ynew, znew

    oneOverZoom = 1./float(2**level)

    if (level.ge.1 .and. currentNode%weight.gt.0) then
       currentNode%xpos = x1 + currentNode%xpos * oneOverZoom
       currentNode%ypos = y1 + currentNode%ypos * oneOverZoom
       currentNode%zpos = z1 + currentNode%zpos * oneOverZoom
    endif

    if (currentNode%n.ge.2) then
       do i = 1, 2
          if (i.eq.1) then
             xnew = x1
          else
             xnew = x1 + 0.5 * oneOverZoom
          endif
          do j = 1, 2
             if (j.eq.1) then
                ynew = y1
             else
                ynew = y1 + 0.5 * oneOverZoom
             endif
             do k = 1, 2
                if (k.eq.1) then
                   znew = z1
                else
                   znew = z1 + 0.5 * oneOverZoom
                endif
                call convertPosToAbsoluteCoordinates(currentNode%node(i,j,k),level+1,xnew,ynew,znew)
             enddo
          enddo
       enddo
    endif

  end subroutine convertPosToAbsoluteCoordinates

  recursive subroutine computeGasPDF(currentCell)

    use definitions

    implicit none
    type(zoneType), target :: currentCell
    integer :: i, j, k, itmp
    real(kind=RealKind) :: tmp

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call computeGasPDF(currentCell%cell(i,j,k))
             enddo
          enddo
       enddo
    else
       tmp = dlog10(currentCell%rho/msun*pc**3)
       if (tmp.gt.apdf .and. tmp.lt.bpdf) then
          itmp = int((tmp-apdf)/(bpdf-apdf)*float(npdf)) + 1
          pdfGas(itmp) = pdfGas(itmp) + 2.**(-3*currentCell%level)
       else
          pdfGasOutside = pdfGasOutside + 2.**(-3*currentCell%level)
       endif
    endif

  end subroutine computeGasPDF

  recursive subroutine computeClumping(currentCell)

    use definitions

    implicit none
    type(zoneType), target :: currentCell
    integer :: i, j, k
    real(kind=RealKind) :: nh

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call computeClumping(currentCell%cell(i,j,k))
             enddo
          enddo
       enddo
    else
       nh = psi * currentCell%rho / mh
       volSum = volSum + 1. / 2.**(3*currentCell%level)
       volSumDensity = volSumDensity + nh / 2.**(3*currentCell%level)
       volSumDensity2 = volSumDensity2 + nh**2 / 2.**(3*currentCell%level)
    endif

  end subroutine computeClumping


  recursive subroutine readLatestIonization(currentCell,level)

    use definitions

    implicit none
    type(zoneType), target :: currentCell
    integer, intent(in) :: level
    integer :: i, j, k
    real(kind=RealKind) :: nh, nhe, tmp

    icosmic = icosmic + 1
    if (cellArrayLevel(icosmic).eq.level) then
       currentCell%HI = dble(cellArrayHI(icosmic))
       currentCell%HeI = dble(cellArrayHeI(icosmic))
       currentCell%HeII = dble(cellArrayHeII(icosmic))
       nh = psi*currentCell%rho / mh
       nhe = (1.-psi)*currentCell%rho / mhe

       ! if (currentCell%HI.gt.1.001*nh .or. (currentCell%HeI+currentCell%HeII).gt.1.001*nhe .or. &
       !      min(currentCell%HI,currentCell%HeI,currentCell%HeII).lt.0.d0) then
       !    print*, '-----------------------'
       !    print*, 'nh =', nh
       !    print*, 'nhe =', nhe
       !    print*, currentCell%HI/nh, currentCell%HeI/nhe, currentCell%HeII/nhe
       !    print*, 'rho/temp =', currentCell%rho, cellArrayTemp(icosmic)
       ! endif

       currentCell%HI = max(currentCell%HI,0.d0)
       currentCell%HeI = max(currentCell%HeI,0.d0)
       currentCell%HeII = max(currentCell%HeII,0.d0)
       if (currentCell%HI.gt.nh) currentCell%HI = nh
       if ((currentCell%HeI+currentCell%HeII).gt.nhe) then
          tmp = (currentCell%HeI+currentCell%HeII)/nhe
          currentCell%HeI = currentCell%HeI / tmp
          currentCell%HeII = currentCell%HeII / tmp
       endif
       currentCell%tgas = dble(cellArrayTemp(icosmic))
    else
       if (cellArrayLevel(icosmic).gt.level) then
          if (.not.associated(currentCell%cell)) then
             write(*,*) 'error in readLatestIonization'
             stop
          endif
          icosmic = icosmic - 1
          do i = 1, 2
             do j = 1, 2
                do k = 1, 2
                   call readLatestIonization(currentCell%cell(i,j,k),level+1)
                enddo
             enddo
          enddo
       else
          write(*,*) 'error in levels', icosmic, cellArrayLevel(icosmic), level
          stop
       endif
    endif

  end subroutine readLatestIonization

  subroutine writeIonization(itime,nx,ny,nz)

    use definitions

    implicit none
    integer, intent(in) :: itime, nx, ny, nz
    integer :: i, j, k
    character(4) :: suffix
    integer :: sfstart, sffinfo, sfselect, sfginfo, sfrdata, sfendacc, sfend, sfcreate, sfwdata
    integer, dimension(3) :: baseGridSize

    ! ---- store hierarchy with space-filling curve in HDF4: Jmean, HI ----

    icosmic = 0
    do i = 1, nx
       do j = 1, ny
          do k = 1, nz
             call countCells(baseGrid%cell(i,j,k))
          enddo
       enddo
    enddo
    write(*,*) icosmic

    allocate (cellArrayLevel(icosmic))
    allocate (cellArrayHI(icosmic))
    allocate (cellArrayHeI(icosmic))
    allocate (cellArrayHeII(icosmic))
    allocate (cellArrayTemp(icosmic))
    allocate (cellArrayDensity(icosmic))
    if (readKinematics) allocate (cellArrayVelx(icosmic),cellArrayVely(icosmic),cellArrayVelz(icosmic))
    if (readMetals) allocate (cellArrayAbun2(icosmic))

    icosmic = 0
    do i = 1, nx
       do j = 1, ny
          do k = 1, nz
             call writeCell(baseGrid%cell(i,j,k),0)
          enddo
       enddo
    enddo

    write(suffix,fmt="(i4)") itime
    do i = 1, 3
       if(suffix(i:i) == " ") suffix(i:i) = "0"
    end do

    sd_id = sfstart('cellArray'//suffix//'.h4',dfacc_create)

    baseGridSize = (/nx,ny,nz/)

    edges(1) = 3
    start = 0
    stride = 1

    sds_id = sfcreate(sd_id,'base grid size',dfnt_int32,1,edges)
    status = sfwdata(sds_id,start,stride,edges,baseGridSize)
    status = sfendacc(sds_id)

    edges(1) = icosmic

    sds_id = sfcreate(sd_id,'level',dfnt_int32,1,edges)
    status = sfwdata(sds_id,start,stride,edges,cellArrayLevel)
    status = sfendacc(sds_id)

    sds_id = sfcreate(sd_id,'HI',dfnt_float32,1,edges)
    status = sfwdata(sds_id,start,stride,edges,cellArrayHI)
    status = sfendacc(sds_id)

    sds_id = sfcreate(sd_id,'HeI',dfnt_float32,1,edges)
    status = sfwdata(sds_id,start,stride,edges,cellArrayHeI)
    status = sfendacc(sds_id)

    sds_id = sfcreate(sd_id,'HeII',dfnt_float32,1,edges)
    status = sfwdata(sds_id,start,stride,edges,cellArrayHeII)
    status = sfendacc(sds_id)

    sds_id = sfcreate(sd_id,'temperature',dfnt_float32,1,edges)
    status = sfwdata(sds_id,start,stride,edges,cellArrayTemp)
    status = sfendacc(sds_id)

    sds_id = sfcreate(sd_id,'density',dfnt_float32,1,edges)
    status = sfwdata(sds_id,start,stride,edges,cellArrayDensity)
    status = sfendacc(sds_id)

    if (readKinematics) then

       sds_id = sfcreate(sd_id,'velx',dfnt_float32,1,edges)
       status = sfwdata(sds_id,start,stride,edges,cellArrayVelx)
       status = sfendacc(sds_id)

       sds_id = sfcreate(sd_id,'vely',dfnt_float32,1,edges)
       status = sfwdata(sds_id,start,stride,edges,cellArrayVely)
       status = sfendacc(sds_id)

       sds_id = sfcreate(sd_id,'velz',dfnt_float32,1,edges)
       status = sfwdata(sds_id,start,stride,edges,cellArrayVelz)
       status = sfendacc(sds_id)

    endif

    if (readMetals) then

       sds_id = sfcreate(sd_id,'abun2',dfnt_float32,1,edges)
       status = sfwdata(sds_id,start,stride,edges,cellArrayAbun2)
       status = sfendacc(sds_id)

    endif

    status = sfend(sd_id)

    deallocate(cellArrayLevel,cellArrayHI,cellArrayHeI,cellArrayHeII, &
         cellArrayTemp,cellArrayDensity)
    if (readKinematics) deallocate(cellArrayVelx,cellArrayVely,cellArrayVelz)
    if (readMetals) deallocate(cellArrayAbun2)

  end subroutine writeIonization

  recursive subroutine projectVariableToMap(currentCell,currentPixel, &
       x0,y0,cellSizeAbsoluteUnits,is,js,ks)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    real*4, intent(inout) :: currentPixel
    real(kind=RealKind), intent(in) :: x0, y0
    real(kind=RealKind), intent(in) :: cellSizeAbsoluteUnits
    integer, dimension(2,2,2), intent(in) :: is, js, ks
    integer :: i, j, k
    real(kind=RealKind) :: xnew, ynew, tmp

    if (currentCell%refined) then
       if (x0.lt.0.5) then
          i = 1
          xnew = 2.*x0
       else
          i = 2
          xnew = 2.*x0 - 1.
       endif
       if (y0.lt.0.5) then
          j = 1
          ynew = 2.*y0
       else
          j = 2
          ynew = 2.*y0 - 1.
       endif
       do k = 1, 2
          call projectVariableToMap(currentCell%cell(is(i,j,k),js(i,j,k),ks(i,j,k)),currentPixel, &
               xnew,ynew,cellSizeAbsoluteUnits/2.,is,js,ks)
       enddo
    else
       tmp = (currentCell%rho*1.e25) * (cellSizeAbsoluteUnits/1.e21)**3
       currentPixel = currentPixel + currentCell%abun2 * tmp
       totalMass = totalMass + tmp
!        print*, currentCell%abun2
    endif

  end subroutine projectVariableToMap

  recursive subroutine computeOpacities(currentCell)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    integer :: i, j, k

    currentCell%Jmean1 = 0.
    currentCell%Jmean2 = 0.
    currentCell%Jmean3 = 0.

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call computeOpacities(currentCell%cell(i,j,k))
             enddo
          enddo
       enddo
    else
       currentCell%kappa1 = currentCell%HI * group1%beta24
       currentCell%kappa2 = currentCell%HI * group2%beta24 + currentCell%HeI * group2%beta26
       currentCell%kappa3 = currentCell%HI * group3%beta24 + &
            currentCell%HeI * group3%beta26 + currentCell%HeII * group3%beta25
    endif

  end subroutine computeOpacities
  
  subroutine powerSpectrumIndex(uvb1,alpha1, &
       uvb2,alpha2,uvbTotal,alphaTotal,nug,nugplus,bound)

    use definitions, only : RealKind

    implicit none
    real(kind=RealKind), intent(in)  :: uvb1, alpha1, &
         uvb2, alpha2, nug, nugplus
    logical, intent(in) :: bound
    real(kind=RealKind), intent(out) :: uvbTotal, alphaTotal
    real(kind=RealKind) :: tmp, tmp1, tmp2, tmpold, fun, fun1, fun2, funToMatch

    uvbTotal = uvb1 + uvb2

    if (bound) then
       funToMatch = uvb1/(alpha1-1.)*(1.-(nug/nugplus)**(alpha1-1.)) + &
            uvb2/(alpha2-1.)*(1.-(nug/nugplus)**(alpha2-1.))
    else
       funToMatch = uvb1/(alpha1-1.) + uvb2/(alpha2-1.)
    endif

    tmp1 = 1.1 * alpha1 - 0.1 * alpha2
    tmp2 = 1.1 * alpha2 - 0.1 * alpha1
    if (bound) then
       fun1 = uvbTotal/(tmp1-1.)*(1.-(nug/nugplus)**(tmp1-1.)) - funToMatch
       fun2 = uvbTotal/(tmp2-1.)*(1.-(nug/nugplus)**(tmp2-1.)) - funToMatch
    else
       fun1 = uvbTotal/(tmp1-1.) - funToMatch
       fun2 = uvbTotal/(tmp2-1.) - funToMatch
    endif

    if (.not.opposite(fun1,fun2)) then
       print*, 'wrong sign', fun1, fun2
       stop
    endif

    tmpold = tmp1
    tmp = tmp2
    do while (abs(tmp-tmpold).ge.1e-8)
       tmpold = tmp
       tmp = (tmp1*abs(fun2)+tmp2*abs(fun1))/(abs(fun1)+abs(fun2))
       if (bound) then
          fun = uvbTotal/(tmp-1.)*(1.-(nug/nugplus)**(tmp-1.)) - funToMatch
       else
          fun = uvbTotal/(tmp-1.) - funToMatch
       endif
       if (opposite(fun,fun1)) then
          tmp2 = tmp
          fun2 = fun
       else
          tmp1 = tmp
          fun1 = fun
       endif
    enddo

    alphaTotal = tmp

  end subroutine powerSpectrumIndex

  function opposite(a,b)

    use definitions, only : RealKind

    implicit none
    real(kind=RealKind) :: a, b
    logical :: opposite

    if (((a > 0.).and.(b < 0.)).or.((a < 0.).and.(b > 0.))) then
       opposite = .true.
    else
       opposite = .false.
    endif

  end function opposite

end program pointTransfer
