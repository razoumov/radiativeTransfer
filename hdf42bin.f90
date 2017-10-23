module localDefinitions

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
  real(kind=RealKind), parameter :: psi = 0.76

  integer :: icosmic, ncosmic, itmp

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
  real*4, dimension(:), pointer :: cellArrayXpos, cellArrayYpos, cellArrayZpos, &
       cellArrayHI, cellArrayHeI, cellArrayHeII, &
       cellArrayTemp, cellArrayDensity, cellArrayVelx, cellArrayVely, cellArrayVelz

end module localDefinitions

program hdf42bin

  use localDefinitions

  implicit none
  integer :: i, j, k, icell, jcell, kcell, nx, ny, nz, &
       i0, j0, k0, kstart, kend
  real(kind=RealKind) :: xnew, ynew, znew, x0, y0, z0, tmp, physicalBoxSize, &
       xa, xb, ya, yb, za, zb, xpos, ypos, zpos
  integer, dimension(3) :: baseGridSize
  integer, dimension(2,2,2) :: is, js, ks
  integer, dimension(0:30) :: numberCells

  integer sfstart, sffinfo, sfselect, sfginfo, sfrdata, sfendacc, sfend, sfcreate, sfwdata
  real*4, dimension(:,:), pointer :: mapArray
  character(60) :: dirname, basename

  dirname = '11.63/b10_0.9/'
  physicalBoxSize = 1200.*kpc
  basename = 'cellArray0047'

  sd_id = sfstart(trim(dirname)//trim(basename)//'.h4', dfacc_read)
  status = sffinfo(sd_id, n_datasets, n_file_attrs)

  start = 0
  edges(1) = 3
  stride = 1

  sds_id = sfselect(sd_id, 0)
  status = sfrdata(sds_id, start, stride, edges, baseGridSize)
  status = sfendacc(sds_id)

  write(*,*) 'base grid size =', baseGridSize
  nx = baseGridSize(1)
  ny = baseGridSize(2)
  nz = baseGridSize(3)

  kstart = 1
  kend = nz

  sds_id = sfselect(sd_id, 1)
  status = sfginfo(sds_id,sds_name,rank,dim_sizes,data_type,n_attrs)
  if (status.ne.0 .or. rank.ne.1) then
     write(*,*) 'error reading ...'
     stop
  endif
  ncosmic = dim_sizes(1)
  edges(1) = dim_sizes(1)

  print*, 'number of cells =', ncosmic

  allocate (cellArrayLevel(ncosmic))
  allocate (cellArrayXpos(ncosmic))
  allocate (cellArrayYpos(ncosmic))
  allocate (cellArrayZpos(ncosmic))
  allocate (cellArrayHI(ncosmic))
  allocate (cellArrayHeI(ncosmic))
  allocate (cellArrayHeII(ncosmic))
  allocate (cellArrayTemp(ncosmic)) ! can erase this one
  allocate (cellArrayDensity(ncosmic)) ! can erase this one

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

  sds_id = sfselect(sd_id, 6)
  status = sfrdata(sds_id, start, stride, edges, cellArrayDensity)
  status = sfendacc(sds_id)

  status = sfend(sd_id)

  icosmic = 0
  do i = 1, nx
     xpos = (float(i)-0.5) / float(nx)
     do j = 1, ny
        ypos = (float(j)-0.5) / float(ny)
        do k = 1, nz
           zpos = (float(k)-0.5) / float(nz)
           call computeCellCoordinates(xpos,ypos,zpos,0,1.d0/dfloat(nx))
        enddo
     enddo
  enddo

  xa = - 0.5*physicalBoxSize/kpc
  xb = 0.5*physicalBoxSize/kpc

  ya = - 0.5*physicalBoxSize/kpc
  yb = 0.5*physicalBoxSize/kpc

  za = - 0.5*physicalBoxSize/kpc
  zb = 0.5*physicalBoxSize/kpc

  itmp = 0
  do icosmic = 1, ncosmic
     cellArrayXpos(icosmic) = cellArrayXpos(icosmic)*(xb-xa) + xa
     cellArrayYpos(icosmic) = cellArrayYpos(icosmic)*(yb-ya) + ya
     cellArrayZpos(icosmic) = cellArrayZpos(icosmic)*(zb-za) + za
     itmp = max(itmp,cellArrayLevel(icosmic))
  enddo
  print*, 'max level of refinement =', itmp
  print*, 'base level =', cellArrayLevel(1)

!   do icosmic = 1, ncosmic
!      if (cellArrayLevel(icosmic).eq.7) then
!         write(*,1010) icosmic, cellArrayLevel(icosmic), &
!              cellArrayXpos(icosmic), cellArrayYpos(icosmic), cellArrayZpos(icosmic)
!      endif
!   enddo
! 1010 format(i10,i5,3f23.14)

  numberCells = 0
  do icosmic = 1, ncosmic
     numberCells(cellArrayLevel(icosmic)) = numberCells(cellArrayLevel(icosmic)) + 1
  enddo
  print*, numberCells(0:15)

  open(14,file=trim(dirname)//trim(basename)//'.dat',status='replace',form='unformatted')
  write(14) (cellArrayLevel(i),i=1,ncosmic)
  write(14) (cellArrayXpos(i),i=1,ncosmic)
  write(14) (cellArrayYpos(i),i=1,ncosmic)
  write(14) (cellArrayZpos(i),i=1,ncosmic)
  write(14) (cellArrayHI(i),i=1,ncosmic)
  write(14) (cellArrayHeI(i),i=1,ncosmic)
  write(14) (cellArrayHeII(i),i=1,ncosmic)
  write(14) (cellArrayTemp(i),i=1,ncosmic)
  write(14) (cellArrayDensity(i),i=1,ncosmic)
  close(14)

end program hdf42bin

recursive subroutine computeCellCoordinates(x0,y0,z0,level,cellSize)

  use localDefinitions

  implicit none
  integer, intent(in) :: level
  real(kind=RealKind), intent(in) :: x0, y0, z0, cellSize
  integer :: i, j, k
  real(kind=RealKind) :: xnew, ynew, znew

  icosmic = icosmic + 1

  if (cellArrayLevel(icosmic).eq.level) then
     cellArrayXpos(icosmic) = x0
     cellArrayYpos(icosmic) = y0
     cellArrayZpos(icosmic) = z0
  else
     if (cellArrayLevel(icosmic).gt.level) then
        icosmic = icosmic - 1
        do i = 1, 2
           if (i.eq.1) then
              xnew = x0 - 0.25*cellSize
           else
              xnew = x0 + 0.25*cellSize
           endif
           do j = 1, 2
              if (j.eq.1) then
                 ynew = y0 - 0.25*cellSize
              else
                 ynew = y0 + 0.25*cellSize
              endif
              do k = 1, 2
                 if (k.eq.1) then
                    znew = z0 - 0.25*cellSize
                 else
                    znew = z0 + 0.25*cellSize
                 endif
                 call computeCellCoordinates(xnew,ynew,znew,level+1,cellSize/2.d0)
              enddo
           enddo
        enddo
     else
        write(*,*) 'error in levels', icosmic, cellArrayLevel(icosmic), level
        stop
     endif
  endif

end subroutine computeCellCoordinates
