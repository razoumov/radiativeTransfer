module localDefinitions

  implicit none
  integer, parameter :: SingleKind = kind(1.0)
  integer, parameter :: RealKind = kind(1.0d0)
  real(kind=RealKind), parameter :: pc = 3.08568025e18
  real(kind=RealKind), parameter :: kpc = 1.e3*pc
  real(kind=RealKind), parameter :: Mpc = 1.e6*pc
  real(kind=RealKind), parameter :: mp = 1.6726231e-24
  real(kind=RealKind), parameter :: mn = 1.67492728e-24
  real(kind=RealKind), parameter :: mh = mp
  real(kind=RealKind), parameter :: mhe = 2.*(mp+mn)
  real(kind=RealKind), parameter :: psi = 0.76

  integer :: icosmic, ncosmic

  type :: zoneType
     real*4 :: rho, tgas, HI, HeI, HeII
     logical(1) :: refined
     integer(1) :: level ! maximum 127 levels of cell refinement
     type(zoneType), pointer :: parent
     type(zoneType), dimension(:,:,:), pointer :: cell
  end type zoneType

  type(zoneType), target :: baseGrid

  integer*4, dimension(:), pointer :: cellArrayLevel
  real*4, dimension(:), pointer :: cellArrayXpos, cellArrayYpos, cellArrayZpos, &
       cellArrayHI, cellArrayHeI, cellArrayHeII, cellArrayTemp, cellArrayDensity

  type :: readLevelType
     integer :: ncell
     real*4, dimension(:,:), pointer :: pos, vel, abun
     real*4, dimension(:), pointer :: lT, lnH, lx
  end type readLevelType

end module localDefinitions

program convertFormats

  use localDefinitions

  implicit none
  integer :: i, j, k, icell, nx, ny, nz, i0, j0, k0, lmax, nlevels, level, itmp
  real(kind=RealKind) :: x0, y0, z0, xa, xb, ya, yb, za, zb, tmp1, tmp2, &
       xnew, ynew, znew, xpos, ypos, zpos
  type(readLevelType), dimension(:), pointer :: readLevel
  character(90) :: dirname, filein, fileout

  dirname = './lyalpha/'
  filein = 'S29COSMO-055-8.2200kpc_128_z5.76_velmet.h4'
  fileout = 'cellArray.dat'

  ! --------------- read binary format

  open(14,file=trim(dirname)//trim(filein), status='old',form='unformatted')

  read(14) nlevels
  print*, 'nlevels =', nlevels

  allocate(readLevel(nlevels))

  do level = 1, nlevels

     read(14) readLevel(level)%ncell
     write(*,*) 'level =', level, readLevel(level)%ncell

     allocate(readLevel(level)%pos(readLevel(level)%ncell,3))
     allocate(readLevel(level)%lT(readLevel(level)%ncell))
     allocate(readLevel(level)%lnH(readLevel(level)%ncell))
     allocate(readLevel(level)%lx(readLevel(level)%ncell))
!     allocate(readLevel(level)%abun(readLevel(level)%ncell,4))
!     allocate(readLevel(level)%vel(readLevel(level)%ncell,3))

     read(14) (readLevel(level)%pos(icell,1),icell=1,readLevel(level)%ncell)
     read(14) (readLevel(level)%pos(icell,2),icell=1,readLevel(level)%ncell)
     read(14) (readLevel(level)%pos(icell,3),icell=1,readLevel(level)%ncell)
     read(14) (readLevel(level)%lT(icell),icell=1,readLevel(level)%ncell)
     read(14) (readLevel(level)%lnH(icell),icell=1,readLevel(level)%ncell)
     read(14) (readLevel(level)%lx(icell),icell=1,readLevel(level)%ncell)
!     read(14) (readLevel(level)%abun(icell,1),icell=1,readLevel(level)%ncell)
!     read(14) (readLevel(level)%abun(icell,2),icell=1,readLevel(level)%ncell)
!     read(14) (readLevel(level)%abun(icell,3),icell=1,readLevel(level)%ncell)
!     read(14) (readLevel(level)%abun(icell,4),icell=1,readLevel(level)%ncell)
!     read(14) (readLevel(level)%vel(icell,1),icell=1,readLevel(level)%ncell)
!     read(14) (readLevel(level)%vel(icell,2),icell=1,readLevel(level)%ncell)
!     read(14) (readLevel(level)%vel(icell,3),icell=1,readLevel(level)%ncell)

     if (level .eq. 2) then
        print*, readLevel(level)%pos(1,:)
     endif

  enddo

  close(14)

  ! --------------- set up a fully nested 3D grid

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
           baseGrid%cell(i,j,k)%level = 0
           baseGrid%cell(i,j,k)%parent => baseGrid
           nullify(baseGrid%cell(i,j,k)%cell)
        enddo
     enddo
  enddo

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

        call placeCellProject(baseGrid%cell(i0,j0,k0),level,xnew,ynew,znew, &
             readLevel(level)%lT(icell),readLevel(level)%lnH(icell),readLevel(level)%lx(icell))

     enddo
  enddo

  do level = 1, nlevels
     deallocate(readLevel(level)%pos)
     deallocate(readLevel(level)%lT)
     deallocate(readLevel(level)%lnH)
     deallocate(readLevel(level)%lx)
!     deallocate(readLevel(level)%vel)
  enddo
  deallocate(readLevel)

  ! --------------- store hierarchy with space-filling curve in binary format

  icosmic = 0
  do i = 1, nx
     do j = 1, ny
        do k = 1, nz
           call countCells(baseGrid%cell(i,j,k))
        enddo
     enddo
  enddo
  ncosmic = icosmic
  write(*,*) 'total number of cells =', ncosmic

  allocate (cellArrayLevel(ncosmic))
  allocate (cellArrayXpos(ncosmic))
  allocate (cellArrayYpos(ncosmic))
  allocate (cellArrayZpos(ncosmic))
  allocate (cellArrayHI(ncosmic))
  allocate (cellArrayHeI(ncosmic))
  allocate (cellArrayHeII(ncosmic))
  allocate (cellArrayTemp(ncosmic))
  allocate (cellArrayDensity(ncosmic))

  icosmic = 0
  do i = 1, nx
     do j = 1, ny
        do k = 1, nz
           call writeCell(baseGrid%cell(i,j,k),0)
        enddo
     enddo
  enddo

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

  itmp = 0
  do icosmic = 1, ncosmic
     cellArrayXpos(icosmic) = cellArrayXpos(icosmic)*(xb-xa) + xa
     cellArrayYpos(icosmic) = cellArrayYpos(icosmic)*(yb-ya) + ya
     cellArrayZpos(icosmic) = cellArrayZpos(icosmic)*(zb-za) + za
     itmp = max(itmp,cellArrayLevel(icosmic))
  enddo
  print*, 'base grid level =', cellArrayLevel(1)
  print*, 'max level of refinement =', itmp

  open(14,file=trim(dirname)//trim(fileout),status='replace',form='unformatted')
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

  do icosmic = 1, ncosmic
     if (cellArrayLevel(icosmic).eq.1) then
        print*, icosmic, cellArrayXpos(icosmic), cellArrayYpos(icosmic), cellArrayZpos(icosmic)
        stop
     endif
  enddo

contains

  recursive subroutine placeCellProject(parentCell,level,x0,y0,z0,logtgas,lognh,logxneu)

    use localDefinitions

    implicit none
    integer, intent(in) :: level
    integer :: inew, jnew, knew, i, j, k
    real(kind=RealKind), intent(in) :: x0, y0, z0
    real*4, intent(in) :: logtgas, lognh, logxneu
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
                   parentCell%cell(i,j,k)%parent => parentCell
                   nullify(parentCell%cell(i,j,k)%cell)
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
       ! nh = 1.e-5 ! uniform for testing
       xneu = 10.**logxneu
       parentCell%rho = nh * mh/psi
       parentCell%HI = nh * xneu
       nhe = (1.-psi) * parentCell%rho / mhe
       parentCell%HeI = nhe * 1.
       parentCell%HeII = nhe * 0.

       ! print*, 'rho,temp =', parentCell%rho, parentCell%tgas
       ! print*, 'HI,xneu =', parentCell%HI, xneu
       ! print*, 'HeI,HeII =', parentCell%HeI, parentCell%HeII

    endif

  end subroutine placeCellProject

  recursive subroutine countCells(currentCell)

    use localDefinitions

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
    endif

  end subroutine countCells

  recursive subroutine writeCell(currentCell,level)

    use localDefinitions

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
    endif

  end subroutine writeCell

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
                   call computeCellCoordinates(xnew,ynew,znew,level+1,cellSize/2.)
                enddo
             enddo
          enddo
       else
          write(*,*) 'error in levels', icosmic, cellArrayLevel(icosmic), level
          stop
       endif
    endif

  end subroutine computeCellCoordinates

end program convertFormats
