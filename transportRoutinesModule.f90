module transportRoutinesModule

  implicit none

contains

  subroutine setPattern(pattern,phi,theta)

    use definitions

    implicit none
    real(kind=RealKind), intent(in) :: phi, theta
    real(kind=RealKind) :: tmp1, tmp2, tmp3, tmpa1, tmpa2, tmpb1, tmpb2
    type(patternType) :: pattern

    tmp1 = 1./sin(theta)
    tmp2 = (1.-pattern%xyRay%x0)/(cos(phi)*cos(theta))
    tmp3 = (1.-pattern%xyRay%y0)/(sin(phi)*cos(theta))

    if (tmp1.lt.min(tmp2,tmp3)) then
       pattern%xyRay%len = tmp1
       pattern%xzRayActive = .false.
       pattern%yzRayActive = .false.
       pattern%xyTop = xyEnd
       pattern%xzTop = 0
       pattern%yzTop = 0
    else
       if (tmp2.lt.min(tmp1,tmp3)) then
          pattern%xyRay%len = tmp2
          pattern%yzRayActive = .true.
          pattern%yzRay%y0 = (1.-pattern%xyRay%x0)*tan(phi) + pattern%xyRay%y0
          pattern%yzRay%z0 = pattern%xyRay%len*sin(theta)
          if (pattern%yzRay%y0.gt.1. .or. pattern%yzRay%z0.gt.1.) then
             write(*,*) 'Error in yzRay%y0, z0', pattern%yzRay%y0, pattern%yzRay%z0
             stop
          endif
          tmpa1 = (1.-pattern%yzRay%z0)/sin(theta)
          tmpa2 = (1.-pattern%yzRay%y0)/(sin(phi)*cos(theta))
          if (tmpa1.lt.tmpa2) then
             pattern%yzRay%len = tmpa1
             pattern%xzRayActive = .false.
             pattern%xyTop = yzEnd
             pattern%xzTop = 0
             pattern%yzTop = xyEnd
          else
             pattern%yzRay%len = tmpa2
             pattern%xzRayActive = .true.
             pattern%xzRay%x0 = (1.-pattern%yzRay%y0)/tan(phi)
             pattern%xzRay%z0 = pattern%yzRay%z0 + tmpa2*sin(theta)
             pattern%xzRay%len = (1.-pattern%xzRay%z0)/sin(theta)
             pattern%xyTop = xzEnd
             pattern%xzTop = yzEnd
             pattern%yzTop = xyEnd
          endif
       else
          pattern%xyRay%len = tmp3
          pattern%xzRayActive = .true.
          pattern%xzRay%x0 = (1.-pattern%xyRay%y0)/tan(phi) + pattern%xyRay%x0
          pattern%xzRay%z0 = tmp3*sin(theta)
          if (pattern%xzRay%x0.gt.1. .or. pattern%xzRay%z0.gt.1.) then
             write(*,*) 'Error in xzRay%x0,z0', pattern%xzRay%x0, pattern%xzRay%z0
             stop
          endif
          tmpb1 = (1.-pattern%xzRay%z0)/sin(theta)
          tmpb2 = (1.-pattern%xzRay%x0)/(cos(phi)*cos(theta))
          if (tmpb1.lt.tmpb2) then
             pattern%xzRay%len = tmpb1
             pattern%yzRayActive = .false.
             pattern%xyTop = xzEnd
             pattern%xzTop = xyEnd
             pattern%yzTop = 0
          else
             pattern%xzRay%len = tmpb2
             pattern%yzRayActive = .true.
             pattern%yzRay%y0 = (1.-pattern%xzRay%x0)*tan(phi)
             pattern%yzRay%z0 = pattern%xzRay%len*sin(theta) + pattern%xzRay%z0
             pattern%yzRay%len = (1.-pattern%yzRay%z0)/sin(theta)
             pattern%xyTop = yzEnd
             pattern%xzTop = xyEnd
             pattern%yzTop = xzEnd
          endif
       endif
    endif

  end subroutine setPattern

  subroutine copyPatternToCell(pattern,cell)

    use definitions

    implicit none
    type(patternType), target :: pattern
    type(zoneType) :: cell

!  cell%rt%xyRay%x0 = pattern%xyRay%x0
!  cell%rt%xyRay%y0 = pattern%xyRay%y0
!  cell%rt%xyRay%len = pattern%xyRay%len

!  cell%rt%xzRayActive = pattern%xzRayActive
!  if (pattern%xzRayActive) then
!     cell%rt%xzRay%x0 = pattern%xzRay%x0
!     cell%rt%xzRay%z0 = pattern%xzRay%z0
!     cell%rt%xzRay%len = pattern%xzRay%len
!  endif

!  cell%rt%yzRayActive = pattern%yzRayActive
!  if (pattern%yzRayActive) then
!     cell%rt%yzRay%y0 = pattern%yzRay%y0
!     cell%rt%yzRay%z0 = pattern%yzRay%z0
!     cell%rt%yzRay%len = pattern%yzRay%len
!  endif

!  cell%xyTop = pattern%xyTop
!  cell%xzTop = pattern%xzTop
!  cell%yzTop = pattern%yzTop

    cell%pattern => pattern

  end subroutine copyPatternToCell

  recursive subroutine setRaysRefined(parentCell,parentPattern, &
       is,js,ks,phi,theta)

    use definitions

    implicit none
    type(zoneType) :: parentCell
    type(patternType) :: parentPattern
    integer, dimension(2,2,2), intent(in) :: is, js, ks
    real(kind=RealKind), intent(in) :: phi, theta
    integer :: i, j, k
    type(zoneType), pointer :: currentCell
    type(patternType), pointer :: patternBelow, currentPattern

    if (.not.parentPattern%refined) then

       allocate(parentPattern%cell(2,2,2))
       parentPattern%refined = .true.

       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                currentPattern => parentPattern%cell(is(i,j,k),js(i,j,k),ks(i,j,k))
                currentPattern%refined = .false.
                nullify(currentPattern%cell)
             enddo
          enddo
       enddo

       currentPattern => parentPattern%cell(is(1,1,1),js(1,1,1),ks(1,1,1))
       if (parentPattern%xyRay%x0.lt.0.5) then
          currentPattern%xyRay%x0 = 2.*parentPattern%xyRay%x0
       else
          currentPattern%xyRay%x0 = 2.*parentPattern%xyRay%x0 - 1.
       endif
       if (parentPattern%xyRay%y0.lt.0.5) then
          currentPattern%xyRay%y0 = 2.*parentPattern%xyRay%y0
       else
          currentPattern%xyRay%y0 = 2.*parentPattern%xyRay%y0 - 1.
       endif
       call setPattern(currentPattern,phi,theta)

!     call checkPattern(currentPattern,phi,theta)

       patternBelow => parentPattern%cell(is(1,1,1),js(1,1,1),ks(1,1,1))
       currentPattern => parentPattern%cell(is(2,1,1),js(2,1,1),ks(2,1,1))
       select case (patternBelow%xyTop)
       case(xyEnd) ! use patternBelow%xyRay
          currentPattern%xyRay%x0 = patternBelow%xyRay%x0 + cos(phi)/tan(theta)
          currentPattern%xyRay%y0 = patternBelow%xyRay%y0 + sin(phi)/tan(theta)
       case(xzEnd) ! use patternBelow%xzRay
          currentPattern%xyRay%x0 = patternBelow%xzRay%x0 + &
               patternBelow%xzRay%len*cos(theta)*cos(phi)
          currentPattern%xyRay%y0 = patternBelow%xzRay%len*cos(theta)*sin(phi)
       case(yzEnd) ! use patternBelow%yzRay
          currentPattern%xyRay%x0 = patternBelow%yzRay%len*cos(theta)*cos(phi)
          currentPattern%xyRay%y0 = patternBelow%yzRay%y0 + &
               patternBelow%yzRay%len*cos(theta)*sin(phi)
       case(0)
          write(*,*) 'error in xyTop'
          stop
       end select
       if (currentPattern%xyRay%x0.gt.1. .or. currentPattern%xyRay%y0.gt.1.) then
          write(*,*) '2) Error: xyRay%x0, y0', currentPattern%xyRay%x0, currentPattern%xyRay%y0
          stop
       endif
       call setPattern(currentPattern,phi,theta)

!     call checkPattern(currentPattern,phi,theta)

       do i = 1, 2
          currentPattern => parentPattern%cell(is(i,1,1),js(i,1,1),ks(i,1,1))
          parentPattern%cell(is(i,1,2),js(i,1,2),ks(i,1,2)) = currentPattern
          parentPattern%cell(is(i,2,1),js(i,2,1),ks(i,2,1)) = currentPattern
          parentPattern%cell(is(i,2,2),js(i,2,2),ks(i,2,2)) = currentPattern
       enddo

    endif

    do i = 1, 2
       do j = 1, 2
          do k = 1, 2
             currentPattern => parentPattern%cell(is(i,j,k),js(i,j,k),ks(i,j,k))
             currentCell => parentCell%cell(is(i,j,k),js(i,j,k),ks(i,j,k))
!           call copyPatternToCell(currentPattern,currentCell)
             currentCell%pattern => currentPattern
             if (currentCell%refined) call setRaysRefined(currentCell,currentPattern, &
                  is,js,ks,phi,theta)
          enddo
       enddo
    enddo

!  write(*,*) parentPattern%cell(is(1,1,1),js(1,1,1),ks(1,1,1))%xyRay%x0, &
!       parentPattern%cell(is(1,1,1),js(1,1,1),ks(1,1,1))%xyRay%y0
!  write(*,*) parentPattern%cell(is(2,1,1),js(2,1,1),ks(2,1,1))%xyRay%x0, &
!       parentPattern%cell(is(2,1,1),js(2,1,1),ks(2,1,1))%xyRay%y0

  end subroutine setRaysRefined

  subroutine checkPattern(pattern,phi,theta)

    use definitions

    implicit none
    type(patternType), intent(in) :: pattern
    real(kind=RealKind), intent(in) :: phi, theta
    real(kind=RealKind) :: totalRayLength, upperx0, uppery0

    select case (pattern%xyTop)
    case(xyEnd) ! use pattern%xyRay
       upperx0 = pattern%xyRay%x0 + cos(phi)/tan(theta)
       uppery0 = pattern%xyRay%y0 + sin(phi)/tan(theta)
    case(xzEnd) ! use pattern%xzRay
       upperx0 = pattern%xzRay%x0 + pattern%xzRay%len*cos(theta)*cos(phi)
       uppery0 = pattern%xzRay%len*cos(theta)*sin(phi)
    case(yzEnd) ! use pattern%yzRay
       upperx0 = pattern%yzRay%len*cos(theta)*cos(phi)
       uppery0 = pattern%yzRay%y0 + pattern%yzRay%len*cos(theta)*sin(phi)
    case(0)
       write(*,*) 'error in xyTop'
       stop
    end select

    write(*,'("xyRay: ",2f15.10)') pattern%xyRay%x0, pattern%xyRay%y0
    totalRayLength = pattern%xyRay%len
    if (pattern%xzRayActive) then
       write(*,'("xzRay: ",l1,2f15.10)') pattern%xzRayActive, pattern%xzRay%x0, pattern%xzRay%z0
       totalRayLength = totalRayLength + pattern%xzRay%len
    else
       write(*,'("xzRay: ",l1)') pattern%xzRayActive
    endif
    if (pattern%yzRayActive) then
       write(*,'("yzRay: ",l1,2f15.10)') pattern%yzRayActive, pattern%yzRay%y0, pattern%yzRay%z0
       totalRayLength = totalRayLength + pattern%yzRay%len
    else
       write(*,'("yzRay: ",l1)') pattern%yzRayActive
    endif
    write(*,*) 'len =', totalRayLength, 1./sin(theta)
    write(*,*) 'xyTop =', pattern%xyTop, upperx0, uppery0
    write(*,*) '----------------------------'

  end subroutine checkPattern

  recursive subroutine findNeighbours(finestLevelCell,level,callSequence,is,js,ks,nx,ny,nz,izone)

    use definitions
    use rotateIndicesModule

    implicit none
    type(zoneType), target :: finestLevelCell
    integer, intent(in) :: level
    integer, intent(in) :: callSequence(3*level+3)
    integer, dimension(2,2,2), intent(in) :: is, js, ks
    integer, intent(in) :: nx, ny, nz
    integer(1), intent(in) :: izone
    integer :: ilevel, i, j, k, icell, jcell, kcell
    real(kind=RealKind) :: xCoarserUnits, yCoarserUnits, zCoarserUnits
    type(zoneType), pointer :: coarserCell

! find xy-neighbour

!   write(*,*) 'looking for xy-neighbour at level', level
!   write(*,*) '>>>', finestLevelCell%chi
!   write(*,*) 'callSequence =', callSequence
!   stop

    xCoarserUnits = finestLevelCell%pattern%xyRay%x0
    yCoarserUnits = finestLevelCell%pattern%xyRay%y0
    coarserCell => finestLevelCell
    do ilevel = level, 0, -1
       i = callSequence(3*ilevel+1)
       j = callSequence(3*ilevel+2)
       k = callSequence(3*ilevel+3)
       coarserCell => coarserCell%parent
       if (i.gt.1) then
          if (ilevel.eq.0) then
             call rotateIndices(i-1,j,k,nx,ny,nz,izone,icell,jcell,kcell)
          else
             icell = is(i-1,j,k)
             jcell = js(i-1,j,k)
             kcell = ks(i-1,j,k)
          endif
          call getXYNeighbour(finestLevelCell,coarserCell%cell(icell,jcell,kcell), &
               xCoarserUnits,yCoarserUnits,is,js,ks)
          goto 2100
       else
          if (j.eq.1) then
             yCoarserUnits = yCoarserUnits/2.
          else
             yCoarserUnits = yCoarserUnits/2. + 0.5
          endif
          if (k.eq.1) then
             xCoarserUnits = xCoarserUnits/2.
          else
             xCoarserUnits = xCoarserUnits/2. + 0.5
          endif
       endif
    enddo
    finestLevelCell%xyNeighbourPresent = .false.
!  write(*,*) 'must be boundary cell: no xy-neighbour found'
2100 continue

!   if (finestLevelCell%xyNeighbourPresent) then
!      write(*,*) 'xy: ', finestLevelCell%xyNeighbour%refined, finestLevelCell%xyNeighbour%chi
!      stop
!   endif

! find xz-neighbour

    if (finestLevelCell%pattern%xzRayActive) then
       xCoarserUnits = finestLevelCell%pattern%xzRay%x0
       zCoarserUnits = finestLevelCell%pattern%xzRay%z0
       coarserCell => finestLevelCell
       do ilevel = level, 0, -1
          i = callSequence(3*ilevel+1)
          j = callSequence(3*ilevel+2)
          k = callSequence(3*ilevel+3)
          coarserCell => coarserCell%parent
          if (j.gt.1) then
             if (ilevel.eq.0) then
                call rotateIndices(i,j-1,k,nx,ny,nz,izone,icell,jcell,kcell)
             else
                icell = is(i,j-1,k)
                jcell = js(i,j-1,k)
                kcell = ks(i,j-1,k)
             endif
             call getXZNeighbour(finestLevelCell,coarserCell%cell(icell,jcell,kcell), &
                  xCoarserUnits,zCoarserUnits,is,js,ks)
             goto 2200
          else
             if (i.eq.1) then
                zCoarserUnits = zCoarserUnits/2.
             else
                zCoarserUnits = zCoarserUnits/2. + 0.5
             endif
             if (k.eq.1) then
                xCoarserUnits = xCoarserUnits/2.
             else
                xCoarserUnits = xCoarserUnits/2. + 0.5
             endif
          endif
       enddo
       finestLevelCell%xzNeighbourPresent = .false.
!     write(*,*) 'must be boundary cell: no xz-neighbour found'
2200   continue
    endif

!   if (finestLevelCell%xzNeighbourPresent) then
!      write(*,*) 'xz: ', finestLevelCell%xzNeighbour%refined, finestLevelCell%xzNeighbour%chi
!      stop
!   endif

! find yz-neighbour

    if (finestLevelCell%pattern%yzRayActive) then
       yCoarserUnits = finestLevelCell%pattern%yzRay%y0
       zCoarserUnits = finestLevelCell%pattern%yzRay%z0
       coarserCell => finestLevelCell
       do ilevel = level, 0, -1
          i = callSequence(3*ilevel+1)
          j = callSequence(3*ilevel+2)
          k = callSequence(3*ilevel+3)
          coarserCell => coarserCell%parent
          if (k.gt.1) then
             if (ilevel.eq.0) then
                call rotateIndices(i,j,k-1,nx,ny,nz,izone,icell,jcell,kcell)
             else
                icell = is(i,j,k-1)
                jcell = js(i,j,k-1)
                kcell = ks(i,j,k-1)
             endif
             call getYZNeighbour(finestLevelCell,coarserCell%cell(icell,jcell,kcell), &
                  yCoarserUnits,zCoarserUnits,is,js,ks)
             goto 2300
          else
             if (i.eq.1) then
                zCoarserUnits = zCoarserUnits/2.
             else
                zCoarserUnits = zCoarserUnits/2. + 0.5
             endif
             if (j.eq.1) then
                yCoarserUnits = yCoarserUnits/2.
             else
                yCoarserUnits = yCoarserUnits/2. + 0.5
             endif
          endif
       enddo
       finestLevelCell%yzNeighbourPresent = .false.
!     write(*,*) 'must be boundary cell: no yz-neighbour found'
2300   continue
    endif

!   if (finestLevelCell%yzNeighbourPresent) then
!      write(*,*) 'yz: ', finestLevelCell%yzNeighbour%refined, finestLevelCell%yzNeighbour%chi
!      stop
!   endif

  end subroutine findNeighbours


  recursive subroutine localizeCellFindNeighbours(parentCell,level,callSequence,is,js,ks,nx,ny,nz,izone)

    use definitions

    implicit none
    type(zoneType) :: parentCell
    integer, intent(in) :: level
    integer, intent(in) :: callSequence(3*level+3)
    integer, dimension(2,2,2), intent(in) :: is, js, ks
    integer, intent(in) :: nx, ny, nz
    integer(1), intent(in) :: izone
    integer :: i, j, k, callLength
    integer :: newCallSequence(3*level+6)

    callLength = 3*level + 3

    if (parentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
!              write(*,'("entered cell",3i4,"   of level",i5)') is(i,j,k), js(i,j,k), ks(i,j,k), level+1
                newCallSequence(1:callLength) = callSequence
                newCallSequence(callLength+1:callLength+3) = (/i,j,k/)
                call localizeCellFindNeighbours(parentCell%cell(is(i,j,k),js(i,j,k),ks(i,j,k)), &
                     level+1,newCallSequence,is,js,ks,nx,ny,nz,izone)
             enddo
          enddo
       enddo
    else
       call findNeighbours(parentCell,level,callSequence,is,js,ks,nx,ny,nz,izone)
    endif

  end subroutine localizeCellFindNeighbours

  recursive subroutine getXYNeighbour(cellToGetXYNeighbour,neighbourContainerCell,x0,y0,is,js,ks)

    use definitions

    implicit none
    type(zoneType), target :: cellToGetXYNeighbour
    type(zoneType), target :: neighbourContainerCell
    real(kind=RealKind), intent(in) :: x0, y0
    integer, dimension(2,2,2), intent(in) :: is, js, ks

    if (neighbourContainerCell%refined) then
       if (x0.le.0.5) then
          if (y0.le.0.5) then
             call getXYNeighbour(cellToGetXYNeighbour, &
                  neighbourContainerCell%cell(is(2,1,1),js(2,1,1),ks(2,1,1)),2.*x0,2.*y0,is,js,ks)
          else
             call getXYNeighbour(cellToGetXYNeighbour, &
                  neighbourContainerCell%cell(is(2,2,1),js(2,2,1),ks(2,2,1)),2.*x0,2.*y0-1.,is,js,ks)
          endif
       else
          if (y0.le.0.5) then
             call getXYNeighbour(cellToGetXYNeighbour, &
                  neighbourContainerCell%cell(is(2,1,2),js(2,1,2),ks(2,1,2)),2.*x0-1.,2.*y0,is,js,ks)
          else
             call getXYNeighbour(cellToGetXYNeighbour, &
                  neighbourContainerCell%cell(is(2,2,2),js(2,2,2),ks(2,2,2)),2.*x0-1.,2.*y0-1.,is,js,ks)
          endif
       endif
    else
       cellToGetXYNeighbour%xyNeighbourPresent = .true.
       cellToGetXYNeighbour%xyNeighbour => neighbourContainerCell
    endif

  end subroutine getXYNeighbour

  recursive subroutine getXZNeighbour(cellToGetXZneighbour,neighbourContainerCell,x0,z0,is,js,ks)

    use definitions

    implicit none
    type(zoneType), target :: cellToGetXZneighbour
    type(zoneType), target :: neighbourContainerCell
    real(kind=RealKind), intent(in) :: x0, z0
    integer, dimension(2,2,2), intent(in) :: is, js, ks

    if (neighbourContainerCell%refined) then
       if (x0.le.0.5) then
          if (z0.le.0.5) then
             call getXZNeighbour(cellToGetXZNeighbour, &
                  neighbourContainerCell%cell(is(1,2,1),js(1,2,1),ks(1,2,1)),2.*x0,2.*z0,is,js,ks)
          else
             call getXZNeighbour(cellToGetXZNeighbour, &
                  neighbourContainerCell%cell(is(2,2,1),js(2,2,1),ks(2,2,1)),2.*x0,2.*z0-1.,is,js,ks)
          endif
       else
          if (z0.le.0.5) then
             call getXZNeighbour(cellToGetXZNeighbour, &
                  neighbourContainerCell%cell(is(1,2,2),js(1,2,2),ks(1,2,2)),2.*x0-1.,2.*z0,is,js,ks)
          else
             call getXZNeighbour(cellToGetXZNeighbour, &
                  neighbourContainerCell%cell(is(2,2,2),js(2,2,2),ks(2,2,2)),2.*x0-1.,2.*z0-1.,is,js,ks)
          endif
       endif
    else
       cellToGetXZneighbour%xzNeighbourPresent = .true.
       cellToGetXZneighbour%xzNeighbour => neighbourContainerCell
    endif

  end subroutine getXZNeighbour

  recursive subroutine getYZNeighbour(cellToGetYZneighbour,neighbourContainerCell,y0,z0,is,js,ks)

    use definitions

    implicit none
    type(zoneType), target :: cellToGetYZneighbour
    type(zoneType), target :: neighbourContainerCell
    real(kind=RealKind), intent(in) :: y0, z0
    integer, dimension(2,2,2), intent(in) :: is, js, ks

    if (neighbourContainerCell%refined) then
       if (y0.le.0.5) then
          if (z0.le.0.5) then
             call getYZNeighbour(cellToGetYZNeighbour, &
                  neighbourContainerCell%cell(is(1,1,2),js(1,1,2),ks(1,1,2)),2.*y0,2.*z0,is,js,ks)
          else
             call getYZNeighbour(cellToGetYZNeighbour, &
                  neighbourContainerCell%cell(is(2,1,2),js(2,1,2),ks(2,1,2)),2.*y0,2.*z0-1.,is,js,ks)
          endif
       else
          if (z0.le.0.5) then
             call getYZNeighbour(cellToGetYZNeighbour, &
                  neighbourContainerCell%cell(is(1,2,2),js(1,2,2),ks(1,2,2)),2.*y0-1.,2.*z0,is,js,ks)
          else
             call getYZNeighbour(cellToGetYZNeighbour, &
                  neighbourContainerCell%cell(is(2,2,2),js(2,2,2),ks(2,2,2)),2.*y0-1.,2.*z0-1.,is,js,ks)
          endif
       endif
    else
       cellToGetYZneighbour%yzNeighbourPresent = .true.
       cellToGetYZneighbour%yzNeighbour => neighbourContainerCell
    endif

  end subroutine getYZNeighbour

  recursive subroutine transport(transportCell,weight,is,js,ks,cellSizeAbsoluteUnits)

    use definitions

    implicit none
    type(zoneType) :: transportCell
    real(kind=RealKind), intent(in) :: weight
    integer, dimension(2,2,2), intent(in) :: is, js, ks
    real(kind=RealKind), intent(in) :: cellSizeAbsoluteUnits
    real(kind=RealKind) :: Iin1, Iin2, Iin3, Itmp1, Itmp2, Itmp3, &
         Jmean1, Jmean2, Jmean3, &
         dpath, tau1, tau2, tau3, tmpabs1, tmpabs2, tmpabs3, &
         tmpemi1, tmpemi2, tmpemi3, nemi1, nemi2, nemi3, tmp
    integer :: i, j, k, imean
    type(patternType), pointer :: neighbourPattern
    type(radiativeTransportType), pointer :: neighbourTransport

    if (transportCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
!              write(*,'("entered cell",3i4)') is(i,j,k), js(i,j,k), ks(i,j,k)
                call transport(transportCell%cell(is(i,j,k),js(i,j,k),ks(i,j,k)), &
                     weight,is,js,ks,cellSizeAbsoluteUnits/2.)
             enddo
          enddo
       enddo
    else
       Jmean1 = 0.
       Jmean2 = 0.
       Jmean3 = 0.
       imean = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!! xy-ray transport !!!!!!!!!!!!!!!!!!!!!!!!!!!
!     write(*,*) 'xy-ray transport'
       if (.not.transportCell%xyNeighbourPresent) then
          Iin1 = uvb1
          Iin2 = uvb2
          Iin3 = uvb3
       else
          select case (transportCell%xyNeighbour%pattern%xyTop)
          case(xyEnd)
             Itmp1 = transportCell%xyNeighbour%rt%xyRay%Iout1
             Itmp2 = transportCell%xyNeighbour%rt%xyRay%Iout2
             Itmp3 = transportCell%xyNeighbour%rt%xyRay%Iout3
          case(xzEnd)
             Itmp1 = transportCell%xyNeighbour%rt%xzRay%Iout1
             Itmp2 = transportCell%xyNeighbour%rt%xzRay%Iout2
             Itmp3 = transportCell%xyNeighbour%rt%xzRay%Iout3
          case(yzEnd)
             Itmp1 = transportCell%xyNeighbour%rt%yzRay%Iout1
             Itmp2 = transportCell%xyNeighbour%rt%yzRay%Iout2
             Itmp3 = transportCell%xyNeighbour%rt%yzRay%Iout3
          case(0)
             if (transportCell%level.le.transportCell%xyNeighbour%level) then
                write(*,*) 'error in xyTop'
                stop
             endif

             neighbourPattern => transportCell%xyNeighbour%pattern
             neighbourTransport => transportCell%xyNeighbour%rt
             if (neighbourPattern%xzRayActive) then
                Itmp1 = 0.5 * (neighbourTransport%xzRay%Iout1+neighbourTransport%xyRay%Iout1)
                Itmp2 = 0.5 * (neighbourTransport%xzRay%Iout2+neighbourTransport%xyRay%Iout2)
                Itmp3 = 0.5 * (neighbourTransport%xzRay%Iout3+neighbourTransport%xyRay%Iout3)
             else
                if (neighbourPattern%yzRayActive) then
                   Itmp1 = 0.5 * (neighbourTransport%yzRay%Iout1+neighbourTransport%xyRay%Iout1)
                   Itmp2 = 0.5 * (neighbourTransport%yzRay%Iout2+neighbourTransport%xyRay%Iout2)
                   Itmp3 = 0.5 * (neighbourTransport%yzRay%Iout3+neighbourTransport%xyRay%Iout3)
                else
                   Itmp1 = neighbourTransport%xyRay%Iout1
                   Itmp2 = neighbourTransport%xyRay%Iout2
                   Itmp3 = neighbourTransport%xyRay%Iout3
                endif
             endif

          end select
          if (transportCell%level.le.transportCell%xyNeighbour%level) then
             Iin1 = Itmp1
             Iin2 = Itmp2
             Iin3 = Itmp3
          else
             Iin1 = Itmp1
             Iin2 = Itmp2
             Iin3 = Itmp3
!         write(*,*) 'xy-interpolation not implemented yet', &
!              transportCell%pattern%xyRay%x0, transportCell%pattern%xyRay%y0
!         stop
          endif
       endif

       dpath = cellSizeAbsoluteUnits * transportCell%pattern%xyRay%len ! [cm]
       tau1 = transportCell%kappa1 * dpath
       tau2 = transportCell%kappa2 * dpath
       tau3 = transportCell%kappa3 * dpath
       tmpabs1 = exp(-tau1)
       tmpabs2 = exp(-tau2)
       tmpabs3 = exp(-tau3)
       if (tau1.gt.1.e-10) then
          tmpemi1 = (1.-tmpabs1)/transportCell%kappa1
       else
          tmpemi1 = dpath
       endif
       if (tau2.gt.1.e-10) then
          tmpemi2 = (1.-tmpabs2)/transportCell%kappa2
       else
          tmpemi2 = dpath
       endif
       if (tau3.gt.1.e-10) then
          tmpemi3 = (1.-tmpabs3)/transportCell%kappa3
       else
          tmpemi3 = dpath
       endif
       nemi1 = 0. ! emissivity [erg/s/cm^2/Hz/ster]
       nemi2 = 0.
       nemi3 = 0.
       transportCell%rt%xyRay%Iout1 = Iin1*tmpabs1 + nemi1*tmpemi1/dpath
       transportCell%rt%xyRay%Iout2 = Iin2*tmpabs2 + nemi2*tmpemi2/dpath
       transportCell%rt%xyRay%Iout3 = Iin3*tmpabs3 + nemi3*tmpemi3/dpath

       tmp = transportCell%rt%xyRay%Iout1+transportCell%rt%xyRay%Iout2+transportCell%rt%xyRay%Iout3
       if (tmp.lt.1.d-20 .and. tmp.gt.-1.d-20) then
          continue
       else
          print*, tmpabs1, tau1
          print*, transportCell%kappa1, dpath
          print*, transportCell%HI, transportCell%HeI, transportCell%HeII
          stop ! abc
       endif

!     Jmean1 = Jmean1 + Iin1 + transportCell%rt%xyRay%Iout1
!     Jmean2 = Jmean2 + Iin2 + transportCell%rt%xyRay%Iout2
!     Jmean3 = Jmean3 + Iin3 + transportCell%rt%xyRay%Iout3
!     imean = imean + 2

       call computeCellIntensity(Jmean1,Iin1,transportCell%rt%xyRay%Iout1)
       call computeCellIntensity(Jmean2,Iin2,transportCell%rt%xyRay%Iout2)
       call computeCellIntensity(Jmean3,Iin3,transportCell%rt%xyRay%Iout3)
       imean = imean + 1

!      if (transportCell%rt%xyRay%Iout1.eq.0.) then
!         write(*,*) 'a) zero Iout1 ...'
!         write(*,*) Iin1, transportCell%xyNeighbourPresent, transportCell%xyNeighbour%pattern%xyTop
!         stop
!      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!! xz-ray transport !!!!!!!!!!!!!!!!!!!!!!!!!!!
!     write(*,*) 'xz-ray transport', transportCell%rt%xzRayActive
       if (transportCell%pattern%xzRayActive) then
          if (.not.transportCell%xzNeighbourPresent) then
             Iin1 = uvb1
             Iin2 = uvb2
             Iin3 = uvb3
          else
             select case (transportCell%xzNeighbour%pattern%xzTop)
             case(xyEnd)
                Itmp1 = transportCell%xzNeighbour%rt%xyRay%Iout1
                Itmp2 = transportCell%xzNeighbour%rt%xyRay%Iout2
                Itmp3 = transportCell%xzNeighbour%rt%xyRay%Iout3
             case(xzEnd)
                if (.not.transportCell%xzNeighbour%pattern%xzRayActive) then
                   write(*,*) '3) Error: xzRay should be active ...', i, j, k
                   stop
                endif
                Itmp1 = transportCell%xzNeighbour%rt%xzRay%Iout1
                Itmp2 = transportCell%xzNeighbour%rt%xzRay%Iout2
                Itmp3 = transportCell%xzNeighbour%rt%xzRay%Iout3
             case(yzEnd)
                if (.not.transportCell%xzNeighbour%pattern%yzRayActive) then
                   write(*,*) '3) Error: yzRay should be active ...', i, j, k
                   stop
                endif
                Itmp1 = transportCell%xzNeighbour%rt%yzRay%Iout1
                Itmp2 = transportCell%xzNeighbour%rt%yzRay%Iout2
                Itmp3 = transportCell%xzNeighbour%rt%yzRay%Iout3
             case(0)
                if (transportCell%level.le.transportCell%xzNeighbour%level) then
                   write(*,*) 'error in xzTop'
                   stop
                endif

                neighbourPattern => transportCell%xzNeighbour%pattern
                neighbourTransport => transportCell%xzNeighbour%rt
                if (neighbourPattern%xzRayActive) then
                   Itmp1 = 0.5 * (neighbourTransport%xzRay%Iout1+neighbourTransport%xyRay%Iout1)
                   Itmp2 = 0.5 * (neighbourTransport%xzRay%Iout2+neighbourTransport%xyRay%Iout2)
                   Itmp3 = 0.5 * (neighbourTransport%xzRay%Iout3+neighbourTransport%xyRay%Iout3)
                else
                   if (neighbourPattern%yzRayActive) then
                      Itmp1 = 0.5 * (neighbourTransport%yzRay%Iout1+neighbourTransport%xyRay%Iout1)
                      Itmp2 = 0.5 * (neighbourTransport%yzRay%Iout2+neighbourTransport%xyRay%Iout2)
                      Itmp3 = 0.5 * (neighbourTransport%yzRay%Iout3+neighbourTransport%xyRay%Iout3)
                   else
                      Itmp1 = neighbourTransport%xyRay%Iout1
                      Itmp2 = neighbourTransport%xyRay%Iout2
                      Itmp3 = neighbourTransport%xyRay%Iout3
                   endif
                endif

             end select
             if (transportCell%level.le.transportCell%xzNeighbour%level) then
                Iin1 = Itmp1
                Iin2 = Itmp2
                Iin3 = Itmp3
             else
                Iin1 = Itmp1
                Iin2 = Itmp2
                Iin3 = Itmp3
!         write(*,*) 'xz-interpolation not implemented yet', &
!              transportCell%pattern%xyRay%x0, transportCell%pattern%xyRay%z0
!         stop
             endif
          endif

          dpath = cellSizeAbsoluteUnits * transportCell%pattern%xzRay%len ! [cm]
          tau1 = transportCell%kappa1 * dpath
          tau2 = transportCell%kappa2 * dpath
          tau3 = transportCell%kappa3 * dpath
          tmpabs1 = exp(-tau1)
          tmpabs2 = exp(-tau2)
          tmpabs3 = exp(-tau3)
          if (tau1.gt.1.e-10) then
             tmpemi1 = (1.-tmpabs1)/transportCell%kappa1
          else
             tmpemi1 = dpath
          endif
          if (tau2.gt.1.e-10) then
             tmpemi2 = (1.-tmpabs2)/transportCell%kappa2
          else
             tmpemi2 = dpath
          endif
          if (tau3.gt.1.e-10) then
             tmpemi3 = (1.-tmpabs3)/transportCell%kappa3
          else
             tmpemi3 = dpath
          endif
          nemi1 = 0. ! emissivity [erg/s/cm^2/Hz/ster]
          nemi2 = 0.
          nemi3 = 0.
          transportCell%rt%xzRay%Iout1 = Iin1*tmpabs1 + nemi1*tmpemi1/dpath
          transportCell%rt%xzRay%Iout2 = Iin2*tmpabs2 + nemi2*tmpemi2/dpath
          transportCell%rt%xzRay%Iout3 = Iin3*tmpabs3 + nemi3*tmpemi3/dpath

          tmp = transportCell%rt%xyRay%Iout1+transportCell%rt%xyRay%Iout2+transportCell%rt%xyRay%Iout3
          if (tmp.lt.1.d-20 .and. tmp.gt.-1.d-20) then
             continue
          else
             print*, '-5-', transportCell%rt%xyRay%Iout1, transportCell%rt%xyRay%Iout2, transportCell%rt%xyRay%Iout3
             stop
          endif

!        Jmean1 = Jmean1 + Iin1 + transportCell%rt%xzRay%Iout1
!        Jmean2 = Jmean2 + Iin2 + transportCell%rt%xzRay%Iout2
!        Jmean3 = Jmean3 + Iin3 + transportCell%rt%xzRay%Iout3
!        imean = imean + 2

          call computeCellIntensity(Jmean1,Iin1,transportCell%rt%xzRay%Iout1)
          call computeCellIntensity(Jmean2,Iin2,transportCell%rt%xzRay%Iout2)
          call computeCellIntensity(Jmean3,Iin3,transportCell%rt%xzRay%Iout3)
          imean = imean + 1

!         if (transportCell%rt%xzRay%Iout1.eq.0.) then
!            write(*,*) 'b) zero Iout1 ...'
!            write(*,*) Iin1
!            stop
!         endif

       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!! yz-ray transport !!!!!!!!!!!!!!!!!!!!!!!!!!!
!     write(*,*) 'yz-ray transport'
       if (transportCell%pattern%yzRayActive) then
          if (.not.transportCell%yzNeighbourPresent) then
             Iin1 = uvb1
             Iin2 = uvb2
             Iin3 = uvb3
          else
             select case (transportCell%yzNeighbour%pattern%yzTop)
             case(xyEnd)
                Itmp1 = transportCell%yzNeighbour%rt%xyRay%Iout1
                Itmp2 = transportCell%yzNeighbour%rt%xyRay%Iout2
                Itmp3 = transportCell%yzNeighbour%rt%xyRay%Iout3
             case(xzEnd)
                if (.not.transportCell%yzNeighbour%pattern%xzRayActive) then
                   write(*,*) '4) Error: xzRay should be active ...', i, j, k
                   stop
                endif
                Itmp1 = transportCell%yzNeighbour%rt%xzRay%Iout1
                Itmp2 = transportCell%yzNeighbour%rt%xzRay%Iout2
                Itmp3 = transportCell%yzNeighbour%rt%xzRay%Iout3
             case(yzEnd)
                if (.not.transportCell%yzNeighbour%pattern%yzRayActive) then
                   write(*,*) '4) Error: yzRay should be active ...', i, j, k
                   stop
                endif
                Itmp1 = transportCell%yzNeighbour%rt%yzRay%Iout1
                Itmp2 = transportCell%yzNeighbour%rt%yzRay%Iout2
                Itmp3 = transportCell%yzNeighbour%rt%yzRay%Iout3
             case(0)
                if (transportCell%level.le.transportCell%yzNeighbour%level) then
                   write(*,*) '2) error in yzTop'
                   stop
                endif

                neighbourPattern => transportCell%yzNeighbour%pattern
                neighbourTransport => transportCell%yzNeighbour%rt
                if (neighbourPattern%xzRayActive) then
                   Itmp1 = 0.5 * (neighbourTransport%xzRay%Iout1+neighbourTransport%xyRay%Iout1)
                   Itmp2 = 0.5 * (neighbourTransport%xzRay%Iout2+neighbourTransport%xyRay%Iout2)
                   Itmp3 = 0.5 * (neighbourTransport%xzRay%Iout3+neighbourTransport%xyRay%Iout3)
                else
                   if (neighbourPattern%yzRayActive) then
                      Itmp1 = 0.5 * (neighbourTransport%yzRay%Iout1+neighbourTransport%xyRay%Iout1)
                      Itmp2 = 0.5 * (neighbourTransport%yzRay%Iout2+neighbourTransport%xyRay%Iout2)
                      Itmp3 = 0.5 * (neighbourTransport%yzRay%Iout3+neighbourTransport%xyRay%Iout3)
                   else
                      Itmp1 = neighbourTransport%xyRay%Iout1
                      Itmp2 = neighbourTransport%xyRay%Iout2
                      Itmp3 = neighbourTransport%xyRay%Iout3
                   endif
                endif

             end select
             if (transportCell%level.le.transportCell%yzNeighbour%level) then
                Iin1 = Itmp1
                Iin2 = Itmp2
                Iin3 = Itmp3
             else
                Iin1 = Itmp1
                Iin2 = Itmp2
                Iin3 = Itmp3
!         write(*,*) 'yz-interpolation not implemented yet', &
!              transportCell%rt%xyRay%y0, transportCell%rt%xyRay%z0
!         stop
             endif
          endif

          dpath = cellSizeAbsoluteUnits * transportCell%pattern%yzRay%len ! [cm]
          tau1 = transportCell%kappa1 * dpath
          tau2 = transportCell%kappa2 * dpath
          tau3 = transportCell%kappa3 * dpath
          tmpabs1 = exp(-tau1)
          tmpabs2 = exp(-tau2)
          tmpabs3 = exp(-tau3)
          if (tau1.gt.1.e-10) then
             tmpemi1 = (1.-tmpabs1)/transportCell%kappa1
          else
             tmpemi1 = dpath
          endif
          if (tau2.gt.1.e-10) then
             tmpemi2 = (1.-tmpabs2)/transportCell%kappa2
          else
             tmpemi2 = dpath
          endif
          if (tau3.gt.1.e-10) then
             tmpemi3 = (1.-tmpabs3)/transportCell%kappa3
          else
             tmpemi3 = dpath
          endif
          nemi1 = 0. ! emissivity [erg/s/cm^2/Hz/ster]
          nemi2 = 0.
          nemi3 = 0.
          transportCell%rt%yzRay%Iout1 = Iin1*tmpabs1 + nemi1*tmpemi1/dpath
          transportCell%rt%yzRay%Iout2 = Iin2*tmpabs2 + nemi2*tmpemi2/dpath
          transportCell%rt%yzRay%Iout3 = Iin3*tmpabs3 + nemi3*tmpemi3/dpath

          tmp = transportCell%rt%xyRay%Iout1+transportCell%rt%xyRay%Iout2+transportCell%rt%xyRay%Iout3
          if (tmp.lt.1.d-20 .and. tmp.gt.-1.d-20) then
             continue
          else
             print*, '-6-', transportCell%rt%xyRay%Iout1, transportCell%rt%xyRay%Iout2, transportCell%rt%xyRay%Iout3
             stop
          endif

!        Jmean1 = Jmean1 + Iin1 + transportCell%rt%yzRay%Iout1
!        Jmean2 = Jmean2 + Iin2 + transportCell%rt%yzRay%Iout2
!        Jmean3 = Jmean3 + Iin3 + transportCell%rt%yzRay%Iout3
!        imean = imean + 2

          call computeCellIntensity(Jmean1,Iin1,transportCell%rt%yzRay%Iout1)
          call computeCellIntensity(Jmean2,Iin2,transportCell%rt%yzRay%Iout2)
          call computeCellIntensity(Jmean3,Iin3,transportCell%rt%yzRay%Iout3)
          imean = imean + 1

!         if (transportCell%rt%yzRay%Iout1.eq.0.) then
!            write(*,*) 'c) zero Iout1 ...'
!            write(*,*) Iin1
!            stop
!         endif

       endif
!     write(*,*) 'yz-ray done'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       transportCell%Jmean1 = transportCell%Jmean1 + Jmean1/float(imean) * weight
       transportCell%Jmean2 = transportCell%Jmean2 + Jmean2/float(imean) * weight
       transportCell%Jmean3 = transportCell%Jmean3 + Jmean3/float(imean) * weight

!     transportCell%Jabs1 = transportCell%Jabs1 + Jabs1 * weight
!     transportCell%Jabs2 = transportCell%Jabs2 + Jabs2 * weight
!     transportCell%Jabs3 = transportCell%Jabs3 + Jabs3 * weight

    endif

  end subroutine transport

! recursive subroutine copyBaseToTrans(transCell,baseCell,is,js,ks)

!   use definitions

!   implicit none
!   type(zoneType), target :: transCell
!   type(zoneType) :: baseCell
!   integer, dimension(2,2,2), intent(in) :: is, js, ks
!   integer :: i, j, k

!   transCell = baseCell
!   if (baseCell%refined) then
!      allocate(transCell%cell(2,2,2))
!      do i = 1, 2
!         do j = 1, 2
!            do k = 1, 2
!               call copyBaseToTrans(transCell%cell(i,j,k), &
!                    baseCell%cell(is(i,j,k),js(i,j,k),ks(i,j,k)), &
!                    is,js,ks)
!               transCell%cell(i,j,k)%parent => transCell
!            enddo
!         enddo
!      enddo
!   endif

! end subroutine copyBaseToTrans

! recursive subroutine copyTransToBase(transCell,baseCell,is,js,ks)

!   use definitions

!   implicit none
!   type(zoneType), target :: transCell
!   type(zoneType) :: baseCell
!   integer, dimension(2,2,2), intent(in) :: is, js, ks
!   integer :: i, j, k

!   if (transCell%refined) then
!      do i = 1, 2
!         do j = 1, 2
!            do k = 1, 2
!               call copyTransToBase(transCell%cell(i,j,k), &
!                    baseCell%cell(is(i,j,k),js(i,j,k),ks(i,j,k)),is,js,ks)
!            enddo
!         enddo
!      enddo
!   else
!      baseCell%Jmean1 = baseCell%Jmean1 + transCell%Jmean1
!   endif

! end subroutine copyTransToBase

  recursive subroutine patternNullify(pattern)

    use definitions

    implicit none
    type(patternType) :: pattern
    integer :: i, j, k

    do i = 1, 2
       do j = 1, 2
          do k = 1, 2
             if (pattern%cell(i,j,k)%refined) call patternNullify(pattern%cell(i,j,k))
          enddo
       enddo
    enddo
    deallocate(pattern%cell)

  end subroutine patternNullify

  subroutine computeCellIntensity(Jmean,Iin,Iout)

    use definitions, only: RealKind

    implicit none
    real(kind=RealKind), intent(inout) :: Jmean
    real(kind=RealKind), intent(in) :: Iin, Iout

    if (Iout.lt.Iin) then
       Jmean = Jmean + (Iin - Iout)/log(Iin/Iout)
    else
       Jmean = Jmean + 0.5 * (Iin + Iout)
    endif

!  Jmean = Jmean + Iin

!  Jabs = Jabs + (Iin - Iout)

  end subroutine computeCellIntensity

  recursive subroutine assignUvbRadiation(currentCell)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    integer :: i, j, k
    real(kind=RealKind) :: HI, HeI, HeII, meanFreePathLymanLimit

    currentCell%Jmean1 = 0.
    currentCell%Jmean2 = 0.
    currentCell%Jmean3 = 0.

    if (currentCell%refined) then
       do i = 1, 2
          do j = 1, 2
             do k = 1, 2
                call assignUvbRadiation(currentCell%cell(i,j,k))
             enddo
          enddo
       enddo
    else
       HI = min(currentCell%HI,psi*currentCell%rho/mh)
       HeI = currentCell%HeI
       HeII = currentCell%HeII
       meanFreePathLymanLimit = 1./(HI*6.3e-18 + HeI*7.42e-18 + HeII*1.58e-18)
       if (meanFreePathLymanLimit.ge.selfShieldingThreshold) then
          currentCell%Jmean1 = uvb1
          currentCell%Jmean2 = uvb2
          currentCell%Jmean3 = uvb3
       else
          currentCell%Jmean1 = 0.
          currentCell%Jmean2 = 0.
          currentCell%Jmean3 = 0.
       endif
    endif

  end subroutine assignUvbRadiation

end module transportRoutinesModule
