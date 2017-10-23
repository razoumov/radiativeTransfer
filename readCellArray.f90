program readCellArray

  use definitions
  use rotateIndicesModule

  implicit none
  integer :: i, j, k, icell, jcell, kcell, nxmap, nymap, i0, j0, k0, lmax, nx, ny, nz
  integer(1) :: izone
  real(kind=RealKind) :: xnew, ynew, znew, x0, y0, zslice, z0
  integer, dimension(3) :: baseGridSize
  integer, dimension(2,2,2) :: is, js, ks

  integer sfstart, sffinfo, sfselect, sfginfo, sfrdata, sfendacc, sfend, sfcreate, sfwdata
  real*4, dimension(:,:), pointer :: mapArray
  character(60) :: dirname

  dirname = '/path/to/K15-06-64/3.6-rot/'
  sd_id = sfstart(trim(dirname)//'cellArray0016.h4', dfacc_read)
  nxmap = 4096
  nymap = nxmap
  izone = 2
  zslice = 0.506896972656250

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

  sds_id = sfselect(sd_id, 1)
  status = sfginfo(sds_id,sds_name,rank,dim_sizes,data_type,n_attrs)
  if (status.ne.0 .or. rank.ne.1) then
     write(*,*) 'error reading ...'
     stop
  endif
  ncosmic = dim_sizes(1)
  edges(1) = dim_sizes(1)

  allocate (cellArrayLevel(ncosmic))
  allocate (cellArrayHI(ncosmic))
  allocate(cellArrayTemp(ncosmic))
  allocate(cellArrayDensity(ncosmic))

  status = sfrdata(sds_id, start, stride, edges, cellArrayLevel)
  status = sfendacc(sds_id)

  sds_id = sfselect(sd_id, 2)
  status = sfrdata(sds_id, start, stride, edges, cellArrayHI)
  status = sfendacc(sds_id)

  sds_id = sfselect(sd_id, 5)
  status = sfrdata(sds_id, start, stride, edges, cellArrayTemp)
  status = sfendacc(sds_id)

  sds_id = sfselect(sd_id, 6)
  status = sfrdata(sds_id, start, stride, edges, cellArrayDensity)
  status = sfendacc(sds_id)

  status = sfend(sd_id)

  print*, ncosmic
  write(*,*) cellArrayLevel(1)
  write(*,*) cellArrayHI(1)
  write(*,*) cellArrayTemp(1)
  write(*,*) cellArrayDensity(1)

!   open(14,file=trim(dirname)//'cellArray.dat',status='replace',form='unformatted')
!   write(14) (cellArrayLevel(i),i=1,ncosmic)
!   write(14) (cellArrayHI(i),i=1,ncosmic)
!   write(14) (cellArrayTemp(i),i=1,ncosmic)
!   write(14) (cellArrayDensity(i),i=1,ncosmic)
!   close(14)

  allocate(baseGrid%cell(nx,ny,nz))

  icosmic = 0
  lmax = 0
  do i = 1, nx
     do j = 1, ny
        do k = 1, nz
           baseGrid%cell(i,j,k)%refined = .false.
           call createFullyThreadedStructure(baseGrid%cell(i,j,k),0,lmax)
        enddo
     enddo
  enddo

  write(*,*) 'maximum level of refinement =', lmax

  deallocate(cellArrayLevel,cellArrayHI)
  stop

!   densestCell => baseGrid%cell(1,1,1)
!   do i = 1, nx
!      x0 = (float(i)-0.5)/float(nx)
!      do j = 1, ny
!         y0 = (float(j)-0.5)/float(ny)
!         do k = 1, nz
!            z0 = (float(k)-0.5)/float(nz)
!            call findDensestCell(baseGrid%cell(i,j,k),x0,y0,z0,1./float(nx))
!         enddo
!      enddo
!   enddo
!   print*, 'peak:', densestCell%HI, xbase, ybase, zbase

! ------------------------ plot HI physical density in a slice ----------------

  do i = 1, 2
     do j = 1, 2
        do k = 1, 2
           call rotateIndices(i,j,k,2,2,2,izone,is(i,j,k),js(i,j,k),ks(i,j,k))
        enddo
     enddo
  enddo
  allocate(mapArray(nxmap,nymap))
  do i = 1, nxmap
     print*, 'doing', i, ' out of', nxmap
     x0 = (float(i)-0.5)/float(nxmap)
     i0 = int(x0*nx) + 1
     do j = 1, nymap
        y0 = (float(j)-0.5)/float(nymap)
        j0 = int(y0*ny) + 1
        xnew = x0*float(nx) - float(i0-1)
        ynew = y0*float(ny) - float(j0-1)
        mapArray(i,j) = 0.
        k0 = int(zslice*nz) + 1
        znew = zslice*float(nz) - float(k0-1)
        call rotateIndices(i0,j0,k0,nx,ny,nz,izone,icell,jcell,kcell)
        call sliceCell(baseGrid%cell(icell,jcell,kcell),mapArray(i,j), &
             xnew,ynew,znew,is,js,ks)
     enddo
  enddo
  sd_id = sfstart(trim(dirname)//'map.h4',dfacc_create)
  edges(1) = nxmap
  edges(2) = nymap
  sds_id = sfcreate(sd_id,'map',dfnt_float32,2,edges)
  start = 0
  stride = 1
  status = sfwdata(sds_id,start,stride,edges,mapArray)
  status = sfendacc(sds_id)
  status = sfend(sd_id)
  deallocate(mapArray)

contains

  recursive subroutine createFullyThreadedStructure(currentCell,level,lmax)

    use definitions

    implicit none
    type(zoneType), target :: currentCell
    integer, intent(in) :: level
    integer, intent(inout) :: lmax
    integer :: i, j, k

    icosmic = icosmic + 1
    lmax = max(cellArrayLevel(icosmic),lmax)
    if (cellArrayLevel(icosmic).eq.level) then
       currentCell%HI = dble(cellArrayHI(icosmic))
    else
       if (cellArrayLevel(icosmic).gt.level) then
          allocate(currentCell%cell(2,2,2))
          currentCell%refined = .true.
          icosmic = icosmic - 1
          do i = 1, 2
             do j = 1, 2
                do k = 1, 2
                   currentCell%cell(i,j,k)%refined = .false.
                   call createFullyThreadedStructure(currentCell%cell(i,j,k),level+1,lmax)
                enddo
             enddo
          enddo
       else
          write(*,*) 'error in levels', icosmic, cellArrayLevel(icosmic), level
          stop
       endif
    endif

  end subroutine createFullyThreadedStructure

  recursive subroutine sliceCell(currentCell,currentPixel, &
       x0,y0,z0,is,js,ks)

    use definitions

    implicit none
    type(zoneType) :: currentCell
    real*4, intent(inout) :: currentPixel
    real(kind=RealKind), intent(in) :: x0, y0, z0
    integer, dimension(2,2,2), intent(in) :: is, js, ks
    integer :: i, j, k
    real(kind=RealKind) :: xnew, ynew, znew

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
       if (z0.lt.0.5) then
          k = 1
          znew = 2.*z0
       else
          k = 2
          znew = 2.*z0 - 1.
       endif
       call sliceCell(currentCell%cell(is(i,j,k),js(i,j,k),ks(i,j,k)),currentPixel, &
            xnew,ynew,znew,is,js,ks)
    else
       currentPixel = currentCell%HI
    endif

  end subroutine sliceCell

  recursive subroutine findDensestCell(currentCell,x0,y0,z0,cellSize)

    use definitions

    implicit none
    type(zoneType), target :: currentCell
    real(kind=RealKind), intent(in) :: x0, y0, z0, cellSize
    integer :: i, j, k
    real(kind=RealKind) :: xnew, ynew, znew

    if (currentCell%refined) then
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
                call findDensestCell(currentCell%cell(i,j,k),xnew,ynew,znew,cellSize/2.)
             enddo
          enddo
       enddo
    else
       if (currentCell%HI.gt.densestCell%HI) then
          densestCell => currentCell
          xbase = x0
          ybase = y0
          zbase = z0
       endif
    endif

  end subroutine findDensestCell

end program readCellArray
