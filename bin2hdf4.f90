program bin2hdf4

  use definitions

  implicit none
  integer :: nlevels, icell, level, i, j, k, itmp, i0, j0, k0
  real(kind=RealKind) :: tmp, tmp1, tmp2, x0, y0, z0, xnew, ynew, znew, &
       xa, xb, ya, yb, za, zb, phi, theta, weight, time, dt, xneu, nh, nhe

  type(readLevelType), dimension(:), pointer :: readLevel

  type(zoneType), pointer :: currentCell

  type(pixelType), target :: sphere
  type(pixelType), pointer :: currentPixel, tmpPixel, leafPixel

  integer :: sfstart, sffinfo, sfselect, sfginfo, sfrdata, sfendacc, sfend, sfcreate, sfwdata
  integer, dimension(3) :: baseGridSize

  character(60), parameter :: dirname = '11.63/'
  character(60), parameter :: filename = 'S29COSMO-055-8.1200kpc_128_z11.63_reionz10_0.9_velmet'

  readMetals = .false.
  itmp = len(trim(filename))
  do i = 1, itmp
     if (filename(i:i).eq.'m'.and.i.lt.itmp-1) then
        if (filename(i:i+2).eq.'met') readMetals = .true.
     endif
  enddo

  readKinematics = .false.
  itmp = len(trim(filename))
  do i = 1, itmp
     if (filename(i:i).eq.'v'.and.i.lt.itmp-1) then
        if (filename(i:i+2).eq.'vel') readKinematics = .true.
     endif
  enddo

  ! --------------- read data in binary format

  open(14,file=trim(dirname)//trim(filename)//'.dat',status='old',form='unformatted')

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
     if (readMetals) allocate(readLevel(level)%abun(readLevel(level)%ncell,4))
     if (readKinematics) allocate(readLevel(level)%vel(readLevel(level)%ncell,3))

     read(14) (readLevel(level)%pos(icell,1),icell=1,readLevel(level)%ncell)
     read(14) (readLevel(level)%pos(icell,2),icell=1,readLevel(level)%ncell)
     read(14) (readLevel(level)%pos(icell,3),icell=1,readLevel(level)%ncell)
     read(14) (readLevel(level)%lT(icell),icell=1,readLevel(level)%ncell)
     read(14) (readLevel(level)%lnH(icell),icell=1,readLevel(level)%ncell)
     read(14) (readLevel(level)%lx(icell),icell=1,readLevel(level)%ncell)
     if (readMetals) then
        read(14) (readLevel(level)%abun(icell,1),icell=1,readLevel(level)%ncell)
        read(14) (readLevel(level)%abun(icell,2),icell=1,readLevel(level)%ncell)
        read(14) (readLevel(level)%abun(icell,3),icell=1,readLevel(level)%ncell)
        read(14) (readLevel(level)%abun(icell,4),icell=1,readLevel(level)%ncell)

!         do icell = 1, readLevel(level)%ncell
!            if (readLevel(level)%abun(icell,2).gt.0.) then
!               print*, '---', level, readLevel(level)%abun(icell,1:4)
!            endif
!         enddo

     endif
     if (readKinematics) then
        read(14) (readLevel(level)%vel(icell,1),icell=1,readLevel(level)%ncell)
        read(14) (readLevel(level)%vel(icell,2),icell=1,readLevel(level)%ncell)
        read(14) (readLevel(level)%vel(icell,3),icell=1,readLevel(level)%ncell)
     endif

  enddo

  close(14)

!   level = 1
!   do icell = 1, readLevel(level)%ncell
!      print*, readLevel(level)%pos(icell,1), readLevel(level)%pos(icell,2), readLevel(level)%pos(icell,3)
!      print*, readLevel(level)%vel(icell,1), readLevel(level)%vel(icell,2), readLevel(level)%vel(icell,3)
!      print*, readLevel(level)%abun(icell,1:4)
!   enddo

  do while (readLevel(nlevels)%ncell.eq.0)
     nlevels = nlevels - 1
  enddo

  do level = 1, nlevels
     write(*,*) 'actual level =', level, readLevel(level)%ncell
  enddo

  ! --------------- write data in hdf4 format

  sd_id = sfstart(trim(dirname)//trim(filename)//'.h4',dfacc_create)

  edges(1) = 1
  start = 0
  stride = 1
  sds_id = sfcreate(sd_id,'nlevels',dfnt_int32,1,edges)
  status = sfwdata(sds_id,start,stride,edges,nlevels)
  status = sfendacc(sds_id)

  do level = 1, nlevels

     edges(1:2) = (/ readLevel(level)%ncell, 3 /)
     start = 0
     stride = 1
     sds_id = sfcreate(sd_id,'pos',dfnt_float32,2,edges)
     status = sfwdata(sds_id,start,stride,edges,readLevel(level)%pos)
     status = sfendacc(sds_id)

     edges(1) = readLevel(level)%ncell
     start = 0
     stride = 1
     sds_id = sfcreate(sd_id,'lT',dfnt_float32,1,edges)
     status = sfwdata(sds_id,start,stride,edges,readLevel(level)%lT)
     status = sfendacc(sds_id)

     edges(1) = readLevel(level)%ncell
     start = 0
     stride = 1
     sds_id = sfcreate(sd_id,'lnH',dfnt_float32,1,edges)
     status = sfwdata(sds_id,start,stride,edges,readLevel(level)%lnH)
     status = sfendacc(sds_id)

     edges(1) = readLevel(level)%ncell
     start = 0
     stride = 1
     sds_id = sfcreate(sd_id,'lx',dfnt_float32,1,edges)
     status = sfwdata(sds_id,start,stride,edges,readLevel(level)%lx)
     status = sfendacc(sds_id)

     if (readMetals) then
        edges(1:2) = (/ readLevel(level)%ncell, 4 /)
        start = 0
        stride = 1
        sds_id = sfcreate(sd_id,'abun',dfnt_float32,2,edges)
        status = sfwdata(sds_id,start,stride,edges,readLevel(level)%abun)
        status = sfendacc(sds_id)
     endif

     if (readKinematics) then
        edges(1:2) = (/ readLevel(level)%ncell, 3 /)
        start = 0
        stride = 1
        sds_id = sfcreate(sd_id,'vel',dfnt_float32,2,edges)
        status = sfwdata(sds_id,start,stride,edges,readLevel(level)%vel)
        status = sfendacc(sds_id)
     endif

  enddo

  status = sfend(sd_id)

  do level = 1, nlevels
     deallocate(readLevel(level)%pos)
     deallocate(readLevel(level)%lT)
     deallocate(readLevel(level)%lnH)
     deallocate(readLevel(level)%lx)
     if (readMetals) deallocate(readLevel(level)%abun)
     if (readKinematics) deallocate(readLevel(level)%vel)
  enddo
  deallocate(readLevel)

end program bin2hdf4
