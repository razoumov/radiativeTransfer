module rotateIndicesModule

  use definitions

contains

  subroutine rotateIndices(i,j,k,nx,ny,nz,izone,icell,jcell,kcell)

    implicit none
    integer, intent(in) :: i, j, k, nx, ny, nz
    integer(1), intent(in) :: izone
    integer, intent(out) :: icell, jcell, kcell

    select case (izone)
    case(1)
       icell = i
       jcell = j
       kcell = k
    case(2)
       icell = j
       jcell = k
       kcell = i
    case(3)
       icell = k
       jcell = i
       kcell = j
    case(4)
       icell = i
       jcell = k
       kcell = nz + 1 - j
    case(5)
       icell = j
       jcell = i
       kcell = nz + 1 - k
    case(6)
       icell = k
       jcell = j
       kcell = nz + 1 - i
    case(7)
       icell = i
       jcell = ny + 1 - j
       kcell = nz + 1 - k
    case(8)
       icell = j
       jcell = ny + 1 - k
       kcell = nz + 1 - i
    case(9)
       icell = k
       jcell = ny + 1 - i
       kcell = nz + 1 - j
    case(10)
       icell = i
       jcell = ny + 1 - k
       kcell = j
    case(11)
       icell = j
       jcell = ny + 1 - i
       kcell = k
    case(12)
       icell = k
       jcell = ny + 1 - j
       kcell = i
    case(13)
       icell = nx + 1 - i
       jcell = j
       kcell = k
    case(14)
       icell = nx + 1 - j
       jcell = k
       kcell = i
    case(15)
       icell = nx + 1 - k
       jcell = i
       kcell = j
    case(16)
       icell = nx + 1 - i
       jcell = k
       kcell = nz + 1 - j
    case(17)
       icell = nx + 1 - j
       jcell = i
       kcell = nz + 1 - k
    case(18)
       icell = nx + 1 - k
       jcell = j
       kcell = nz + 1 - i
    case(19)
       icell = nx + 1 - i
       jcell = ny + 1 - j
       kcell = nz + 1 - k
    case(20)
       icell = nx + 1 - j
       jcell = ny + 1 - k
       kcell = nz + 1 - i
    case(21)
       icell = nx + 1 - k
       jcell = ny + 1 - i
       kcell = nz + 1 - j
    case(22)
       icell = nx + 1 - i
       jcell = ny + 1 - k
       kcell = j
    case(23)
       icell = nx + 1 - j
       jcell = ny + 1 - i
       kcell = k
    case(24)
       icell = nx + 1 - k
       jcell = ny + 1 - j
       kcell = i
    end select

  end subroutine rotateIndices

end module rotateIndicesModule
