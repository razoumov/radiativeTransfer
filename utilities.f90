module utilities

  use definitions

  interface swap
     module procedure swap_i, swap_r, swap_s
  end interface

contains

  subroutine sortStellarParticles(arr)

! based on heapsort from numerical recipes

    implicit none
    type(starType), dimension(:), intent(inout) :: arr
    integer :: i, n

    n = size(arr)
    do i=n/2,1,-1
       call sift_down(i,n,arr)
    end do
    do i=n,2,-1
       call swap(arr(1),arr(i))
       call sift_down(1,i-1,arr)
    end do

  end subroutine sortStellarParticles

  subroutine sift_down(l,r,arr)

    implicit none
    integer, intent(in) :: l,r
    type(starType), dimension(:), intent(inout) :: arr
    integer :: j,jold
    type(starType) :: a

    a=arr(l)
    jold=l
    j=l+l
    do
       if (j > r) exit
       if (j < r) then
          if (arr(j)%location < arr(j+1)%location) j=j+1
       end if
       if (a%location >= arr(j)%location) exit
       arr(jold)=arr(j)
       jold=j
       j=j+j
    end do
    arr(jold)=a

  end subroutine sift_down

  subroutine swap_i(a,b)
    integer, intent(inout) :: a,b
    integer :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_i

  subroutine swap_r(a,b)
    real(kind=RealKind), intent(inout) :: a,b
    real(kind=RealKind) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_r

  subroutine swap_s(a,b)
    type(starType), intent(inout) :: a,b
    type(starType) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_s

  subroutine ratint(xa,ya,x,y,dy)

    implicit none
    real(kind=RealKind), dimension(:), intent(in) :: xa,ya
    real(kind=RealKind), intent(in) :: x
    real(kind=RealKind), intent(out) :: y,dy
    integer :: m,n,ns
    real(kind=RealKind), dimension(size(xa)) :: c,d,dd,h,t
    real(kind=RealKind), parameter :: tiny=1.e-25

    n=assert_eq(size(xa),size(ya),'ratint')
    h=xa-x
    ns=iminloc(abs(h))
    y=ya(ns)
    if (x == xa(ns)) then
       dy=0.0
       return
    end if
    c=ya
    d=ya+tiny
    ns=ns-1
    do m=1,n-1
       t(1:n-m)=(xa(1:n-m)-x)*d(1:n-m)/h(1+m:n)
       dd(1:n-m)=t(1:n-m)-c(2:n-m+1)
       if (any(dd(1:n-m) == 0.0)) &
            call nrerror('failure in ratint')
       dd(1:n-m)=(c(2:n-m+1)-d(1:n-m))/dd(1:n-m)
       d(1:n-m)=c(2:n-m+1)*dd(1:n-m)
       c(1:n-m)=t(1:n-m)*dd(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do

  end subroutine ratint

  subroutine nrerror(string)
    implicit none
    character(len=*), intent(in) :: string
    write (*,*) 'nrerror: ',string
    stop 'program terminated by nrerror'
  end subroutine nrerror

  function iminloc(arr)
    implicit none
    real(kind=RealKind), dimension(:), intent(in) :: arr
    integer, dimension(1) :: imin
    integer :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  end function iminloc

  function assert_eq(n1,n2,string)
    implicit none
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1,n2
    integer :: assert_eq
    if (n1 == n2) then
       assert_eq=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', string
       stop 'program terminated by assert_eq'
    end if
  end function assert_eq

end module utilities
