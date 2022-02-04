! ###################################################################
! Copyright (c) 2013-2022, Marc De Graef Research Group/Carnegie Mellon University
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are
! permitted provided that the following conditions are met:
!
!     - Redistributions of source code must retain the above copyright notice, this list
!        of conditions and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright notice, this
!        list of conditions and the following disclaimer in the documentation and/or
!        other materials provided with the distribution.
!     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names
!        of its contributors may be used to endorse or promote products derived from
!        this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ###################################################################

!--------------------------------------------------------------------------
! EMsoft:math.f90
!--------------------------------------------------------------------------
!
! MODULE: math
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief collection of mathematical/numerical routines that don't fit anywhere else
!
!> @date 10/13/98 MDG 1.0 original
!> @date 05/19/01 MDG 2.0 f90
!> @date 11/27/01 MDG 2.1 added kind support
!> @date 03/19/13 MDG 3.0 updated all routines
!> @date 11/13/13 MDG 4.0 added MatrixExponential routine
!> @date 11/23/15 MDG 4.1 moved several routines from other mods into this one
!> @date 10/24/17 MDG 4.2 added infty()/inftyd() functions to return the IEEE infinity value
!> @date 08/23/19 MDG 4.3 removed spaces around "kind" statements to facilitate f90wrap python wrapper generation
!> @date 10/04/19 MDG 4.4 adds vecnorm to replace non-standard NORM2 calls  (F2003 compliance)
!> @date 10/04/19 MDG 4.5 adds nan() function, returning a single or double precision IEEE NaN value
!> @date 11/01/19 MDG 4.6 adds Jaccard_Distance routine (moved from Indexingmod)
!--------------------------------------------------------------------------
! ###################################################################
!

module mod_math
  !! author: MDG
  !! version: 1.0
  !! date: 01/17/20
  !!
  !! collection of mathematical/numerical routines that don't fit anywhere else (no classes, just routines)

use mod_kinds
use mod_global

! public :: mInvert, cross3, infty, inftyd, nan, nan_d, 

interface mInvert
        module procedure mInvert
        module procedure mInvert_d
end interface

interface cross3
        module procedure cross3
        module procedure cross3_d
end interface

interface vecnorm
        module procedure vecnorm
        module procedure vecnorm_d
        module procedure vecnorm2
        module procedure vecnorm2_d
end interface

contains

!--------------------------------------------------------------------------
recursive function vecnorm(vec) result(veclen)
!DEC$ ATTRIBUTES DLLEXPORT :: vecnorm
  !! author: MDG
  !! version: 1.0
  !! date: 01/17/20
  !!
  !! return the single precision length of a 1D vector

real(kind=sgl),INTENT(IN)        :: vec(:)
real(kind=sgl)                   :: veclen

integer(kind=irg)                :: sz(1)

sz = size(vec)

veclen = sqrt(sum(vec(1:sz(1))*vec(1:sz(1))))

end function vecnorm

!--------------------------------------------------------------------------
recursive function vecnorm_d(vec) result(veclen)
  !DEC$ ATTRIBUTES DLLEXPORT :: vecnorm_d
  !! author: MDG
  !! version: 1.0
  !! date: 01/17/20
  !!
  !! return the double precision length of a 1D vector

real(kind=dbl),INTENT(IN)       :: vec(:)
real(kind=dbl)                  :: veclen

integer(kind=irg)               :: sz(1)

sz = size(vec)

veclen = sqrt(sum(vec(1:sz(1))*vec(1:sz(1))))

end function vecnorm_d

!--------------------------------------------------------------------------
recursive function vecnorm2(vec) result(veclen)
!DEC$ ATTRIBUTES DLLEXPORT :: vecnorm2
  !! author: MDG
  !! version: 1.0
  !! date: 01/17/20
  !!
  !! return the single precision length of a 2D array

real(kind=sgl),INTENT(IN)        :: vec(:,:)
real(kind=sgl)                   :: veclen

integer(kind=irg)                :: sz(2)

sz = size(vec)

veclen = sqrt(sum(vec(1:sz(1),1:sz(2))*vec(1:sz(1),1:sz(2))))

end function vecnorm2

!--------------------------------------------------------------------------
recursive function vecnorm2_d(vec) result(veclen)
!DEC$ ATTRIBUTES DLLEXPORT :: vecnorm2_d
  !! author: MDG
  !! version: 1.0
  !! date: 01/17/20
  !!
  !! return the double precision length of a 2D array

real(kind=dbl),INTENT(IN)        :: vec(:,:)
real(kind=dbl)                   :: veclen

integer(kind=irg)                :: sz(2)

sz = size(vec)

veclen = sqrt(sum(vec(1:sz(1),1:sz(2))*vec(1:sz(1),1:sz(2))))

end function vecnorm2_d

!--------------------------------------------------------------------------
recursive function infty() result(infinity)
!DEC$ ATTRIBUTES DLLEXPORT :: infty
  !! author: MDG
  !! version: 1.0
  !! date: 01/17/20
  !!
  !! return the single precision IEEE value for infinity

real(kind=sgl)      :: infinity
real(kind=sgl)      :: big

big = HUGE(1.0)
infinity = big + HUGE(1.0)

end function infty

!--------------------------------------------------------------------------
recursive function inftyd() result(infinity)
!DEC$ ATTRIBUTES DLLEXPORT :: inftyd
  !! author: MDG
  !! version: 1.0
  !! date: 01/17/20
  !!
  !! return the double precision IEEE value for infinity

real(kind=dbl)      :: infinity
real(kind=dbl)      :: big

big = HUGE(1.D0)
infinity = big + HUGE(1.D0)

end function inftyd

!--------------------------------------------------------------------------
recursive function nan() result(x)
!DEC$ ATTRIBUTES DLLEXPORT :: nan
  !! author: MDG
  !! version: 1.0
  !! date: 01/17/20
  !!
  !! return the single precision IEEE value for nan

 use, intrinsic :: iso_fortran_env
 use, intrinsic :: ieee_arithmetic

 IMPLICIT NONE

real(kind=sgl)        :: x

x = ieee_value(x, ieee_quiet_nan)

end function nan

!--------------------------------------------------------------------------
recursive function nan_d() result(x)
!DEC$ ATTRIBUTES DLLEXPORT :: nan_d
  !! author: MDG
  !! version: 1.0
  !! date: 01/17/20
  !!
  !! return the double precision IEEE value for nan

 use, intrinsic :: iso_fortran_env
 use, intrinsic :: ieee_arithmetic

 IMPLICIT NONE

real(kind=dbl)        :: x

x = ieee_value(x, ieee_quiet_nan)

end function nan_d

!--------------------------------------------------------------------------
!
! FUNCTION: cross3
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief  cross product of two 3D vector in the order of input
!
!
!> @param u input vector 1
!> @param v input vector 2
!
!> @date 03/03/16   SS 1.0 original
!> @date 12/01/16  MDG 1.1 split in single and double precision versions
!--------------------------------------------------------------------------
recursive function cross3(u, v) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: cross3

IMPLICIT NONE

real(kind=sgl),INTENT(IN)      :: u(3)
real(kind=sgl),INTENT(IN)      :: v(3)
real(kind=sgl)                 :: res(3)

res(1) = u(2)*v(3) - u(3)*v(2)
res(2) = u(3)*v(1) - u(1)*v(3)
res(3) = u(1)*v(2) - u(2)*v(1)

end function cross3

!--------------------------------------------------------------------------
!
! FUNCTION: cross3_d
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief  cross product of two 3D vector in the order of input (double precision)
!
!
!> @param u input vector 1
!> @param v input vector 2
!
!> @date 03/03/16   SS 1.0 original
!--------------------------------------------------------------------------
recursive function cross3_d(u, v) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: cross3_d

IMPLICIT NONE

real(kind=dbl),INTENT(IN)      :: u(3)
real(kind=dbl),INTENT(IN)      :: v(3)
real(kind=dbl)                 :: res(3)

res(1) = u(2)*v(3) - u(3)*v(2)
res(2) = u(3)*v(1) - u(1)*v(3)
res(3) = u(1)*v(2) - u(2)*v(1)

end function cross3_d


!--------------------------------------------------------------------------
recursive subroutine mInvert_d(a,b,uni)
!DEC$ ATTRIBUTES DLLEXPORT :: mInvert_d
  !! author: MDG
  !! version: 1.0
  !! date: 01/17/20
  !!
  !! Invert a 3x3 matrix; if unitary, simply transpose

use mod_io

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: a(3,3)
 !! input matrix
real(kind=dbl),INTENT(OUT)              :: b(3,3)
 !! output matrix
logical,INTENT(IN)                      :: uni
 !! unitary logical

type(IO_T)                              :: Message
real(kind=dbl)                          :: d
integer(kind=irg)                       :: i, j

! it is a regular (non-unitary) matrix
 if (.not.uni) then
  d = a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+ &
         a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)- &
         a(1,2)*a(2,1)*a(3,3)-a(1,1)*a(2,3)*a(3,2)
  if (d.ne.0.D0) then
   b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
   b(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
   b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
   b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
   b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
   b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
   b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
   b(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)
   b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
   b = b/d
  else
    do i=1,3
      write (*,*) (a(i,j),j=1,3)
    end do
!  call Message%printError('mInvert','matrix has zero determinant')
   call Message%printMessage('mInvert: matrix has zero determinant')
   b = a
  end if
 else
! it is a unitary matrix, so simply get the transpose
  b = transpose(a)
 endif

end subroutine mInvert_d

!--------------------------------------------------------------------------
recursive subroutine mInvert(a,b,uni)
!DEC$ ATTRIBUTES DLLEXPORT :: mInvert
  !! author: MDG
  !! version: 1.0
  !! date: 01/17/20
  !!
  !! Invert a single precision 3x3 matrix; if unitary, simply transpose

use mod_io

IMPLICIT NONE

real(kind=sgl),INTENT(IN)               :: a(3,3)
 !! input matrix
real(kind=sgl),INTENT(OUT)              :: b(3,3)
 !! output matrix
logical,INTENT(IN)                      :: uni
 !! unitary logical

type(IO_T)                              :: Message
real(kind=sgl)                          :: d                    !< auxiliary variable
integer(kind=irg)                       :: i, j
! it is a regular (non-unitary) matrix
 if (.not.uni) then
  d = a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+ &
         a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)- &
         a(1,2)*a(2,1)*a(3,3)-a(1,1)*a(2,3)*a(3,2)
  if (d.ne.0.0) then
   b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
   b(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
   b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
   b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
   b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
   b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
   b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
   b(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)
   b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
   b = b/d
  else
    do i=1,3
      write (*,*) (a(i,j),j=1,3)
    end do
!  call Message%printError('mInvert','matrix has zero determinant')
   call Message%printMessage('mInvert: matrix has zero determinant')
   b = a
  end if
 else
! it is a unitary matrix, so simply get the transpose
  b = transpose(a)
 endif

end subroutine mInvert


end module mod_math
