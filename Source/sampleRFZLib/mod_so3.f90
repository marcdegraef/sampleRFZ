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

module mod_so3
  !! author: MDG
  !! version: 1.0
  !! date: 01/21/20
  !!
  !! everything that has to do with sampling of rotation space SO(3)

use mod_kinds
use mod_global
use mod_rotations

IMPLICIT NONE
private

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! the following table is used for two-phase disorientation fundamental zones
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! this table encodes Figure 1 of the paper  "Representation of Orientation and
! Disorientation data for Cubic, Hexagonal, Tetragonal, and Orthorhombic Crystals", A. Heinz
! and P. Neumann, Acta Cryst. A47, 780-789 (1991)
! The following conversions are used
! 0 -> x  (no symmetry)
! 1 -> a  mixed cubic-hexagonal FZ
! 2 -> b  mixed FZ
! 3 -> c  octahedral FZ
! 4 -> d  tetrahedral FZ
! 5 -> e  24-sided prismatic FZ
! 6 -> f  622 hexagonal dihedral FZ
! 7 -> g  422 octagonal dihedral FZ
! 8 -> h  32 trigonal dihedral FZ
! 9 -> i  222 dihedral FZ
! This table is used in the so3.f90 module to figure out which FZ should be used for a single phase
! or two phase FZ computation; all FZs are also available in the povray.f90 module for 3D visualization.
! The new routine getFZtypeandorder in so3.f90 will take two point group numbers, possibly identical,
! and return the class FZtype and FZorder parameters that are currently used already in other routines.
integer(kind=irg), dimension(32,32) :: FZtypeTable = reshape( (/ &
 0, 0, 0, 0, 0, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 4, 4, 3, 4, 3, &
 0, 0, 0, 0, 0, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 4, 4, 3, 4, 3, &
 0, 0, 0, 0, 0, 9, 9, 9, 7, 9, 7, 7, 7, 7, 7, 8, 8, 6, 8, 6, 6, 8, 6, 6, 6, 6, 6, 4, 4, 3, 4, 3, &
 0, 0, 0, 0, 0, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 4, 4, 3, 4, 3, &
 0, 0, 0, 0, 0, 9, 9, 9, 7, 9, 7, 7, 7, 7, 7, 8, 8, 6, 8, 6, 6, 8, 6, 6, 6, 6, 6, 4, 4, 3, 4, 3, &
 9, 9, 9, 9, 9, 9, 9, 9, 7, 9, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 3, 4, 3, &
 0, 0, 9, 0, 9, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 4, 4, 3, 4, 3, &
 9, 9, 9, 9, 9, 9, 9, 9, 7, 9, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 3, 4, 3, &
 0, 0, 7, 0, 7, 7, 0, 7, 0, 0, 0, 7, 0, 7, 7, 0, 0, 5, 0, 5, 0, 0, 0, 5, 0, 5, 5, 3, 3, 3, 3, 3, &
 0, 0, 9, 0, 9, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 4, 4, 3, 4, 3, &
 0, 0, 7, 0, 7, 7, 0, 7, 0, 0, 0, 7, 0, 7, 7, 0, 0, 5, 0, 5, 0, 0, 0, 5, 0, 5, 5, 3, 3, 3, 3, 3, &
 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3, &
 0, 0, 7, 0, 7, 7, 0, 7, 0, 0, 0, 7, 0, 7, 7, 0, 0, 5, 0, 5, 0, 0, 0, 5, 0, 5, 5, 3, 3, 3, 3, 3, &
 9, 9, 7, 9, 7, 7, 9, 7, 7, 9, 7, 7, 7, 9, 7, 6, 6, 5, 6, 5, 6, 6, 6, 5, 6, 5, 5, 3, 3, 3, 3, 3, &
 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3, &
 0, 0, 8, 0, 8, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 2, 2, 1, 2, 1, &
 0, 0, 8, 0, 8, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 2, 2, 1, 2, 1, &
 8, 8, 6, 8, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 6, 8, 6, 6, 6, 8, 6, 2, 2, 1, 2, 1, &
 0, 0, 8, 0, 8, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 2, 2, 1, 2, 1, &
 8, 8, 6, 8, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 6, 8, 6, 6, 6, 8, 6, 2, 2, 1, 2, 1, &
 0, 0, 6, 0, 6, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 2, 2, 1, 2, 1, &
 0, 0, 8, 0, 8, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 6, 6, 2, 2, 1, 2, 1, &
 0, 0, 6, 0, 6, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 2, 2, 1, 2, 1, &
 6, 6, 6, 6, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 1, 2, 1, &
 0, 0, 6, 0, 6, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 2, 2, 1, 2, 1, &
 8, 8, 6, 8, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 6, 6, 6, 6, 6, 8, 6, 2, 2, 1, 2, 1, &
 6, 6, 6, 6, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 1, 2, 1, &
 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 3, 4, 3, &
 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 3, 4, 3, &
 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, &
 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 3, 4, 3, &
 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3 &
 /), (/ 32, 32/) )

! The following two arrays are used to determine the FZtype (FZtarray) and primary rotation axis order (FZoarray)
! for each of the 32 crystallographic point group symmetries (in the order of the International Tables)
!
!                                       '    1','   -1','    2','    m','  2/m','  222', &
!                                       '  mm2','  mmm','    4','   -4','  4/m','  422', &
!                                       '  4mm',' -42m','4/mmm','    3','   -3','   32', &
!                                       '   3m','  -3m','    6','   -6','  6/m','  622', &
!                                       '  6mm',' -6m2','6/mmm','   23','   m3','  432', &
!                                       ' -43m',' m-3m'/
!
! 1 (C1), -1 (Ci), [triclinic]
! 2 (C2), m (Cs), 2/m (C2h), [monoclinic]
! 222 (D2), mm2 (C2v), mmm (D2h), [orthorhombic]
! 4 (C4), -4 (S4), 4/m (C4h), 422 (D4), 4mm (C4v), -42m (D2d), 4/mmm (D4h), [tetragonal]
! 3 (C3), -3 (C3i), 32 (D3), 3m (C3v), -3m (D3d), [trigonal]
! 6 (C6), -6 (C3h), 6/m (C6h), 622 (D6), 6mm (C6v), -6m2 (D3h), 6/mmm (D6h), [hexagonal]
! 23 (T), m3 (Th), 432 (O), -43m (Td), m-3m (Oh) [cubic]
!
! FZtype
! 0        no symmetry at all
! 1        cyclic symmetry
! 2        dihedral symmetry
! 3        tetrahedral symmetry
! 4        octahedral symmetry
!
integer(kind=irg),dimension(36)     :: FZtarray = (/ 0,0,1,1,1,2,2,2,1,1,1,2,2,2,2,1,1,2, &
                                                     2,2,1,1,1,2,2,2,2,3,3,4,3,4,5,2,2,2 /)

integer(kind=irg),dimension(36)     :: FZoarray = (/ 0,0,2,2,2,2,2,2,4,4,4,4,4,4,4,3,3,3, &
                                                     3,3,6,6,6,6,6,6,6,0,0,0,0,0,0,8,10,12 /)



! public :: SampleRFZ, IsinsideFZ, CubochoricNeighbors

! logical functions to determine if point is inside specific FZ
!private :: insideCyclicFZ, insideDihedralFZ, insideCubicFZ
! public:: insideCyclicFZ, insideDihedralFZ, insideCubicFZ

type, public :: FZpointd
  type(r_T)               :: rod       ! Rodrigues-Frank vector [nx, ny, nz, tan(omega/2) ]
  type(r_T)               :: trod      ! second Rodrigues-Frank vector; can be used for coordinate transformations
  integer(kind=irg)       :: gridpt(3) ! coordinates of grid point ! added on 06/19/18 by SS
  type(FZpointd),pointer  :: next      ! link to next point
end type FZpointd


type, public :: so3_T
  private
    integer(kind=irg)       :: FZtype
    integer(kind=irg)       :: FZ2type
    integer(kind=irg)       :: FZorder
    integer(kind=irg)       :: pgnum
    integer(kind=irg)       :: gridtype
    integer(kind=irg)       :: FZcnt
    type(FZpointd),pointer  :: FZlist
  contains
  private

    procedure, pass(self) :: getFZtypeandorder_
    procedure, pass(self) :: setFZtypeandorder_
    procedure, pass(self) :: IsinsideFZ_
    procedure, pass(self) :: insideCyclicFZ_
    procedure, pass(self) :: insideDihedralFZ_
    procedure, pass(self) :: insideCubicFZ_
    procedure, pass(self) :: getListHead_
    procedure, pass(self) :: getListCount_
    procedure, pass(self) :: setGridType_

    procedure, pass(self) :: delete_FZlist_
    procedure, pass(self) :: nullifyList_
    procedure, pass(self) :: SampleRFZ_
    procedure, pass(self) :: writeOrientationstoFile_
! some other related routines
    final :: so3_destructor

    generic, public :: getFZtypeandorder => getFZtypeandorder_
    generic, public :: setFZtypeandorder => setFZtypeandorder_
    generic, public :: IsinsideFZ => IsinsideFZ_
    generic, public :: insideCyclicFZ => insideCyclicFZ_
    generic, public :: insideDihedralFZ => insideDihedralFZ_
    generic, public :: insideCubicFZ => insideCubicFZ_
    generic, public :: getListHead => getListHead_
    generic, public :: getListCount => getListCount_
    generic, public :: setGridType => setGridType_

    generic, public :: delete_FZlist => delete_FZlist_
    generic, public :: nullifyList => nullifyList_
    generic, public :: SampleRFZ => SampleRFZ_
    generic, public :: writeOrientationstoFile => writeOrientationstoFile_

end type so3_T

! the constructor routine for this class
interface so3_T
  module procedure so3_constructor
end interface so3_T

contains

!--------------------------------------------------------------------------
type(so3_T) function so3_constructor( pgnum) result(SO)
!DEC$ ATTRIBUTES DLLEXPORT :: so3_constructor
!! author: MDG
!! version: 1.0
!! date: 01/21/20
!!
!! constructor for the so3_T Class

IMPLICIT NONE

integer(kind=irg), INTENT(IN)             :: pgnum
 !! primary point group

call SO%setFZtypeandorder(pgnum)
SO%pgnum = pgnum

call SO%nullifyList()

end function so3_constructor

!--------------------------------------------------------------------------
subroutine so3_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: so3_destructor
!! author: MDG
!! version: 1.0
!! date: 02/02/20
!!
!! destructor for the so3_T Class

IMPLICIT NONE

type(so3_T), INTENT(INOUT)  :: self

call reportDestructor('so3_T')

! nothing to do for now...

end subroutine so3_destructor

!--------------------------------------------------------------------------
subroutine nullifyList_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: nullifyList_
!! author: MDG
!! version: 1.0
!! date: 02/18/20
!!
!! nullify the selected list

IMPLICIT NONE

class(so3_T), INTENT(INOUT)           :: self

nullify(self%FZlist)
self%FZcnt = 0

end subroutine nullifyList_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Routine to return the FZtype and FZorder parameters for single or two-phase
! fundamental zone (FZ) computations; this includes all the FZ types from the
! following paper:
!
! "Representation of Orientation and Disorientation data for Cubic, Hexagonal,
! Tetragonal, and Orthorhombic Crystals", A. Heinz and P. Neumann, Acta Cryst. A47,
! 780-789 (1991)
!
! this routine also allows for icosahedral symmetry, although this is not part
! of the paper above.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine setFZtypeandorder_(self, pgnum)
!DEC$ ATTRIBUTES DLLEXPORT :: setFZtypeandorder_
!! author: MDG
!! version: 1.0
!! date: 01/21/20
!!
!! set the point group number(s) and the fundamental zone type and order

IMPLICIT NONE

class(so3_T),INTENT(INOUT)                :: self
integer(kind=irg),INTENT(IN)              :: pgnum

self%FZtype = FZtarray(pgnum)
self%FZorder = FZoarray(pgnum)

end subroutine setFZtypeandorder_

!--------------------------------------------------------------------------
recursive subroutine getFZtypeandorder_(self, FZtype, FZorder)
!DEC$ ATTRIBUTES DLLEXPORT :: getFZtypeandorder_
!! author: MDG
!! version: 1.0
!! date: 01/21/20
!!
!! set the point group number(s) and the fundamental zone type and order

IMPLICIT NONE

class(so3_T),INTENT(INOUT)                :: self
integer(kind=irg), INTENT(OUT)            :: FZtype
integer(kind=irg), INTENT(OUT)            :: FZorder

FZtype = self%FZtype
FZorder = self%FZorder

end subroutine getFZtypeandorder_


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! We define a number of logical routines, that decide whether or not
! a point in Rodrigues representation lies inside the fundamental zone (FZ)
! for a given crystal symmetry. This follows the Morawiec@Field paper:
!
! A. Morawiec & D. P. Field (1996) Rodrigues parameterization for orientation
! and misorientation distributions, Philosophical Magazine A, 73:4, 1113-1130,
! DOI: 10.1080/01418619608243708
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive function IsinsideFZ_(self, rod, qFZ) result(insideFZ)
!DEC$ ATTRIBUTES DLLEXPORT :: IsinsideFZ_
  !! author: MDG
  !! version: 1.0
  !! date: 01/21/20
  !!
  !! does Rodrigues point lie inside the relevant FZ?

use mod_math
use mod_quaternions

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self
type(r_T), INTENT(INOUT)                :: rod
 !! input Rodrigues vector
type(q_T), INTENT(INOUT), OPTIONAL      :: qFZ
 !! quaternion that rotates the FZ into a new orientation (optional)

logical                                 :: insideFZ

type(r_T)                               :: newrod
type(q_T)                               :: qu
type(quaternion_T)                      :: qu1, qu2, qq
real(kind=dbl)                          :: x(4)

! do we need to rotate the FZ ? (we do this by rotating rod in the opposite way)
if (present(qFZ)) then
  qu = rod%rq()
  qu1 = quaternion_T( qd = qu%q_copyd() )
  qu2 = quaternion_T( qd = qFZ%q_copyd() )
  qq = qu2 * ( qu1 * conjg(qu2) )
  qu = q_T( qdinp = qq%get_quatd() )
  newrod = qu%qr()
else
  newrod = rod
end if

insideFZ = .FALSE.

! dealing with 180 rotations is needed only for
! FZtypes 0 and 1; the other FZs are always finite.
x = newrod%r_copyd()
  select case (self%FZtype)
    case (0)
      insideFZ = .TRUE.   ! all points are inside the FZ
    case (1)
      insideFZ = self%insideCyclicFZ(newrod)        ! infinity is checked inside this function
    case (2)
      if (x(4).ne.inftyd()) insideFZ = self%insideDihedralFZ(newrod, self%FZorder)
    case (3)
      if (x(4).ne.inftyd()) insideFZ = self%insideCubicFZ(newrod,'tet')
    case (4)
      if (x(4).ne.inftyd()) insideFZ = self%insideCubicFZ(newrod,'oct')
  end select

end function IsinsideFZ_

!--------------------------------------------------------------------------
recursive function insideCyclicFZ_(self, rod, M) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: insideCyclicFZ_
  !! author: MDG
  !! version: 1.0
  !! date: 01/21/20
  !!
  !! does Rodrigues point lie inside cyclic FZ (for 2, 3, 4, and 6-fold)?

use mod_math

IMPLICIT NONE

class(so3_T),INTENT(INOUT)        :: self
type(r_T), INTENT(INOUT)          :: rod
logical,INTENT(IN),OPTIONAL       :: M

logical                           :: res, doM
real(kind=dbl)                    :: x(4)

res = .FALSE.
doM = .FALSE.
if (present(M)) then
  if (M.eqv..TRUE.) doM = .TRUE.
end if

x = rod%r_copyd()
if (x(4).ne.inftyd()) then
    if ((self%FZtype.eq.1.).and.(self%FZorder.eq.2)) then
! check the y-component vs. tan(pi/2n)
      res = dabs(x(2)*x(4)).le.LPs%BP(self%FZorder)
    else
! check the z-component vs. tan(pi/2n)
      res = dabs(x(3)*x(4)).le.LPs%BP(self%FZorder)
    end if
else
    if ((self%FZtype.eq.1.).and.(self%FZorder.eq.2)) then
      if(x(2) .eq. 0.D0) res = .TRUE.
    else
      if (x(3).eq.0.D0) res = .TRUE.
    end if
endif

end function insideCyclicFZ_

!--------------------------------------------------------------------------
recursive function insideDihedralFZ_(self, rod, order) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: insideDihedralFZ_
  !! author: MDG
  !! version: 1.0
  !! date: 01/21/20
  !!
  !! does Rodrigues point lie inside cyclic FZ (for 2, 3, 4, and 6-fold)?

IMPLICIT NONE

class(so3_T),INTENT(INOUT)        :: self
type(r_T), INTENT(INOUT)          :: rod
integer(kind=irg), INTENT(IN)     :: order

logical                           :: res, c1, c2
real(kind=dbl)                    :: r(3), x(4)
real(kind=dbl),parameter          :: r1 = 1.00D0
real(kind=dbl),allocatable        :: polygonvertex(:,:)
integer(kind=irg)                 :: inout

x = rod%r_copyd()
if (x(4).gt.sqrt(3.D0)) then
  res = .FALSE.
else
  r(1:3) = x(1:3) * x(4)

  ! first, check the z-component vs. tan(pi/2n)  (same as insideCyclicFZ)
  c1 = dabs(r(3)).le.LPs%BP(order)
  res = .FALSE.

  ! check the square boundary planes if c1=.TRUE.
  if (c1) then
    select case (order)
      case (2)
        c2 = (dabs(r(1)).le.r1).and.(dabs(r(2)).le.r1)
      case (3)
        c2 =          dabs( LPs%srt*r(1)+0.5D0*r(2)).le.r1
        c2 = c2.and.( dabs( LPs%srt*r(1)-0.5D0*r(2)).le.r1 )
        c2 = c2.and.( dabs(r(2)).le.r1 )
      case (4)
        c2 = (dabs(r(1)).le.r1).and.(dabs(r(2)).le.r1)
        c2 = c2.and.((LPs%r22*dabs(r(1)+r(2)).le.r1).and.(LPs%r22*dabs(r(1)-r(2)).le.r1))
      case (6)
        c2 =          dabs( 0.5D0*r(1)+LPs%srt*r(2)).le.r1
        c2 = c2.and.( dabs( LPs%srt*r(1)+0.5D0*r(2)).le.r1 )
        c2 = c2.and.( dabs( LPs%srt*r(1)-0.5D0*r(2)).le.r1 )
        c2 = c2.and.( dabs( 0.5D0*r(1)-LPs%srt*r(2)).le.r1 )
        c2 = c2.and.( dabs(r(2)).le.r1 )
        c2 = c2.and.( dabs(r(1)).le.r1 )
    end select
    res = c2
  end if
end if

end function insideDihedralFZ_

!--------------------------------------------------------------------------
recursive function insideCubicFZ_(self, rod, ot) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: insideCubicFZ_
  !! author: MDG
  !! version: 1.0
  !! date: 01/21/20
  !!
  !! does Rodrigues point lie inside cubic FZ (octahedral or tetrahedral)?

IMPLICIT NONE

class(so3_T),INTENT(INOUT)        :: self
type(r_T), INTENT(INOUT)          :: rod
character(3), INTENT(IN)          :: ot

logical                           :: res, c1, c2
real(kind=dbl)                    :: r(3), x(4)
real(kind=dbl),parameter          :: r1  = 1.0D0
real(kind=dbl),parameter          :: eps = 1.0D-8

x = rod%r_copyd()
r(1:3) = x(1:3) * x(4)

res = .FALSE.

! primary cube planes (only needed for octahedral case)
if (ot.eq.'oct') then
  c1 = (maxval(dabs(r)) - LPS%BP(4) .le. eps)
else
  c1 = .TRUE.
end if

! octahedral truncation planes, both for tetrahedral and octahedral point groups
c2 = ((dabs(r(1))+dabs(r(2))+dabs(r(3))) - r1 .le. eps)

! if both c1 and c2, then the point is inside
if (c1.and.c2) res = .TRUE.

end function insideCubicFZ_

!--------------------------------------------------------------------------
recursive subroutine delete_FZlist_(self, l)
!DEC$ ATTRIBUTES DLLEXPORT :: delete_FZlist_
  !! author: MDG
  !! version: 1.0
  !! date: 01/21/20
  !!
  !! delete a linked list of rodrigues vectors

class(so3_T),INTENT(INOUT)        :: self
character(2), INTENT(IN),OPTIONAL :: l

type(FZpointd),pointer            :: ltail, ltmp

ltail => self%FZlist
self%FZcnt = 0

! deallocate the entire linked list before returning, to prevent memory leaks
ltmp => ltail % next
do
  if (associated(ltail)) deallocate(ltail)
  if (.not. associated(ltmp)) EXIT
  ltail => ltmp
  ltmp => ltail % next
end do

end subroutine delete_FZlist_

!--------------------------------------------------------------------------
recursive subroutine SampleRFZ_(self, nsteps, qFZ)
!DEC$ ATTRIBUTES DLLEXPORT :: SampleRFZ_
  !! author: MDG
  !! version: 1.0
  !! date: 01/21/20
  !!
  !! Generate a uniform sampling of a Rodriguess FZ
  !!
  !! This routine fills in a linked list FZlist of Rodrigues points that
  !! are inside a specific fundamental zone determined by the sample point group;
  !! this list can then be further dealt with in the calling program.
  !!
  !! Here's how you would use this routine in a main program:
  !!
  !!    use mod_so3
  !!    use mod_rotations
  !!
  !!    IMPLICIT NONE
  !!
  !!    type(so3_T)             :: SO
  !!    integer(kind=irg)       :: i, pgnum, nsteps, FZcnt
  !!    type(FZpointd), pointer :: FZtmp
  !!    type(e_T)               :: eu
  !!
  !!    pgnum = 32
  !!    SO = so3_T( pgnum )
  !!
  !!    nsteps = 10
  !!    call SO%sampleRFZ(nsteps)
  !!
  !! Then you can access all the entries in the list and, for instance, convert them to Euler angles...
  !!
  !!    FZtmp => SO%getListHead('FZ')          ! point to the top of the list
  !!    FZcnt = SO%getListCount('FZ')          ! get the number of entries in the list
  !!    do i = 1, FZcnt                        ! loop over all entries
  !!      eu = FZtmp%rod%re()                  ! convert to Euler angles (in radians by default)
  !!    !  do something with eu                ! for instance, write eu to a file
  !!      FZtmp => FZtmp%next                  ! point to the next entry
  !!    end do
  !!
  !! If you just want to look at the first 10 entries on the list and show all other orientation representations:
  !!
  !!    type(orientation_T) :: ot
  !!
  !!    FZtmp => SO%getListHead('FZ')
  !!    do i = 1,10
  !!      ot = orientation_T( FZtmp%rod )
  !!      call ot%print_orientation('d')    ! the argument 'd' means angles will be in degrees ('r' for radians)
  !!      FZtmp => FZtmp%next
  !!    end do

IMPLICIT NONE

class(so3_T),INTENT(INOUT)           :: self

integer(kind=irg), INTENT(IN)        :: nsteps
type(q_T),INTENT(INOUT),OPTIONAL     :: qFZ

type(r_T)                            :: rod
type(c_T)                            :: cu
real(kind=dbl)                       :: x, y, z, delta, shift, sedge, ztmp
type(FZpointd), pointer              :: FZtmp, FZtmp2
integer(kind=irg)                    :: i, j, k
logical                              :: b, rotateFZ = .FALSE.

if (present(qFZ)) rotateFZ = .TRUE.

! cube semi-edge length s = 0.5D0 * LPs%ap
! step size for sampling of grid; total number of samples = (2*nsteps+1)**3
sedge = 0.5D0 * LPs%ap
delta = sedge / dble(nsteps)

if (self%gridtype.eq.0) then
  shift = 0.0D0
else
  shift = 0.5D0
end if

! set the counter to zero
self%FZcnt = 0

! note that when FZtype is cyclic (1) and FZorder is 2, then we must rotate the
! rotation axis to lie along the b (y) direction, not z !!!!

! loop over the cube of volume pi^2; note that we do not want to include
! the opposite edges/facets of the cube, to avoid double counting rotations
! with a rotation angle of 180 degrees.  This only affects the cyclic groups.

 do i=-nsteps+1,nsteps
  x = (dble(i)+shift)*delta
  do j=-nsteps+1,nsteps
   y = (dble(j)+shift)*delta
   do k=-nsteps+1,nsteps
    z = (dble(k)+shift)*delta
! make sure that this point lies inside the cubochoric cell
    if (maxval( (/ abs(x), abs(y), abs(z) /) ).le.sedge) then

! convert to Rodrigues representation
      cu = c_T( cdinp = (/ x, y, z /) )
      rod = cu%cr()

! If insideFZ=.TRUE., then add this point to the linked list FZlist and keep
! track of how many points there are on this list
       if (rotateFZ.eqv..TRUE.) then
         b = self%IsinsideFZ(rod, qFZ)
       else
         b = self%IsinsideFZ(rod)
       end if
       if (b) then
        if (.not.associated(self%FZlist)) then
          allocate(self%FZlist)
          FZtmp => self%FZlist
        else
          allocate(FZtmp%next)
          FZtmp => FZtmp%next
        end if
        nullify(FZtmp%next)
! if monoclinic, then reorder the components !!!
!        if ((FZtype.eq.1).and.(FZorder.eq.2)) then
!          ztmp = rod(3)
!          rod(3) = rod(1)
!          rod(1) = rod(2)
!          rod(2) = ztmp
!        end if
        FZtmp%rod = rod
        FZtmp%gridpt(1:3) = (/i, j, k/)
        self%FZcnt = self%FZcnt + 1
       end if
    end if
  end do
 end do
end do

end subroutine SampleRFZ_

!--------------------------------------------------------------------------
recursive subroutine writeOrientationstoFile_(self, filename, mode, list, trod)
!DEC$ ATTRIBUTES DLLEXPORT :: writeOrientationstoFile_
  !! author: MDG
  !! version: 1.0
  !! date: 01/22/20
  !!
  !! write a list of orientations from a linked list to a text file

use mod_io
use mod_math

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self

character(fnlen),INTENT(IN)             :: filename
 !! complete path to output file name
character(2), INTENT(IN)                :: mode
 !! output orientation representation  (eu, ro, ho, ...)
character(2), INTENT(IN), OPTIONAL      :: list
 !! list from which to write
logical, INTENT(IN), OPTIONAL           :: trod
 !! list from which to write

type(e_T)                               :: e
type(o_T)                               :: o
type(q_T)                               :: q
type(s_T)                               :: s
type(v_T)                               :: v
type(h_T)                               :: h
type(c_T)                               :: c
type(r_T)                               :: r
type(a_T)                               :: a

type(IO_T)                              :: Message
type(FZpointd), pointer                 :: FZtmp
integer(kind=irg)                       :: cnt, i
real(kind=dbl)                          :: io_real(9)
logical                                 :: dotrod

dotrod = .FALSE.
if (present(trod)) then
  if (trod.eqv..TRUE.) dotrod = .TRUE.
end if


FZtmp => self%FZlist
cnt = self%FZcnt

open(unit=53, file=trim(filename), status='unknown', form='formatted')
write (53,"(A2)") mode
write (53,"(I6)") cnt

if (dotrod.eqv..TRUE.) then
  do i=1, cnt
    select case(mode)
      case('eu')
        e = FZtmp%trod%re()
        io_real(1:3) = e%e_copyd() / dtor
        call Message%WriteValue('', io_real, 3, frm="(2(F17.9,' '),F17.9)",redirect=53)
      case('ro')
        io_real(1:4) = FZtmp%trod%r_copyd()
        if (io_real(4).eq.inftyd()) then
          call Message%WriteValue('', io_real, 3, frm="(2(F17.9,' '),'infinity')",redirect=53)
        else
          call Message%WriteValue('', io_real, 4, frm="(2(F17.9,' '),F17.9)",redirect=53)
        end if
      case('om')
        o = FZtmp%trod%ro()
        io_real(1:9) = reshape(o%o_copyd(), (/ 9 /) )
        call Message%WriteValue('', io_real, 9, frm="(8(F17.9,' '),F17.9)",redirect=53)
      case('ho')
        h = FZtmp%trod%rh()
        io_real(1:3) = h%h_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F17.9,' '),F17.9)",redirect=53)
      case('cu')
        c = FZtmp%trod%rc()
        io_real(1:3) = c%c_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F17.9,' '),F17.9)",redirect=53)
      case('rv')
        v = FZtmp%trod%rv()
        io_real(1:3) = v%v_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F17.9,' '),F17.9)",redirect=53)
      case('st')
        s = FZtmp%trod%rs()
        io_real(1:3) = s%s_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F17.9,' '),F17.9)",redirect=53)
      case('ax')
        a = FZtmp%trod%ra()
        io_real(1:4) = a%a_copyd()
        io_real(4) = io_real(4) / dtor
        call Message%WriteValue('', io_real, 4, frm="(3(F17.9,' '),F17.9)",redirect=53)
      case('qu')
        q = FZtmp%trod%rq()
        io_real(1:4) = q%q_copyd()
        call Message%WriteValue('', io_real, 4, frm="(3(F17.9,' '),F17.9)",redirect=53)
      case default
    end select
    FZtmp => FZtmp%next
  end do
else
  do i=1, cnt
    select case(mode)
      case('eu')
        e = FZtmp%rod%re()
        io_real(1:3) = e%e_copyd() / dtor
        call Message%WriteValue('', io_real, 3, frm="(2(F17.9,' '),F17.9)",redirect=53)
      case('ro')
        io_real(1:4) = FZtmp%rod%r_copyd()
        if (io_real(4).eq.inftyd()) then
          call Message%WriteValue('', io_real, 3, frm="(3(F17.9,' '),'infinity')",redirect=53)
        else
          call Message%WriteValue('', io_real, 4, frm="(3(F17.9,' '),F17.9)",redirect=53)
        end if
      case('om')
        o = FZtmp%rod%ro()
        io_real(1:9) = reshape(o%o_copyd(), (/ 9 /) )
        call Message%WriteValue('', io_real, 9, frm="(8(F17.9,' '),F17.9)",redirect=53)
      case('ho')
        h = FZtmp%rod%rh()
        io_real(1:3) = h%h_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F17.9,' '),F17.9)",redirect=53)
      case('cu')
        c = FZtmp%rod%rc()
        io_real(1:3) = c%c_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F17.9,' '),F17.9)",redirect=53)
      case('rv')
        v = FZtmp%rod%rv()
        io_real(1:3) = v%v_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F17.9,' '),F17.9)",redirect=53)
      case('st')
        s = FZtmp%rod%rs()
        io_real(1:3) = s%s_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F17.9,' '),F17.9)",redirect=53)
      case('ax')
        a = FZtmp%rod%ra()
        io_real(1:4) = a%a_copyd()
        io_real(4) = io_real(4) / dtor
        call Message%WriteValue('', io_real, 4, frm="(3(F17.9,' '),F17.9)",redirect=53)
      case('qu')
        q = FZtmp%rod%rq()
        io_real(1:4) = q%q_copyd()
        call Message%WriteValue('', io_real, 4, frm="(3(F17.9,' '),F17.9)",redirect=53)
      case default
    end select
    FZtmp => FZtmp%next
  end do
end if

close(unit=53, status = 'keep')

end subroutine writeOrientationstoFile_


!--------------------------------------------------------------------------
recursive function getListHead_(self, l) result(FZptr)
!DEC$ ATTRIBUTES DLLEXPORT :: getListHead_
  !! author: MDG
  !! version: 1.0
  !! date: 01/22/20
  !!
  !! return the pointer to the selected list

IMPLICIT NONE

class(so3_T),INTENT(INOUT)    :: self
character(2), INTENT(IN)      :: l

type(FZpointd), pointer       :: FZptr

FZptr => self%FZlist

end function getListHead_

!--------------------------------------------------------------------------
recursive function getListCount_(self, l) result(cnt)
!DEC$ ATTRIBUTES DLLEXPORT :: getListCount_
  !! author: MDG
  !! version: 1.0
  !! date: 01/22/20
  !!
  !! return the counter for the selected list

IMPLICIT NONE

class(so3_T),INTENT(INOUT)    :: self
character(2), INTENT(IN)      :: l

integer(kind=irg)             :: cnt

cnt = self%FZcnt

end function getListCount_

!--------------------------------------------------------------------------
recursive subroutine setGridType_(self, g)
!DEC$ ATTRIBUTES DLLEXPORT :: setGridType_
  !! author: MDG
  !! version: 1.0
  !! date: 01/22/20
  !!
  !! set the gridtype parameter

IMPLICIT NONE

class(so3_T),INTENT(INOUT)    :: self
integer(kind=irg)             :: g

self%gridtype = g

end subroutine setGridType_

end module mod_so3
