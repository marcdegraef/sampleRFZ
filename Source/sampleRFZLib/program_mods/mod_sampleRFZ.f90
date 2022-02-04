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

module mod_sampleRFZ
  !! author: MDG
  !! version: 1.0
  !! date: 01/22/20
  !!
  !! class definition for the EMsampleRFZ program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMsampleRFZ program
type, public :: sampleRFZNameListType
    integer(kind=irg) :: pgnum
    integer(kind=irg) :: nsteps
    integer(kind=irg) :: gridtype
    character(fnlen)  :: xtalname
    character(fnlen)  :: euoutname
    character(fnlen)  :: cuoutname
    character(fnlen)  :: hooutname
    character(fnlen)  :: rooutname
    character(fnlen)  :: quoutname
    character(fnlen)  :: omoutname
    character(fnlen)  :: axoutname
    character(fnlen)  :: rvoutname
    character(fnlen)  :: stoutname
end type sampleRFZNameListType

type, public :: sampleRFZ_T
private
  character(fnlen)             :: nmldeffile = 'EMsampleRFZ.nml'
  type(sampleRFZNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: CreateSampling_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: CreateSampling => CreateSampling_

end type sampleRFZ_T

! the constructor routine for this class
interface sampleRFZ_T
  module procedure sampleRFZ_constructor
end interface sampleRFZ_T

contains

!--------------------------------------------------------------------------
type(sampleRFZ_T) function sampleRFZ_constructor( nmlfile ) result(RFZ)
!DEC$ ATTRIBUTES DLLEXPORT :: sampleRFZ_constructor
!! author: MDG
!! version: 1.0
!! date: 01/22/20
!!
!! constructor for the sampleRFZ_T Class; reads the name list file

IMPLICIT NONE

character(fnlen), INTENT(IN)    :: nmlfile

call RFZ%readNameList(nmlfile)

end function sampleRFZ_constructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 01/22/20
!!
!! read the namelist for the sampleRFZ_T Class

use mod_io

IMPLICIT NONE

class(sampleRFZ_T), INTENT(INOUT)  :: self
character(fnlen),INTENT(IN)        :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)        :: initonly
 !! fill in the default values only; do not read the file

logical                            :: skipread = .FALSE.

integer(kind=irg)                  :: pgnum, nsteps, gridtype
character(fnlen)                   :: euoutname
character(fnlen)                   :: cuoutname
character(fnlen)                   :: hooutname
character(fnlen)                   :: rooutname
character(fnlen)                   :: quoutname
character(fnlen)                   :: omoutname
character(fnlen)                   :: axoutname
character(fnlen)                   :: rvoutname
character(fnlen)                   :: stoutname

! namelist components
namelist / RFZlist / pgnum, nsteps, gridtype, euoutname, cuoutname, hooutname, rooutname, quoutname, omoutname, axoutname, &
                     rvoutname, stoutname

! initialize to default values
pgnum = 32
nsteps = 50
gridtype = 0
euoutname = 'undefined'
cuoutname = 'undefined'
hooutname = 'undefined'
rooutname = 'undefined'
quoutname = 'undefined'
omoutname = 'undefined'
axoutname = 'undefined'
rvoutname = 'undefined'
stoutname = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=RFZlist)
close(UNIT=dataunit,STATUS='keep')
end if

! and copy the variables to the namelist variable
self%nml%pgnum  = pgnum
self%nml%nsteps = nsteps
self%nml%gridtype = gridtype
self%nml%euoutname = euoutname
self%nml%cuoutname = cuoutname
self%nml%hooutname = hooutname
self%nml%rooutname = rooutname
self%nml%quoutname = quoutname
self%nml%omoutname = omoutname
self%nml%axoutname = axoutname
self%nml%rvoutname = rvoutname
self%nml%stoutname = stoutname

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 01/22/20
!!
!! pass the namelist for the sampleRFZ_T Class to the calling program

IMPLICIT NONE

class(sampleRFZ_T), INTENT(INOUT)    :: self
type(sampleRFZNameListType)          :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
!
! SUBROUTINE:CreateSampling
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Generate a sampling of the Rodrigues Fundamental Zone for a given xtal symmetry
!
!> @todo add an HDF5 output option
!
!> @param nmlfile namelist file name
!
!> @date 05/29/14 MDG 1.0 original
!> @date 12/09/14 MDG 2.0 changed rfznl handling
!> @date 08/19/15 MDG 2.1 added all rotation representations as output options
!> @date 12/22/16 MDG 2.2 added option to generate reduced sampling inside constant misorientation ball
!> @date 02/01/17 MDG 2.3 added option to generate sampling inside a conical volume in Rodrigues space
!> @date 08/16/17 MDG 2.4 added option to generate uniform fiber texture sampling in Rodrigues space
!--------------------------------------------------------------------------
subroutine CreateSampling_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: CreateSampling_
!! author: MDG
!! version: 1.0
!! date: 01/22/20
!!
!! Generate a sampling of the Rodrigues Fundamental Zone for a given xtal symmetry

use mod_kinds
use mod_global
use mod_symmetry
use mod_io
use mod_so3

IMPLICIT NONE

class(sampleRFZ_T), INTENT(INOUT)  :: self

type(sampleRFZNameListType)        :: rfznl
type(IO_T)                         :: Message
type(so3_T)                        :: SO

integer(kind=irg)                  :: i, j, num, m, io_int(1), FZcnt, FZtype, FZorder
real(kind=dbl)                     :: x, h, k, l, ih, ik, il, idiff, eps = 0.0001D0
real(kind=dbl),allocatable         :: itmp(:,:)
logical                            :: doeu = .FALSE., docu = .FALSE., doho = .FALSE., doqu = .FALSE., dorv = .FALSE., &
                                      dost = .FALSE., doom = .FALSE., doax = .FALSE., doro = .FALSE., newpoint, &
                                      rotateFZ = .FALSE.
character(fnlen)                   :: filename


! first get the name list
rfznl = self%getNameList()

! determine which files to create
if (trim(rfznl%euoutname).ne.'undefined') doeu = .TRUE.
if (trim(rfznl%cuoutname).ne.'undefined') docu = .TRUE.
if (trim(rfznl%hooutname).ne.'undefined') doho = .TRUE.
if (trim(rfznl%quoutname).ne.'undefined') doqu = .TRUE.
if (trim(rfznl%rooutname).ne.'undefined') doro = .TRUE.
if (trim(rfznl%omoutname).ne.'undefined') doom = .TRUE.
if (trim(rfznl%axoutname).ne.'undefined') doax = .TRUE.
if (trim(rfznl%stoutname).ne.'undefined') dost = .TRUE.
if (trim(rfznl%rvoutname).ne.'undefined') dorv = .TRUE.

! a bit of output
call Message%printMessage('Starting computation for point group '//PGTHD(rfznl%pgnum))

! determine which function we should call for this point group symmetry
SO = so3_T( rfznl%pgnum )
call SO%setGridType( rfznl%gridtype )

! get the linked list for the FZ for point group symmetry pgnum for nsteps along the cubic semi-edge
call SO%SampleRFZ(rfznl%nsteps)

FZcnt = SO%getListCount('FZ')
io_int(1) = FZcnt
call Message%WriteValue('Total number of unique orientations generated = ',io_int,1,"(I10)")

! generate a list of all orientations in Euler angle format (if requested)
if (doeu) then
  filename = rfznl%euoutname
  call SO%writeOrientationstoFile(filename, 'eu')
end if

! generate a list of all orientations in cubochoric format (if requested)
if (docu) then
  filename = rfznl%cuoutname
  call SO%writeOrientationstoFile(filename, 'cu')
end if

! generate a list of all orientations in homochoric format (if requested)
if (doho) then
  filename = rfznl%hooutname
  call SO%writeOrientationstoFile(filename, 'ho')
end if

! generate a list of all orientations in quternion format (if requested)
if (doqu) then
  filename = rfznl%quoutname
  call SO%writeOrientationstoFile(filename, 'qu')
end if

! generate a list of all orientations in Rodrigues format (if requested)
if (doro) then
  filename = rfznl%rooutname
  call SO%writeOrientationstoFile(filename, 'ro')
end if

! generate a list of all orientations in orientation matrix format (if requested)
if (doom) then
  filename = rfznl%omoutname
  call SO%writeOrientationstoFile(filename, 'om')
end if

! generate a list of all orientations in axis angle pair format (if requested)
if (doax) then
  filename = rfznl%axoutname
  call SO%writeOrientationstoFile(filename, 'ax')
end if

! generate a list of all orientations in stereographic format (if requested)
if (dost) then
  filename = rfznl%stoutname
  call SO%writeOrientationstoFile(filename, 'st')
end if

! generate a list of all orientations in axis angle pair format (if requested)
if (dorv) then
  filename = rfznl%rvoutname
  call SO%writeOrientationstoFile(filename, 'rv')
end if

if (doeu) call Message%printMessage('Euler angles stored in file '//rfznl%euoutname)
if (docu) call Message%printMessage('Cubochoric representation stored in file '//rfznl%cuoutname)
if (doho) call Message%printMessage('Homochoric representation stored in file '//rfznl%hooutname)
if (doqu) call Message%printMessage('Quaternion representation stored in file '//rfznl%quoutname)
if (doro) call Message%printMessage('Rodrigues vector representation stored in file '//rfznl%rooutname)
if (doom) call Message%printMessage('Orientation matrix representation stored in file '//rfznl%omoutname)
if (doax) call Message%printMessage('Axis-angle pair representation stored in file '//rfznl%axoutname)
if (dost) call Message%printMessage('Stereographic representation stored in file '//rfznl%stoutname)
if (dorv) call Message%printMessage('Rotation vector representation stored in file '//rfznl%rvoutname)

end subroutine CreateSampling_


end module mod_sampleRFZ
