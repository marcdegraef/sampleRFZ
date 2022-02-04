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

module mod_quaternions
  !! author: MDG
  !! version: 1.0
  !! date: 01/03/20
  !!
  !! Quaternion and Quaternion Array arithmetic class
  !!
  !! Quaternions are defined with the scalar part in position 1, and the vector part in positions 2:4.
  !!
  !! There are two class definitions in this file, one for single quaternions, the other for
  !! quaternion array operations (some using OpenMP threads). The program MODQuaternionsTest.f90
  !! can be used as part of ctest to run a test program on this module.


use mod_global
use mod_kinds
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

IMPLICIT NONE
  private

! following are used to define the quaternion symmetry operators
    real(kind=dbl), public, parameter  :: sq22=0.7071067811865475244D0 ! sqrt(2)/2
    real(kind=dbl), public, parameter  :: sq32=0.8660254037844386467D0 ! sqrt(3)/2
    real(kind=dbl), public, parameter  :: half=0.5D0                   ! 1/2

! We define the rotational crystal symmetry operators in terms of quaternions (q0, q1,q2,q3) with q0 the scalar part;
! these are used in the dictmod EBSD dictionary indexing module, and are defined with respect to the standard cartesian
! reference frame.  Note that these are not defined as Quaternion_T type, just as arrays of 4-component doubles.
    real(kind=dbl), public, dimension(4,152) :: SYM_Qsymop = reshape( (/ &
                      1.D0, 0.D0, 0.D0, 0.D0, &       ! 1: identity operator
                      0.D0, 1.D0, 0.D0, 0.D0, &       ! 2: 180@[100]
                      0.D0, 0.D0, 1.D0, 0.D0, &       ! 3: 180@[010]
                      0.D0, 0.D0, 0.D0, 1.D0, &       ! 4: 180@[001]
                      sq22, sq22, 0.D0, 0.D0, &       ! 5: 90@[100]
                      sq22, 0.D0, sq22, 0.D0, &       ! 6: 90@[010]
                      sq22, 0.D0, 0.D0, sq22, &       ! 7: 90@[001]
                      sq22,-sq22, 0.D0, 0.D0, &       ! 8: 270@[100]
                      sq22, 0.D0,-sq22, 0.D0, &       ! 9: 270@[010]
                      sq22, 0.D0, 0.D0,-sq22, &       !10: 270@[001]
                      0.D0, sq22, sq22, 0.D0, &       !11: 180@[110]
                      0.D0,-sq22, sq22, 0.D0, &       !12: 180@[-110]
                      0.D0, 0.D0, sq22, sq22, &       !13: 180@[011]
                      0.D0, 0.D0,-sq22, sq22, &       !14: 180@[0-11]
                      0.D0, sq22, 0.D0, sq22, &       !15: 180@[101]
                      0.D0,-sq22, 0.D0, sq22, &       !16: 180@[-101]
                      half, half, half, half, &       !17: 120@[111]
                      half,-half,-half,-half, &       !18: 120@[-1-1-1]
                      half, half,-half, half, &       !19: 120@[1-11]
                      half,-half, half,-half, &       !20: 120@[-11-1]
                      half,-half, half, half, &       !21: 120@[-111]
                      half, half,-half,-half, &       !22: 120@[1-1-1]
                      half,-half,-half, half, &       !23: 120@[-1-11]
                      half, half, half,-half, &       !24: 120@[11-1]
                      sq32, 0.D0, 0.D0, half, &       !25:  60@[001] (hexagonal/trigonal operators start here)
                      half, 0.D0, 0.D0, sq32, &       !26: 120@[001]
                      0.D0, 0.D0, 0.D0, 1.D0, &       !27: 180@[001] (duplicate from above, but useful to keep it here)
                     -half, 0.D0, 0.D0, sq32, &       !28: 240@[001]
                     -sq32, 0.D0, 0.D0, half, &       !29: 300@[001]
                      0.D0, 1.D0, 0.D0, 0.D0, &       !30: 180@[100]
                      0.D0, sq32, half, 0.D0, &       !31: 180@[xxx]
                      0.D0, half, sq32, 0.D0, &       !32: 180@[xxx]
                      0.D0, 0.D0, 1.D0, 0.D0, &       !33: 180@[010]
                      0.D0,-half, sq32, 0.D0, &       !34: 180@[xxx]
                      0.D0,-sq32, half, 0.D0, &       !35: 180@[xxx]
                      1.0000000000000000D0, 0.00000000000000000D0, 0.00000000000000000D0, 0.00000000000000000D0, & ! icosahedral operators
                      0.0000000000000000D0, 0.68819093704223555D0, 0.50000000000000000D0, 0.52573108673095648D0, &
                      0.0000000000000000D0,-0.26286554336547824D0, 0.80901700258254916D0, 0.52573108673095648D0, &
                      0.0000000000000000D0,-0.85065078735351463D0, 0.00000000000000000D0, 0.52573108673095648D0, &
                      0.0000000000000000D0,-0.26286554336547824D0,-0.80901700258254916D0, 0.52573108673095648D0, &
                      0.0000000000000000D0, 0.68819093704223555D0,-0.50000000000000000D0, 0.52573108673095648D0, &
                      0.0000000000000000D0, 0.52573108673095648D0, 0.00000000000000000D0, 0.85065078735351463D0, &
                      0.0000000000000000D0, 0.16245985031127913D0, 0.50000000000000000D0, 0.85065078735351463D0, &
                      0.0000000000000000D0,-0.42532539367675731D0, 0.30901700258254972D0, 0.85065078735351463D0, &
                      0.0000000000000000D0,-0.42532539367675731D0,-0.30901700258254972D0, 0.85065078735351463D0, &
                      0.0000000000000000D0, 0.16245985031127913D0,-0.50000000000000000D0, 0.85065078735351463D0, &
                      0.0000000000000000D0, 0.95105654001235851D0,-0.30901700258254972D0, 0.00000000000000000D0, &
                      0.0000000000000000D0, 0.95105654001235851D0, 0.30901700258254972D0, 0.00000000000000000D0, &
                      0.0000000000000000D0, 0.58778524398803644D0, 0.80901700258254916D0, 0.00000000000000000D0, &
                      0.0000000000000000D0, 0.00000000000000000D0, 1.00000000000000000D0, 0.00000000000000000D0, &
                      0.0000000000000000D0,-0.58778524398803644D0, 0.80901700258254916D0, 0.00000000000000000D0, &
                     0.50000000000000000D0, 0.42532540417601997D0, 0.30901699437494742D0, 0.68819096193209561D0, &
                     0.50000000000000000D0,-0.42532540417601997D0,-0.30901699437494742D0,-0.68819096193209561D0, &
                     0.50000000000000000D0,-0.16245984737382999D0, 0.50000000000000000D0, 0.68819096193209561D0, &
                     0.50000000000000000D0, 0.16245984737382999D0,-0.50000000000000000D0,-0.68819096193209561D0, &
                     0.50000000000000000D0,-0.52573108874869778D0, 0.00000000000000000D0, 0.68819096193209561D0, &
                     0.50000000000000000D0, 0.52573108874869778D0, 0.00000000000000000D0,-0.68819096193209561D0, &
                     0.50000000000000000D0,-0.16245984737382999D0,-0.50000000000000000D0, 0.68819096193209561D0, &
                     0.50000000000000000D0, 0.16245984737382999D0, 0.50000000000000000D0,-0.68819096193209561D0, &
                     0.50000000000000000D0, 0.42532539174817890D0,-0.30901700049082287D0, 0.68819096193209561D0, &
                     0.50000000000000000D0,-0.42532539174817890D0, 0.30901700049082287D0,-0.68819096193209561D0, &
                     0.50000000000000000D0,-0.85065078349635781D0, 0.00000000000000000D0, 0.16245984737382999D0, &
                     0.50000000000000000D0, 0.85065078349635781D0, 0.00000000000000000D0,-0.16245984737382999D0, &
                     0.50000000000000000D0,-0.68819096193209561D0,-0.50000000000000000D0,-0.16245984737382999D0, &
                     0.50000000000000000D0, 0.68819096193209561D0, 0.50000000000000000D0, 0.16245984737382999D0, &
                     0.50000000000000000D0,-0.26286554437434889D0,-0.80901695670145712D0, 0.16245984737382999D0, &
                     0.50000000000000000D0, 0.26286554437434889D0, 0.80901695670145712D0,-0.16245984737382999D0, &
                     0.50000000000000000D0, 0.26286554437434889D0,-0.80901695670145712D0,-0.16245984737382999D0, &
                     0.50000000000000000D0,-0.26286554437434889D0, 0.80901695670145712D0, 0.16245984737382999D0, &
                     0.50000000000000000D0, 0.68819096193209561D0,-0.50000000000000000D0, 0.16245984737382999D0, &
                     0.50000000000000000D0,-0.68819096193209561D0, 0.50000000000000000D0,-0.16245984737382999D0, &
                     0.80901700537708732D0, 0.00000000000000000D0, 0.00000000000000000D0, 0.58778523714932640D0, &
                     0.30901702997862029D0, 0.00000000000000000D0, 0.00000000000000000D0, 0.95105650472681824D0, &
                     0.80901700537708732D0, 0.00000000000000000D0, 0.00000000000000000D0,-0.58778523714932640D0, &
                     0.30901702997862029D0, 0.00000000000000000D0, 0.00000000000000000D0,-0.95105650472681824D0, &
                     0.80901700537708732D0, 0.52573109227969150D0, 0.00000000000000000D0, 0.26286554613984575D0, &
                     0.30901702997862029D0, 0.85065078781948245D0, 0.00000000000000000D0, 0.42532539390974122D0, &
                     0.80901700537708732D0,-0.52573109227969150D0, 0.00000000000000000D0,-0.26286554613984575D0, &
                     0.30901702997862029D0,-0.85065078781948245D0, 0.00000000000000000D0,-0.42532539390974122D0, &
                     0.80901700537708732D0, 0.16245984550474032D0, 0.50000000000000000D0, 0.26286554613984575D0, &
                     0.30901702997862029D0, 0.26286555540853851D0, 0.80901696456355054D0, 0.42532539390974122D0, &
                     0.80901700537708732D0,-0.16245984550474032D0,-0.50000000000000000D0,-0.26286554613984575D0, &
                     0.30901702997862029D0,-0.26286555540853851D0,-0.80901696456355054D0,-0.42532539390974122D0, &
                     0.80901700537708732D0,-0.42532540916195122D0, 0.30901697149092866D0, 0.26286554613984575D0, &
                     0.30901702997862029D0,-0.68819097766197224D0, 0.50000000000000000D0, 0.42532539390974122D0, &
                     0.80901700537708732D0, 0.42532540916195122D0,-0.30901697149092866D0,-0.26286554613984575D0, &
                     0.30901702997862029D0, 0.68819097766197224D0,-0.50000000000000000D0,-0.42532539390974122D0, &
                     0.80901700537708732D0,-0.42532540916195122D0,-0.30901697149092866D0, 0.26286554613984575D0, &
                     0.30901702997862029D0,-0.68819097766197224D0,-0.50000000000000000D0, 0.42532539390974122D0, &
                     0.80901700537708732D0, 0.42532540916195122D0, 0.30901697149092866D0,-0.26286554613984575D0, &
                     0.30901702997862029D0, 0.68819097766197224D0, 0.50000000000000000D0,-0.42532539390974122D0, &
                     0.80901700537708732D0, 0.16245984550474032D0,-0.50000000000000000D0, 0.26286554613984575D0, &
                     0.30901702997862029D0, 0.26286555540853851D0,-0.80901696456355054D0, 0.42532539390974122D0, &
                     0.80901700537708732D0,-0.16245984550474032D0, 0.50000000000000000D0,-0.26286554613984575D0, &
                     0.30901702997862029D0,-0.26286555540853851D0, 0.80901696456355054D0,-0.42532539390974122D0, &

                     ! octagonal QC group 822
                     0.923879532511287D0, 0.D0, 0.D0, 0.38268343236509D0, &       ! 45@[001]
                     0.707106781186547D0, 0.D0, 0.D0, 0.707106781186547D0, &      ! 90@[001]
                     0.38268343236509D0, 0.D0, 0.D0, 0.923879532511287D0, &       ! 135@[001]
                     0.D0, 0.D0, 0.D0, 1.D0, &                                    ! 180@[001]
                    -0.382683432365090D0, 0.D0, 0.D0, 0.923879532511287D0, &     ! 225@[001]
                    -0.707106781186547D0, 0.0D0, 0.0D0, 0.707106781186548D0, &   ! 270@[001]
                    -0.923879532511287D0, 0.0D0, 0.0D0, 0.382683432365090D0, &   ! 315@[001]
                     ! all 2 fold rotation axes
                     0.D0, 1.D0, 0.D0, 0.D0, &                                    ! 180@[100]
                     0.D0, 0.923879532511287D0, 0.38268343236509D0, 0.D0, &       ! 180@[cos(pi/8) sin(pi/8) 0]
                     0.0D0, 0.923879532511287D0, -0.382683432365090D0, 0.0D0, &
                     0.0D0, 0.707106781186548D0, -0.707106781186547D0, 0.0D0, &
                     0.0D0, 0.382683432365090D0, -0.923879532511287D0, 0.0D0, &
                     0.0D0, 0.0D0, -1.0D0, 0.0D0, &
                     0.0D0, -0.382683432365090D0, -0.923879532511287D0, 0.0D0, &
                     0.0D0, -0.707106781186547D0, -0.707106781186548D0, 0.0D0, &

                     ! decagonal QC group 1022
                     0.951056516295154D0, 0.0D0, 0.0D0, 0.309016994374947D0, &    ! 36@[001]
                     0.809016994374947D0, 0.0D0, 0.0D0, 0.587785252292473D0, &    ! 72@[001]
                     0.587785252292473D0, 0.0D0, 0.0D0, 0.809016994374947D0, &    ! 108@[001]
                     0.309016994374947D0, 0.0D0, 0.0D0, 0.951056516295154D0, &    ! 144@[001]
                     0.0D0, 0.0D0, 0.0D0, 1.0D0, &                                ! 180@[001]
                    -0.309016994374947D0, 0.0D0, 0.0D0, 0.951056516295154D0, &    ! 216@[001]
                    -0.587785252292473D0, 0.0D0, 0.0D0, 0.809016994374947D0, &    ! 252@[001]
                    -0.809016994374947D0, 0.0D0, 0.0D0, 0.587785252292473D0, &    ! 288@[001]
                    -0.951056516295154D0, 0.0D0, 0.0D0, 0.309016994374948D0, &    ! 324@[001]
                    ! all 2-fold rotation axis
                     0.D0, 1.D0, 0.D0, 0.D0, &                                    ! 180@[100]
                     0.D0, 0.951056516295154D0, 0.309016994374947D0, 0.D0,   &    ! 180@[cos(pi/10) sin(pi/10) 0]
                     0.0D0, 0.951056516295154D0, -0.309016994374947D0, 0.0D0, &
                     0.0D0, 0.809016994374947D0, -0.587785252292473D0, 0.0D0, &
                     0.0D0, 0.587785252292473D0, -0.809016994374947D0, 0.0D0, &
                     0.0D0, 0.309016994374947D0, -0.951056516295154D0, 0.0D0, &
                     0.0D0, 0.0D0, -1.0D0, 0.0D0, &
                     0.0D0, -0.309016994374947D0, -0.951056516295154D0, 0.0D0, &
                     0.0D0, -0.587785252292473D0, -0.809016994374947D0, 0.0D0, &
                     0.0D0, -0.809016994374947D0, -0.587785252292473D0, 0.0D0, &

                     ! dodecagonal QC group 1222
                     0.965925826289068D0, 0.0D0, 0.0D0, 0.258819045102521D0, &    ! 30@[001]
                     0.866025403784439D0, 0.0D0, 0.0D0, 0.5D0, &                  ! 60@[001]
                     0.707106781186548D0, 0.0D0, 0.0D0, 0.707106781186547D0, &    ! 90@[001]
                     0.5D0, 0.0D0, 0.0D0, 0.866025403784439D0, &                  ! 120@[001]
                     0.258819045102521D0, 0.0D0, 0.0D0, 0.965925826289068D0, &    ! 150@[001]
                     0.0D0, 0.0D0, 0.0D0, 1.0D0, &                                    ! 180@[001]
                    -0.258819045102521D0, 0.0D0, 0.0D0, 0.965925826289068D0, &      ! 210@[001]
                    -0.5D0, 0.0D0, 0.0D0, 0.866025403784439D0, &                  ! 240@[001]
                    -0.707106781186547D0, 0.0D0, 0.0D0, 0.707106781186548D0, &    ! 270@[001]
                    -0.866025403784439D0, 0.0D0, 0.0D0, 0.5D0, &                  ! 300@[001]
                    -0.965925826289068D0, 0.0D0, 0.0D0, 0.258819045102521D0, &    ! 330@[001]
                    ! all 2-fold rotation axes
                    0.0D0, 1.0D0, 0.0D0, 0.0D0, &                                    ! 180@[100]
                    0.0D0, 0.965925826289068D0, 0.258819045102521D0, 0.0D0, &      ! 180@[cos(pi/12) sin(pi/12) 0]
                    0.0D0, 0.965925826289068D0, -0.258819045102521D0, 0.0D0, &
                    0.0D0, 0.866025403784439D0, -0.5D0, 0.0D0, &
                    0.0D0, 0.707106781186548D0, -0.707106781186547D0, 0.0D0, &
                    0.0D0, 0.5D0, -0.866025403784439D0, 0.0D0, &
                    0.0D0, 0.258819045102521D0, -0.965925826289068D0, 0.0D0, &
                    0.0D0, 0.0D0, -1.0D0, 0.0D0, &
                    0.0D0, -0.258819045102521D0, -0.965925826289068D0, 0.0D0, &
                    0.0D0, -0.5D0, -0.866025403784439D0, 0.0D0, &
                    0.0D0, -0.707106781186547D0, -0.707106781186548D0, 0.0D0, &
                    0.0D0, -0.866025403784439D0, -0.5D0, 0.0D0 &
                    /), (/4,152/) )
!DEC$ ATTRIBUTES DLLEXPORT :: SYM_Qsymop

! we overload the conjg, cabs, and .eq. intrinsics
    intrinsic :: conjg, cabs
    public :: conjg, cabs

    interface conjg
      procedure quatconjg
    end interface conjg

    interface cabs
      procedure quatnorm
    end interface cabs

    interface operator(.eq.)
      procedure quatsequal
    end interface

! definition of the quaternion class
  type, public :: Quaternion_T
    !! Quaternion Class definition
    private
      real(kind=sgl), dimension(4) :: q
       !! single precision quaternion
      real(kind=dbl), dimension(4) :: qd
       !! double precision quaternion
      character(1)                 :: s
       !! precision indicator ('s' or 'd')
      real(kind=sgl), dimension(3,3)  :: mu
       !! simplectic transformation array single precision
      real(kind=dbl), dimension(3,3)  :: mud
       !! simplectic transformation array double precision

    contains
    private
! quaternion IO routines
      procedure, pass(self) :: quatprint
      procedure, pass(self) :: getquats
      procedure, pass(self) :: getquatd
      procedure, pass(self) :: setquats
      procedure, pass(self) :: setquatd
      procedure, pass(self) :: setsimplectics
      procedure, pass(self) :: setsimplecticd
! quaternion arithmetic routines
      procedure, pass(self) :: quatflip
      procedure, pass(self) :: quatpos
      procedure, pass(self) :: quatadd
      procedure, pass(self) :: quatsubtract
      procedure, pass(self) :: quatmult
      procedure, pass(self) :: quatsmult
      procedure, pass(self) :: quatsmultd
      procedure, pass(self) :: quatdiv
      procedure, pass(self) :: quatsdiv
      procedure, pass(self) :: quatsdivd
      procedure, pass(self) :: quatconjg
      procedure, pass(self) :: quatnorm
      procedure, pass(self) :: quatnormalize
! quaternion-based transformations
      procedure, pass(self) :: quatLp
      procedure, pass(self) :: quatLpd
      procedure, pass(self) :: quat2simplectic
      procedure, pass(self) :: quat2simplecticd
      procedure, pass(self) :: simplectic2quat
      procedure, pass(self) :: simplectic2quatd
! routines with two or more input quaternions
      procedure, pass(self) :: quatinnerproduct
      procedure, pass(self) :: quatangle
      procedure, pass(self) :: quatslerp
! miscellaneous routines
      procedure, pass(self), public :: quatsequal

      generic, public :: quat_print => quatprint
      generic, public :: quat_flip => quatflip
      generic, public :: quat_pos => quatpos
      generic, public :: get_quats => getquats
      generic, public :: get_quatd => getquatd
      generic, public :: set_quats => setquats
      generic, public :: set_quatd => setquatd
      generic, public :: set_simplectic => setsimplectics, setsimplecticd
      generic, public :: quat_norm => quatnorm
      generic, public :: operator(+) => quatadd
      generic, public :: operator(-) => quatsubtract
      generic, public :: operator(*) => quatmult
      generic, public :: operator(*) => quatsmult, quatsmultd
      generic, public :: operator(/) => quatdiv
      generic, public :: operator(/) => quatsdiv, quatsdivd
      generic, public :: quat_normalize => quatnormalize
      generic, public :: quat_Lp => quatLp, quatLpd
      generic, public :: quat_to_simplectic => quat2simplectic, quat2simplecticd
      generic, public :: simplectic_to_quat => simplectic2quat, simplectic2quatd
      generic, public :: quat_innerproduct => quatinnerproduct
      generic, public :: quat_angle => quatangle
      generic, public :: quat_slerp => quatslerp

  end type Quaternion_T

! the constructor routines for these classes
  interface Quaternion_T
    module procedure Quaternion_constructor
  end interface Quaternion_T

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! We begin with the functions/subroutines that are public in the
! two classes and pair up functions for individual and quaternion
! arrays for easier module maintenance.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
type(Quaternion_T) function Quaternion_constructor( q, qd ) result(Quat)
!DEC$ ATTRIBUTES DLLEXPORT :: Quaternion_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! constructor for the Quaternion Class

IMPLICIT NONE

  real(kind=sgl), INTENT(IN), OPTIONAL      :: q(4)
  real(kind=dbl), INTENT(IN), OPTIONAL      :: qd(4)

! fill in one or the other quaternion
  if ((.not.present(q)).and.(.not.present(qd))) then
    Quat % q = (/ 0.0, 0.0, 0.0, 0.0 /)
    Quat % qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
    Quat % s = 's'
  else
      if (present(q)) then
        Quat % q = q
        Quat % qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        Quat % s = 's'
      end if

      if (present(qd)) then
        Quat % q = (/ 0.0, 0.0, 0.0, 0.0 /)
        Quat % qd = qd
        Quat % s = 'd'
      end if
  end if

end function Quaternion_constructor

!--------------------------------------------------------------------------
subroutine Quaternion_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: Quaternion_destructor
!! author: MDG
!! version: 1.0
!! date: 02/02/20
!!
!! destructor for the Quaternion_T Class

IMPLICIT NONE

type(Quaternion_T), INTENT(INOUT)     :: self

call reportDestructor('Quaternion_T')

end subroutine Quaternion_destructor

!--------------------------------------------------------------------------
recursive subroutine quatprint(self)
!DEC$ ATTRIBUTES DLLEXPORT :: quatprint
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! print a quaternion

use mod_io

IMPLICIT NONE

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion

  type(IO_T)                        :: Message

  if (self%s.eq.'s') then
    call Message % WriteValue('', self%q, 4, frm="('(',4f12.6,'); precision: '$)")
    call Message % WriteValue('',self%s)
  else
    call Message % WriteValue('', self%qd, 4, frm="('(',4f20.14,'); precision: '$)")
    call Message % WriteValue('',self%s)
  end if

end subroutine quatprint

!--------------------------------------------------------------------------
recursive function getquats(self) result(qs)
!DEC$ ATTRIBUTES DLLEXPORT :: getquats
  !! author: MDG
  !! version: 1.0
  !! date: 01/22/20
  !!
  !! return a quaternion

IMPLICIT NONE

class(Quaternion_T),intent(in)    :: self
 !! input quaternion
real(kind=sgl)                    :: qs(4)

qs = self%q

end function getquats

!--------------------------------------------------------------------------
recursive function getquatd(self) result(qd)
!DEC$ ATTRIBUTES DLLEXPORT :: getquatd
  !! author: MDG
  !! version: 1.0
  !! date: 01/22/20
  !!
  !! return a quaternion

IMPLICIT NONE

class(Quaternion_T),intent(in)    :: self
 !! input quaternion
real(kind=dbl)                    :: qd(4)

qd = self%qd

end function getquatd

!--------------------------------------------------------------------------
recursive subroutine setquats(self, qs)
!DEC$ ATTRIBUTES DLLEXPORT :: setquats
  !! author: MDG
  !! version: 1.0
  !! date: 02/18/20
  !!
  !! set a quaternion

IMPLICIT NONE

class(Quaternion_T),intent(inout)    :: self
real(kind=sgl),intent(in)            :: qs(4)
 !! input quaternion

self%q = qs
self%s = 's'

end subroutine setquats

!--------------------------------------------------------------------------
recursive subroutine setquatd(self, qd)
!DEC$ ATTRIBUTES DLLEXPORT :: setquatd
  !! author: MDG
  !! version: 1.0
  !! date: 01/22/20
  !!
  !! set a quaternion

IMPLICIT NONE

class(Quaternion_T),intent(inout)    :: self
real(kind=dbl),intent(in)            :: qd(4)
 !! input quaternion

self%qd = qd
self%s = 'd'

end subroutine setquatd

!--------------------------------------------------------------------------
recursive subroutine setsimplectics(self, mu1, mu2, mu3)
!DEC$ ATTRIBUTES DLLEXPORT :: setsimplectics
  !! author: MDG
  !! version: 1.0
  !! date: 11/22/21
  !!
  !! set the simplectic transformation array

IMPLICIT NONE

class(Quaternion_T),intent(inout)    :: self
real(kind=sgl),intent(in)            :: mu1(3)
real(kind=sgl),intent(in)            :: mu2(3)
real(kind=sgl),intent(in)            :: mu3(3)

self%mu(1,1:3) = mu1
self%mu(2,1:3) = mu2
self%mu(3,1:3) = mu3

end subroutine setsimplectics

!--------------------------------------------------------------------------
recursive subroutine setsimplecticd(self, mu1, mu2, mu3)
!DEC$ ATTRIBUTES DLLEXPORT :: setsimplecticd
  !! author: MDG
  !! version: 1.0
  !! date: 11/22/21
  !!
  !! set the simplectic transformation array

IMPLICIT NONE

class(Quaternion_T),intent(inout)    :: self
real(kind=dbl),intent(in)            :: mu1(3)
real(kind=dbl),intent(in)            :: mu2(3)
real(kind=dbl),intent(in)            :: mu3(3)

self%mud(1,1:3) = mu1
self%mud(2,1:3) = mu2
self%mud(3,1:3) = mu3

end subroutine setsimplecticd

!--------------------------------------------------------------------------
pure recursive subroutine quatflip(self)
!DEC$ ATTRIBUTES DLLEXPORT :: quatflip
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! change the sign of the complete quaternion

IMPLICIT NONE

class(Quaternion_T),intent(inout) :: self

if (self%s.eq.'s') then
  self%q = -self%q
else
  self%qd = -self%qd
end if

end subroutine quatflip

!--------------------------------------------------------------------------
pure recursive subroutine quatpos(self)
!DEC$ ATTRIBUTES DLLEXPORT :: quatpos
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! convert a quaternion to one with a positive scalar part, if it is negative

IMPLICIT NONE

class(Quaternion_T),intent(inout) :: self

if (self%s.eq.'s') then
  if (self%q(1).lt.0.0) self%q = -self%q
else
  if (self%qd(1).lt.0.D0) self%qd = -self%qd
end if

end subroutine quatpos

!--------------------------------------------------------------------------
pure recursive function quatadd(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: quatadd
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! quaternion addition (single/double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in) :: self, y
  type(Quaternion_T)             :: qres

  if (self%s.eq.'s') then
    qres%q = self%q + y%q
    qres%s = 's'
  else
    qres%qd = self%qd + y%qd
    qres%s = 'd'
  end if

end function quatadd

!--------------------------------------------------------------------------
recursive function quatsubtract(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: quatsubtract
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! quaternion subtraction (single/double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in) :: self, y
  type(Quaternion_T)             :: qres

  if (self%s.eq.'s') then
    qres%q = self%q - y%q
    qres%s = 's'
  else
    qres%qd = self%qd - y%qd
    qres%s = 'd'
  end if

end function quatsubtract

!--------------------------------------------------------------------------
pure recursive function quatmult(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: quatmult
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! quaternion multiplication   (single/double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in) :: self, y
   !! input quaternions
  type(Quaternion_T)             :: qres
   !! output quaternion

! the following is a way to reduce the number of multiplications
! needs to be tested and merged with epsijk approach
!
! QuatMul(QUAT *q1, QUAT *q2, QUAT *res){
! float A, B, C, D, E, F, G, H;
! A = (q1->w + q1->x)*(q2->w + q2->x);
! B = (q1->z - q1->y)*(q2->y - q2->z);
! C = (q1->w - q1->x)*(q2->y + q2->z);
! D = (q1->y + q1->z)*(q2->w - q2->x);
! E = (q1->x + q1->z)*(q2->x + q2->y);
! F = (q1->x - q1->z)*(q2->x - q2->y);
! G = (q1->w + q1->y)*(q2->w - q2->z);
! H = (q1->w - q1->y)*(q2->w + q2->z);
! res->w = B + (-E - F + G + H) /2;
! res->x = A - (E + F + G + H)/2;
! res->y = C + (E - F + G - H)/2;
! res->z = D + (E - F - G + H)/2;
! }

  if (self%s.eq.'s') then
    qres%q = (/ self%q(1)*y%q(1) - self%q(2)*y%q(2) -          ( self%q(3)*y%q(3) + self%q(4)*y%q(4) ), &
                self%q(1)*y%q(2) + self%q(2)*y%q(1) + epsijk * ( self%q(3)*y%q(4) - self%q(4)*y%q(3) ), &
                self%q(1)*y%q(3) + self%q(3)*y%q(1) + epsijk * ( self%q(4)*y%q(2) - self%q(2)*y%q(4) ), &
                self%q(1)*y%q(4) + self%q(4)*y%q(1) + epsijk * ( self%q(2)*y%q(3) - self%q(3)*y%q(2) ) /)
    qres%s = 's'
  else
    qres%qd = (/ self%qd(1)*y%qd(1) - self%qd(2)*y%qd(2) -           ( self%qd(3)*y%qd(3) + self%qd(4)*y%qd(4) ), &
                 self%qd(1)*y%qd(2) + self%qd(2)*y%qd(1) + epsijkd * ( self%qd(3)*y%qd(4) - self%qd(4)*y%qd(3) ), &
                 self%qd(1)*y%qd(3) + self%qd(3)*y%qd(1) + epsijkd * ( self%qd(4)*y%qd(2) - self%qd(2)*y%qd(4) ), &
                 self%qd(1)*y%qd(4) + self%qd(4)*y%qd(1) + epsijkd * ( self%qd(2)*y%qd(3) - self%qd(3)*y%qd(2) ) /)
    qres%s = 'd'
  end if

end function quatmult



!--------------------------------------------------------------------------
pure recursive function quatsmult(self, s) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: quatsmult
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! scalar quaternion multiplication   (single precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in)   :: self
   !! input quaternion
  real(kind=sgl), INTENT(IN)       :: s
   !! scalar input
  type(Quaternion_T)               :: qres
   !! output quaternion

  qres%q = (/ s*self%q(1), s*self%q(2), s*self%q(3), s*self%q(4) /)
  qres%s = 's'

end function quatsmult

!--------------------------------------------------------------------------
pure recursive function quatsmultd(self, s) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: quatsmultd
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! scalar quaternion multiplication   (double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in)   :: self
   !! input quaternion
  real(kind=dbl), INTENT(IN)       :: s
   !! scalar input
  type(Quaternion_T)               :: qres
   !! output quaternion

  qres%qd = (/ s*self%qd(1), s*self%qd(2), s*self%qd(3), s*self%qd(4) /)
  qres%s = 'd'

end function quatsmultd

!--------------------------------------------------------------------------
pure recursive function quatconjg(self) result (qres)
!DEC$ ATTRIBUTES DLLEXPORT :: quatconjg
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! quaternion conjugation (extends intrinsic routine conjg)

IMPLICIT NONE

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion
  type(Quaternion_T)                :: qres
   !! output quaternion

  if (self%s.eq.'s') then
    qres%q = (/ self%q(1), -self%q(2), -self%q(3), -self%q(4) /)
    qres%s = 's'
  else
    qres%qd = (/ self%qd(1), -self%qd(2), -self%qd(3), -self%qd(4) /)
    qres%s = 'd'
  end if

end function quatconjg

!--------------------------------------------------------------------------
pure recursive function quatnorm(self) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: quatnorm
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! quaternion norm (extends intrinsic routine abs)

IMPLICIT NONE

  class(Quaternion_T),intent(in) :: self
   !! input quaternion
  real(kind=dbl)                 :: res
   !! output norm

  real(kind=sgl)                 :: n
  real(kind=dbl)                 :: nd, resd

  if (self%s.eq.'s') then
    n = self%q(1)**2 + self%q(2)**2 + self%q(3)**2 + self%q(4)**2
    resd = dsqrt( dble(n) )
    res = dble(sngl(resd))
  else
    nd = self%qd(1)**2 + self%qd(2)**2 + self%qd(3)**2 + self%qd(4)**2
    res = dsqrt( nd )
  end if

end function quatnorm

!--------------------------------------------------------------------------
recursive subroutine quatnormalize(self)
!DEC$ ATTRIBUTES DLLEXPORT :: quatnormalize
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! normalize the input quaternion

IMPLICIT NONE

  class(Quaternion_T),intent(inout) :: self
   !! input quaternion

  type(Quaternion_T)                :: q
  real(kind=sgl)                    :: n
  real(kind=dbl)                    :: nd

  if (self%s.eq.'s') then
    n = self%q(1)**2 + self%q(2)**2 + self%q(3)**2 + self%q(4)**2
    n = sqrt( n )
    q = self%quatsdiv(n)
    self%q = q%q
  else
    nd = self%qd(1)**2 + self%qd(2)**2 + self%qd(3)**2 + self%qd(4)**2
    nd = sqrt( nd )
    q = self%quatsdivd(nd)
    self%qd = q%qd
  end if

end subroutine quatnormalize

!--------------------------------------------------------------------------
recursive function quatdiv(self, y) result (qres)
!DEC$ ATTRIBUTES DLLEXPORT :: quatdiv
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! quaternion division (single/double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in)    :: self, y
   !! input quaternions
  type(Quaternion_T)                :: qres
   !! output quaternion

  type(Quaternion_T)                :: p, cy
  real(kind=sgl)                    :: q
  real(kind=dbl)                    :: qd

  if (self%s.eq.'s') then
      q = quatnorm(y)
      cy = quatconjg(y)
      p = quatsdiv( cy, q*q )
      qres = quatmult(self,p)
      qres%s = 's'
  else
      qd = quatnorm(y)
      cy = quatconjg(y)
      p = quatsdivd( cy, qd*qd )
      qres = quatmult(self,p)
      qres%s = 'd'
  end if


end function quatdiv

!--------------------------------------------------------------------------
recursive function quatsdiv(self, s) result (qres)
!DEC$ ATTRIBUTES DLLEXPORT :: quatsdiv
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! quaternion division (single precision)

use mod_io

IMPLICIT NONE

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=sgl), INTENT(IN)        :: s
   !! input quaternion (denominator)

  type(Quaternion_T)                :: qres
  type(IO_T)                        :: Message

  if (s.ne.0.0) then
      qres%q = (/ self%q(1)/s, self%q(2)/s, self%q(3)/s, self%q(4)/s /)
      qres%s = 's'
  else
    call Message % printWarning('quatsdiv', (/ 'Attempting to divide quaternion by zero; skipping operation ...' /) )
    qres = self
  end if

end function quatsdiv

!--------------------------------------------------------------------------
recursive function quatsdivd(self, s) result (qres)
!DEC$ ATTRIBUTES DLLEXPORT :: quatsdivd
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! quaternion division (doubgle precision)

use mod_io

IMPLICIT NONE

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=dbl), INTENT(IN)        :: s
   !! input quaternion (denominator)

  type(Quaternion_T)                :: qres
  type(IO_T)                        :: Message

  if (s.ne.0.0) then
    qres%qd = (/ self%qd(1)/s, self%qd(2)/s, self%qd(3)/s, self%qd(4)/s /)
    qres%s = 'd'
  else
    call Message % printWarning('quatsdivd', (/ 'Attempting to divide quaternion by zero; skipping operation ...' /) )
    qres = self
  end if

end function quatsdivd

!--------------------------------------------------------------------------
pure recursive function quatinnerproduct(self, y) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: quatinnerproduct
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! quaternion inner product (single precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in)    :: self, y
   !! input quaternions
  real(kind=dbl)                    :: res
   !! inner product

  if (self%s.eq.'s') then
    res = dble(self%q(1) * y%q(1) + self%q(2) * y%q(2) + self%q(3) * y%q(3) + self%q(4) * y%q(4))
  else
    res = self%qd(1) * y%qd(1) + self%qd(2) * y%qd(2) + self%qd(3) * y%qd(3) + self%qd(4) * y%qd(4)
  end if

end function quatinnerproduct

!--------------------------------------------------------------------------!
pure recursive function quatangle(self, y) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: quatangle
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! interquaternion angle   (single/double precision)
  !! this only has meaning for a unit quaternion, so we test first and return -10000.0
  !! if either of the quaternions is not a unit quaternion.

IMPLICIT NONE

  class(Quaternion_T),intent(in)    :: self, y
   !! input quaternions
  real(kind=dbl)                    :: res
   !! angle (radians)

  real(kind=sgl)                    :: q, nself, ny
  real(kind=dbl)                    :: qd, nselfd, nyd

  if (self%s.eq.'s') then
      nself = self%quatnorm()
      ny = y%quatnorm()
      q = self%quat_innerproduct(y)
      res = dble(acos( q/(nself * ny) ))
  else
      nselfd = self%quatnorm()
      nyd = y%quatnorm()
      qd = self%quat_innerproduct(y)
      res = dacos( qd/(nselfd * nyd) )
  end if

end function quatangle

!--------------------------------------------------------------------------!
recursive function quatLp(self, v) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: quatLp
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! actively rotate a unit vector by a unit quaternion, L_p = p v p* (single precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion
  real(kind=sgl),intent(in)         :: v(3)
   !! input vector to be rotated
  real(kind=sgl)                    :: res(3)
   !! output vector

  type(Quaternion_T)                :: qv, rqv, cq

  qv%q = (/ 0.0, v(1), v(2), v(3) /)
  qv%s = 's'
  cq = quatconjg(self)
  rqv = quatmult(self, quatmult(qv, cq) )
  res(1:3) = rqv%q(2:4)

end function quatLp

!--------------------------------------------------------------------------!
! pure recursive function quatLpd(self, v) result (res)
recursive function quatLpd(self, v) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: quatLpd
  !! author: MDG
  !! version: 1.0
  !! date: 01/03/20
  !!
  !! actively rotate a unit vector by a unit quaternion, L_p = p v p* (double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion
  real(kind=dbl),intent(in)         :: v(3)
   !! input vector to be rotated
  real(kind=dbl)                    :: res(3)
   !! output vector

  type(Quaternion_T)                :: qv, rqv, cq

  qv%qd = (/ 0.D0, v(1), v(2), v(3) /)
  qv%s = 'd'
  cq = quatconjg(self)
  rqv = quatmult(self, quatmult(qv, cq) )
  res(1:3) = rqv%qd(2:4)

end function quatLpd

!--------------------------------------------------------------------------!
recursive function quatLpd_vecarray(self, N, v) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: quatLpd_vecarray
  !! author: MDG
  !! version: 1.0
  !! date: 06/02/21
  !!
  !! actively rotate an array of unit vectors by a unit quaternion, L_p = p v p* (double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion
  integer(kind=irg),intent(in)      :: N
  real(kind=dbl),intent(in)         :: v(3, N)
   !! input vector to be rotated
  real(kind=dbl)                    :: res(3, N)
   !! output vector

  type(Quaternion_T)                :: qv, rqv, cq
  integer(kind=irg)                 :: i

  do i=1,N 
    qv%qd = (/ 0.D0, v(1, i), v(2, i), v(3, i) /)
    qv%s = 'd'
    cq = quatconjg(self)
    rqv = quatmult(self, quatmult(qv, cq) )
    res(1:3, i) = rqv%qd(2:4)
  end do

end function quatLpd_vecarray

!--------------------------------------------------------------------------!
recursive subroutine quat2simplectic(self, c1, c2)
!DEC$ ATTRIBUTES DLLEXPORT :: quat2simplectic
  !! author: MDG
  !! version: 1.0
  !! date: 11/22/21
  !!
  !! split a quaternion into two complex numbers (single precision)
  !! Note: the simplectic transformation array has to be set first !

IMPLICIT NONE

  class(Quaternion_T),intent(in)    :: self
  complex(kind=sgl),INTENT(INOUT)   :: c1
  complex(kind=sgl),INTENT(INOUT)   :: c2

  real(kind=sgl)                    :: v(3)

  v = self%q(2:4)

! simplex part
  c1 = cmplx(self%q(1), DOT_PRODUCT(v,self%mu(1,1:3)))

! perplex part
  c2 = cmplx(DOT_PRODUCT(v,self%mu(2,1:3)), DOT_PRODUCT(v,self%mu(3,1:3)))

end subroutine quat2simplectic

!--------------------------------------------------------------------------!
recursive subroutine quat2simplecticd(self, c1, c2)
!DEC$ ATTRIBUTES DLLEXPORT :: quat2simplecticd
  !! author: MDG
  !! version: 1.0
  !! date: 11/22/21
  !!
  !! split a quaternion into two complex numbers (double precision)
  !! Note: the simplectic transformation array has to be set first !

IMPLICIT NONE

  class(Quaternion_T),intent(in)    :: self
  complex(kind=dbl),INTENT(INOUT)   :: c1
  complex(kind=dbl),INTENT(INOUT)   :: c2

  real(kind=dbl)                    :: v(3)

  v = self%qd(2:4)

! simplex part
  c1 = cmplx(self%qd(1), DOT_PRODUCT(v,self%mud(1,1:3)))

! perplex part
  c2 = cmplx(DOT_PRODUCT(v,self%mud(2,1:3)), DOT_PRODUCT(v,self%mud(3,1:3)))

end subroutine quat2simplecticd

!--------------------------------------------------------------------------!
recursive subroutine simplectic2quat(self, c1, c2)
!DEC$ ATTRIBUTES DLLEXPORT :: simplectic2quat
  !! author: MDG
  !! version: 1.0
  !! date: 11/22/21
  !!
  !! perform an inverse simplectic transformation (single precision)
  !! Note: the simplectic transformation array has to be set first !

use, intrinsic :: iso_c_binding 

IMPLICIT NONE

  class(Quaternion_T),INTENT(INOUT) :: self
  complex(kind=sgl),INTENT(INOUT)   :: c1
  complex(kind=sgl),INTENT(INOUT)   :: c2

  real(kind=sgl)                    :: xyz(3)

  xyz = (/ aimag(c1), real(c2), aimag(c2) /)

  self%q = (/ real(c1), DOT_PRODUCT(xyz, self%mu(1:3,1)), DOT_PRODUCT(xyz, self%mu(1:3,2)), DOT_PRODUCT(xyz, self%mu(1:3,3)) /)

end subroutine simplectic2quat

!--------------------------------------------------------------------------!
recursive subroutine simplectic2quatd(self, c1, c2)
!DEC$ ATTRIBUTES DLLEXPORT :: simplectic2quatd
  !! author: MDG
  !! version: 1.0
  !! date: 11/22/21
  !!
  !! perform an inverse simplectic transformation (double precision)
  !! Note: the simplectic transformation array has to be set first !

IMPLICIT NONE

  class(Quaternion_T),INTENT(INOUT) :: self
  complex(kind=dbl),INTENT(IN)      :: c1
  complex(kind=dbl),INTENT(IN)      :: c2

  real(kind=dbl)                    :: xyz(3)

  xyz = (/ aimag(c1), real(c2), aimag(c2) /)

  self%qd = (/ dble(c1), DOT_PRODUCT(xyz, self%mu(1:3,1)), DOT_PRODUCT(xyz, self%mu(1:3,2)), DOT_PRODUCT(xyz, self%mu(1:3,3)) /)

end subroutine simplectic2quatd



!--------------------------------------------------------------------------!
! pure recursive function quatslerp(self, qb, n) result(res)
recursive function quatslerp(self, qb, n) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: quatslerp
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! return an array of interpolated quaternions

IMPLICIT NONE

  class(Quaternion_T),intent(in)         :: self   ! = qa
   !! input quaternion (start)
  class(Quaternion_T),intent(in)         :: qb
   !! input quaternion (end)
  integer(kind=irg),intent(in)           :: n
   !! number of steps in the interpolation
  type(Quaternion_T)                     :: res(n)
   !! output interpolated quaternion list

  type(Quaternion_T)                     :: cqa
  real(kind=sgl)                         :: theta, phi, dphi, s
  real(kind=dbl)                         :: thetad, phid, dphid, sd
  integer(kind=irg)                      :: i

  if (self%s.eq.'s') then
      do i=1,n
        res(i)%q = (/ 0.0, 0.0, 0.0, 0.0 /)
        res(i)%s = 's'
      end do
      cqa = quatconjg(self)
      theta = acos( quatinnerproduct(qb, cqa ))
      if (theta.ne.0.0) then
        s = 1.0/sin(theta)
        dphi = theta/real(n-1)

        do i=1,n
          phi = real(i-1)*dphi
          res(i) = quatsmult(self, sin(theta-phi)*s ) + quatsmult(qb, sin(phi)*s )
        end do
      else
        do i=1,n
          res(i) = self
        end do
      end if
  else
      do i=1,n
        res(i)%qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        res(i)%s = 'd'
      end do
      cqa = quatconjg(self)
      thetad = acos( quatinnerproduct(qb, cqa ))
      if (thetad.ne.0.D0) then
        sd = 1.D0/sin(theta)
        dphi = theta/dble(n-1)

        do i=1,n
          phi = dble(i-1)*dphi
          res(i) = quatsmultd(self, sin(theta-phi)*sd ) + quatsmultd(qb, sin(phi)*sd )
        end do
      else
        do i=1,n
          res(i) = self
        end do
      end if
  end if

end function quatslerp

!--------------------------------------------------------------------------
recursive function quatsequal(self, qb) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: quatsequal
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! quaternion comparison (double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in)    :: self, qb
   !! input quaternions
  logical                           :: res

  type(Quaternion_T)                :: diff

  real(kind=sgl)                    :: d, eps=1.0e-7
  real(kind=dbl)                    :: dd, epsd=1.0e-12

  res = .TRUE.
  diff = self - qb

  if (self%s.eq.'s') then
    d = maxval( abs( diff%q(:) ) )
    if (d.gt.eps) res = .FALSE.
  else
    dd = maxval( abs( diff%qd(:) ) )
    if (dd.gt.epsd) res = .FALSE.
  end if

end function quatsequal

end module mod_quaternions
