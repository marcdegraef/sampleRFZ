! ###################################################################
! Copyright (c) 2014-2022, Marc De Graef Research Group/Carnegie Mellon University
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

program EMsampleRFZ
  !! author: MDG
  !! version: 1.0 
  !! date: 02/04/22
  !!
  !! Basic program to generate a uniform sampling of Rodrigues Fundamental Zone
  !! Extracted from the Object Oriented EMsoftOO version 6.0 package
  !!
  !! This program calls the SampleRFZ routine of the so3 module to generate
  !! an angle file of orientations for points that uniformly sample an RFZ for a given
  !! crystal symmetry. 

use mod_kinds
use mod_global
use mod_sampleRFZ

IMPLICIT NONE

type(sampleRFZ_T)   :: RFZ 

integer(kind=irg)   :: numarg       ! number of command line arguments
integer(kind=irg)   :: iargc        ! external function for command line
character(fnlen)    :: arg          ! to be read from the command line
character(fnlen)    :: nmlfile      ! nml file name

! get the input file name which is an argument on the command line
numarg = command_argument_count()
if (numarg.ne.1) then 
  write (*,*) 'This program requires a single input file based on sampleRFZ.template'
  stop 
end if 

call get_command_argument(numarg,nmlfile)

 ! deal with the namelist stuff
RFZ = sampleRFZ_T(nmlfile)

! perform the sampling algorithm
call RFZ%CreateSampling()

end program EMsampleRFZ
