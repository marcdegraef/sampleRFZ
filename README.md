# sampleRFZ (extracted from version 6.0 of EMsoftOO)

This code base contains only one single program called "sampleRFZ" which generates a uniform sampling of SO(3), limited to the Fundamental Zone, for arbitrary crystal symmetry. The code is extracted from the much larger EMsoftOO package, and Object Oriented fortran 2018 package for electron microscopy etc.

This program can be run from the command line as follows:
```fortran
sampleRFZ sampleRFZ.nml
```
where the argument is the name of a namelist file with the following structure:
```fortran
 &RFZlist
! template file for the sampleRFZ program
!
! point group number
! #1:   '1', #2:   '-1', #3:     '2', #4:   'm', #5: '2/m', #6:  '222', #7:   'mm2', #8:  'mmm'
! #9:   '4', #10:  '-4', #11:  '4/m', #12:'422', #13:'4mm', #14:'-42m', #15:'4/mmm', #16:   '3'
! #17: '-3', #18:  '32', #19:   '3m', #20:'-3m', #21:  '6', #22:  '-6', #23:  '6/m', #24: '622'
! #25:'6mm', #26:'-6m2', #27:'6/mmm', #28: '23', #29: 'm3', #30: '432', #31: '-43m', #32:'m-3m'
 pgnum = 32,
! number of sampling points along cube semi-edge
 nsteps = 100,
! grid type: 0 contains origin, 1 is has origin at center of grid box
 gridtype = 0,
! euler angle output file name
 euoutname = 'undefined',
! cubochoric output file name
 cuoutname = 'undefined',
! homochoric output file name
 hooutname = 'undefined',
! Rodrigues output file name
 rooutname = 'undefined',
! quaternion output file name
 quoutname = 'undefined',
! orientation matrix output file name
 omoutname = 'undefined',
! axis angle pair output file name
 axoutname = 'undefined',
! stereographic output file name
 stoutname = 'undefined',
! rotation vector output file name
 rvoutname = 'undefined',
 /
```
A template for this file can be found in the NameListTemplates folder.  The user should edit the point group number *pgnum*, the number of sampling steps along the semi-edge of the cubochoric cube *nsteps*, the gridtype *gridtype* and one or more output filenames, depending on how many different orientations representations are needed.  For *nsteps=100*, SO(3) will be uniformly sampled with a step size of about 1.4°; more details can be found in the following two papers: [D. Rosca, A. Morawiec, and M. De Graef. "A new method of constructing a grid in the space of 3D rotations and its applications to texture analysis". Modeling and Simulations in Materials Science and Engineering 22, 075013 (2014)](https://doi.org/10.1088/0965-0393/22/7/075013) and [S. Singh and M. De Graef, "Orientation sampling for dictionary-based diffraction pattern indexing methods". MSMSE 24, 085013 (2016)](https://doi.org/10.1088/0965-0393/24/8/085013)

## Compilation Instructions ##
This is a standalone package that requires a recent fortran compiler (ifort or gfortran), GIT, BLAS-LAPACK (installed in the usual location for your platform) as well as a recent version of the cmake package (e.g., version 3.14.6).  The cmake files should allow for this package to be compiled on Mac OSX, Windows, or any flavor of Linux (but let me know if you run into issues).  It has been tested on Mac OS X 10.14.6.  

After cloning the repository onto your platform into a folder called (say) *sampleRFZ*, create a second folder called *sampleRFZbuild* at the same level in your folder hierarchy. Go in to the *sampleRFZbuild* folder and execute the following cmake command:
```fortran
cmake -DCMAKE_BUILD_TYPE=Release ../sampleRFZ
```
This will generate a lot of output but there should not be any error messages.  Then type
```fortran
make -j
```
(or nmake on Windows) and the *sampleRFZ* program should be compiled along with the *sampleRFZLib* library.



## Licenses ##

	!###################################################################
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

