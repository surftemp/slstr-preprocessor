!------------------------------------------------------------------------
!
!    Copyright 2021 Owen Embury, and Claire Bulgin
!    Department of Meteorology, University of Reading, UK
!
!    This file is part of the GBCS software package.
!
!    The GBCS software package is free software: you can redistribute it
!    and/or modify it under the terms of the GNU General Public License
!    as published by the Free Software Foundation, either version 3 of
!    the License, or (at your option) any later version.
!
!    The GBCS software package is distributed in the hope that it will be
!    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the GBCS software package.
!    If not, see <http://www.gnu.org/licenses/>.
!
!------------------------------------------------------------------------



!+ This module contains the definition of GBCS real data kinds

MODULE GbcsKinds

!
! Description:
!
!  Define byte size limits for generic kinds used in the GBCS code
!
! Method:
!
!  Use the FORTRAN SELECTED KIND functions to define the following kinds
!
!    1 byte integer
!    2 byte integer
!    4 byte integer
!    8 byte integer
!
!    4 byte real
!    8 byte real
!
!
! Owner: Manager of TRICS Project
!
! History:
! Version  Date       Comment
! -------  ----       -------
! 0.0   03/01/2005    Creation                                                   CPO
!
!
! Code Description:
!   Language:            Fortran 90
!
!
! Institute of Atmospheric and Environmental Sciences
! The University of Edinburgh, The Crew Building, The King's Buildings
! West Mains Road, Edinburgh, UK  EH9 3JN
!

! Declarations:

! Modules used:

  IMPLICIT NONE

  CHARACTER(LEN=50), PARAMETER, PRIVATE :: Module_Name = 'GbcsKinds.mod'

! Global Declarations:

  INTEGER, PARAMETER :: GbcsInt1 = SELECTED_INT_KIND( 2 )                     ! 1 byte integer
  INTEGER, PARAMETER :: GbcsInt2 = SELECTED_INT_KIND( 4 )                     ! 2 byte integer
  INTEGER, PARAMETER :: GbcsInt4 = SELECTED_INT_KIND( 9 )                     ! 4 byte integer
  INTEGER, PARAMETER :: GbcsInt8 = SELECTED_INT_KIND( 18 )                    ! 8 byte integer

  INTEGER, PARAMETER :: GbcsReal = SELECTED_REAL_KIND( P =  6 , R =  37 )     ! 4 byte real
  INTEGER, PARAMETER :: GbcsDble = SELECTED_REAL_KIND( P = 13 , R = 300 )     ! 8 byte real

! Global named constants:

! Global Type Definitions:

! Various lengths and limits

END MODULE GbcsKinds
