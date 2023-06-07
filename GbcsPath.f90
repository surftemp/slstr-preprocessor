!------------------------------------------------------------------------------
!M+
! NAME:
!       GbcsPath
!
! PURPOSE:
!>      Module contining some useful functions on pathnames
!
! CATEGORY:
!       Utility
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       USE GbcsPath
!
! PUBLIC DATA:
!       None
!
! MODULES:
!       None
!
! CONTAINS:
!       Path_Join
!
! DERIVED TYPES:
!       None
!
! NOTES:
!
! CREATION HISTORY:
!       Written by:   Owen Embury 11/10/2011
!                     IAES, University of Edinburgh
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
!M-
!------------------------------------------------------------------------------
MODULE GbcsPath
  ! ------------
  ! Modules used
  ! ------------

  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Path_Join

  !-------------------
  ! Module variables
  !-------------------
  SAVE
  CHARACTER(LEN=:), ALLOCATABLE :: gbcs_path, gbcs_home

  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER( * ), PRIVATE, PARAMETER :: MODULE_NAME = &
    'GbcsMod_Misc/GbcsPath.f90'

CONTAINS

!------------------------------------------------------------------------------
!F+
! NAME:
!       Path_Join
!
! PURPOSE:
!>       Join two path segments together.
!>
!>       e.g. path + '/' + filename.
!>       If the first path is empty, just returns second part.
!>       If second part is absolute (starts with '/') then ignores first.
!
! CATEGORY:
!       File Utility
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       fullpath = Path_Join(path1, path2)
!
! ARGUMENTS:
!>@ARG{path1, in, CHARACTER(*)} Path component
!>@ARG{path2, in, CHARACTER(*)} Path component
!
! FUNCTION RESULT:
!>@result two input paths joined using the path separator '/'
!
! CALLS:
!       Path_Join_Calc_Length
!
! SIDE EFFECTS:
!       None
!
! EXAMPLE:
!> @EX
!>      CHARACTER(len=50) :: path
!>      path = Path_Join('/home/user', 'myfile')
!> @endcode
!
! CREATION HISTORY:
!       Written by:   Owen Embury 11/10/2011
!                     IAES, University of Edinburgh
!F-
!------------------------------------------------------------------------------
   PURE FUNCTION Path_Join(path1, path2) RESULT(path)
    IMPLICIT NONE

    ! ---------
    ! Arguments
    ! ---------
    CHARACTER(LEN=*), INTENT(IN) :: path1
    CHARACTER(LEN=*), INTENT(IN) :: path2

    ! ---------
    ! Function result
    ! ---------
    CHARACTER(LEN=Path_Join_Calc_Length(path1,path2)) :: path

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: length

    length = LEN_TRIM(path1)

    IF (length==0 .OR. INDEX(path2,'/')==1) THEN
      path = TRIM(path2)
      RETURN
    END IF

    IF (length == INDEX(path1,'/',.TRUE.)) THEN
      path = TRIM(path1) // TRIM(path2)
    ELSE
      path = TRIM(path1) // "/" // TRIM(path2)
    END IF

  END FUNCTION Path_Join

!> Internal function used by path_join()
  ELEMENTAL FUNCTION Path_Join_Calc_Length(path1, path2) RESULT(length)
    IMPLICIT NONE

    ! ---------
    ! Arguments
    ! ---------
    CHARACTER(LEN=*), INTENT(IN) :: path1
    CHARACTER(LEN=*), INTENT(IN) :: path2

    ! ---------
    ! Function result
    ! ---------
    INTEGER :: length

    length = LEN_TRIM(path1)

    IF (length==0 .OR. INDEX(path2,'/')==1) THEN
      length = LEN_TRIM(path2)
      RETURN
    END IF

    IF (length /= INDEX(path1,'/',.TRUE.)) THEN
      length = length + 1
    END IF

    length = length + LEN_TRIM(path2)

  END FUNCTION Path_Join_Calc_Length

END MODULE GbcsPath
