!------------------------------------------------------------------------------
!P+
! NAME:
!       s3regrid
!
! PURPOSE:
!>      Process a SLSTR SAFE directory and regrid the Vis/NIR channels (S1-6)
!>      onto the i-grid
!
! CATEGORY:
!       SLSTR pre-processing
!
! LANGUAGE:
!       Fortran-2003
!
! MODULES:
!       SLSTR_Preprocessor
!
! NOTES:
!
! CREATION HISTORY:
!       Written by:   Owen Embury 12/11/2021
!                     University of Reading
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
!P-
!------------------------------------------------------------------------------
PROGRAM s3regrid
  ! ------------
  ! Modules used
  ! ------------
  IMPLICIT NONE

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(len=256) :: s3path, arg
  LOGICAL :: skiparg
  LOGICAL :: usebstripe
  INTEGER :: iarg, ilen, ios
  INTEGER :: npixel

  ! Defaults
  s3path = ''
  npixel = 5
  usebstripe = .FALSE.

  skiparg = .FALSE.

  DO iarg = 1, COMMAND_ARGUMENT_COUNT()
    IF (skiparg) THEN
      skiparg = .FALSE.
      CYCLE
    END IF

    CALL GET_COMMAND_ARGUMENT(iarg, arg, ilen)

    SELECT CASE(arg(1:ilen))
    ! -- Set number of pixels to use in regridding
    CASE('-n')
      CALL GET_COMMAND_ARGUMENT(iarg+1, arg, ilen)
      READ(arg, fmt=*, iostat=ios) npixel
      IF (ios/=0) THEN
        WRITE(*,*) "s3regrid -n expects integer, got: ", arg(1:ilen)
        STOP
      END IF
      skiparg = .True.

    ! -- Enable use of b-stripe for S4-S6
    CASE('-b', '--use-bstripe')
      usebstripe = .True.

    CASE DEFAULT
      IF (LEN_TRIM(s3path) > 0) THEN
        WRITE(*,*) "Unexpected command line argument: ", arg(1:ilen)
        STOP
      END IF
      s3path = arg

    END SELECT
  END DO

  IF (LEN_TRIM(s3path) == 0) THEN
    WRITE(*,*) "usage: s3regrid [-n N] [-b,--use-bstripe] path"
    STOP
  END IF

  CALL process_view(s3path, 'n', npixel, usebstripe)
  CALL process_view(s3path, 'o', npixel, usebstripe)

CONTAINS
  SUBROUTINE process_view(path, view, nnear, use_bstripe)
    USE netcdf
    USE GbcsKinds
    USE GbcsNetCDF
    USE GbcsPath
    USE SLSTR_Preprocessor
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    CHARACTER(LEN=*), INTENT(in) :: path
    CHARACTER(LEN=1), INTENT(in) :: view
    INTEGER,          INTENT(in) :: nnear
    LOGICAL,          INTENT(in) :: use_bstripe

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: status, ncid, dim1, dim2, varid
    INTEGER :: iband
    INTEGER, DIMENSION(2) :: irsize
    CHARACTER(len=14) :: name
    CHARACTER(len=64) :: attr_long, attr_name, attr_unit
    INTEGER(KIND=GbcsInt2) :: fill
    TYPE(NEIGHBOURHOOD_MAP) :: nnmap, nnmap2
    REAL :: scale, offset
    REAL, ALLOCATABLE, DIMENSION(:,:) :: rads

    WRITE(*,*) "Calculating nearest neighbours: ", view
    CALL compute_scene_neighbourhood(path, view, 'a', nnmap)
    IF (use_bstripe) THEN
      WRITE(*,*) "Adding b-stripe"
      CALL compute_scene_neighbourhood(path, view, 'b', nnmap2)
      CALL merge_neighbourhoods(nnmap, nnmap2)
    END IF
    irsize = SHAPE(nnmap%entries)
    WRITE(*,*) "Output size: ", irsize
    ALLOCATE(rads(irsize(1), irsize(2)))

    DO iband = 1,6
      ! Read attributes from input radiance
      WRITE(name, '("S", I1, "_radiance_a", A1)') iband, view
      status = nf90_open(Path_Join(path, name//'.nc'), NF90_NOWRITE, ncid)
      status = nf90_inq_varid(ncid, name, varid)
      status = nf90_get_att(ncid, varid, '_FillValue', fill)
      status = nf90_get_att(ncid, varid, 'add_offset', offset)
      status = nf90_get_att(ncid, varid, 'scale_factor', scale)
      status = nf90_get_att(ncid, varid, 'long_name', attr_long)
      status = nf90_get_att(ncid, varid, 'standard_name', attr_name)
      status = nf90_get_att(ncid, varid, 'units', attr_unit)
      status = nf90_close(ncid)

      WRITE(name, '("S", I1, "_radiance_i", A1)') iband, view
      WRITE(*,*) 'Channel: ', name
      IF (use_bstripe .AND. iband > 3) THEN
        CALL process_scene_band(view, path, iband, rads, nnmap2, nnear*2, FUNCTION_MEAN)
      ELSE
        CALL process_scene_band(view, path, iband, rads, nnmap, nnear, FUNCTION_MEAN)
      END IF
      status = nf90_create(Path_Join(path, name//'.nc'), NF90_NETCDF4, ncid)
      status = nf90_def_dim(ncid, "columns", irsize(1), dim1)
      status = nf90_def_dim(ncid, "rows", irsize(2), dim2)
      status = nf90_def_var(ncid, name, nf90_short, [dim1, dim2], varid, &
                            chunksizes=irsize, deflate_level=2, shuffle=.True.)
      status = nf90_put_att(ncid, varid, '_FillValue', fill)
      status = nf90_put_att(ncid, varid, 'add_offset', offset)
      status = nf90_put_att(ncid, varid, 'scale_factor', scale)
      status = nf90_put_att(ncid, varid, 'valid_min', -3000_GbcsInt2)
      status = nf90_put_att(ncid, varid, 'valid_max', 32767_GbcsInt2)
      status = nf90_put_att(ncid, varid, 'long_name', attr_long)
      status = nf90_put_att(ncid, varid, 'standard_name', attr_name)
      status = nf90_put_att(ncid, varid, 'units', attr_unit)
      status = nf90_enddef(ncid)

      status = Write_Var(ncid, varid, rads, [1, 1])
      status = nf90_close(ncid)
    END DO

  END SUBROUTINE process_view

END PROGRAM s3regrid
