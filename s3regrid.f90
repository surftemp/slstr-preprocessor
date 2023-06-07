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
  CHARACTER(len=256) :: s3path

  CALL GET_COMMAND_ARGUMENT(1, s3path)

  CALL process_view(s3path, 'n', 5)
  CALL process_view(s3path, 'o', 5)

CONTAINS
  SUBROUTINE process_view(path, view, nnear)
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

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: status, ncid, dim1, dim2, varid
    INTEGER :: iband
    INTEGER, DIMENSION(2) :: irsize
    CHARACTER(len=14) :: name
    CHARACTER(len=64) :: attr_long, attr_name, attr_unit
    INTEGER(KIND=GbcsInt2) :: fill
    TYPE(NEIGHBOURHOOD_MAP) :: nnmap
    REAL :: scale, offset
    REAL, ALLOCATABLE, DIMENSION(:,:) :: rads

    WRITE(*,*) "Calculating nearest neighbours: ", view
    CALL compute_scene_neighbourhood(view, path, nnmap, 'a')
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
      CALL process_scene_band(view, path, iband, rads, nnmap, nnear, FUNCTION_MEAN)
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
