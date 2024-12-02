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
  INTEGER :: i1, i2, nchan
  INTEGER :: npixel
  INTEGER, DIMENSION(6) :: chans

  ! Defaults
  s3path = ''
  npixel = 5
  usebstripe = .FALSE.
  chans = -1

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
        STOP 1
      END IF
      skiparg = .True.

    ! -- Enable use of b-stripe for S4-S6
    CASE('-b', '--use-bstripe')
      usebstripe = .True.

    CASE('-s', '--standard-deviation')
      skiparg = .True.
      CALL GET_COMMAND_ARGUMENT(iarg+1, arg, ilen)
      IF (arg == 'all') THEN
        chans = [1, 2, 3, 4, 5, 6]
        CYCLE
      END IF
      ! Count number of channels listed
      nchan = 1
      i1 = 1
      DO
        i2 = INDEX(arg(i1:), ',')
        IF (i2 == 0) EXIT
        i1 = i1 + i2
        nchan = nchan + 1
      END DO
      READ(arg, *, iostat=ios) chans(1:nchan)
      IF (nchan > 6 .OR. ios /= 0) THEN
        WRITE(*,*) "Error parsing channel list: ", arg(1:ilen)
        STOP 1
      END IF
      IF (ANY(chans(1:nchan) < 1 .OR. chans(1:nchan) > 6)) THEN
        WRITE(*,*) "Invalid channel numbers: ", chans(1:nchan)
        STOP 1
      END IF

    CASE DEFAULT
      IF (arg(1:1) == '-') THEN
        WRITE(*,*) "invalid option: ", TRIM(arg)
        STOP 1
      END IF
      IF (LEN_TRIM(s3path) > 0) THEN
        WRITE(*,*) "Unexpected command line argument: ", arg(1:ilen)
        STOP 1
      END IF
      s3path = arg

    END SELECT
  END DO

  IF (LEN_TRIM(s3path) == 0) THEN
    WRITE(*,*) "usage: s3regrid [-n N] [-b,--use-bstripe] [-s,--standard-deviation CHANNELS] path"
    STOP 1
  END IF

  CALL process_view(s3path, 'n', npixel, usebstripe, chans)
  CALL process_view(s3path, 'o', npixel, usebstripe, chans)

CONTAINS
  SUBROUTINE process_view(path, view, nnear, use_bstripe, channels)
    USE netcdf
    USE GbcsKinds
    USE GbcsPath
    USE SLSTR_Preprocessor
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    CHARACTER(LEN=*), INTENT(in) :: path        !< Directory containing SLSTR files
    CHARACTER(LEN=1), INTENT(in) :: view        !< View to process ('n' or 'o')
    INTEGER,          INTENT(in) :: nnear       !< Number of pixels to use in average
    LOGICAL,          INTENT(in) :: use_bstripe !< Include b-stripe
    INTEGER, DIMENSION(:), INTENT(in) :: channels !< Channels to calculate subpixel texture

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: status, ncid, dim1, dim2, v_rad, v_std
    INTEGER :: iband
    INTEGER, DIMENSION(2) :: irsize
    CHARACTER(len=14) :: name
    CHARACTER(len=80) :: attr_long, attr_name, attr_unit
    INTEGER(KIND=GbcsInt2) :: fill
    TYPE(NEIGHBOURHOOD_MAP) :: nnmap, nnmap2
    REAL :: scale, offset
    REAL, ALLOCATABLE, DIMENSION(:,:) :: rads, stdd
    LOGICAL :: calc_std

    WRITE(*,'(A,I2,A,A)') "Calculating ", nnear, " nearest neighbours: ", view
    CALL compute_scene_neighbourhood(path, view, 'a', nnmap)
    IF (use_bstripe) THEN
      WRITE(*,*) "Adding b-stripe"
      CALL compute_scene_neighbourhood(path, view, 'b', nnmap2)
      CALL merge_neighbourhoods(nnmap, nnmap2)
    END IF
    irsize = SHAPE(nnmap%entries)
    WRITE(*,'(A,I4," x ",I4)') " Output size: ", irsize
    ALLOCATE(rads(irsize(1), irsize(2)))

    DO iband = 1,6
      ! Do we need to write subpixel texture for this channel
      calc_std = ANY(channels == iband)
      IF (calc_std .AND. .NOT. ALLOCATED(stdd)) ALLOCATE(stdd(irsize(1), irsize(2)))
      IF (ALLOCATED(stdd) .AND. .NOT. calc_std) DEALLOCATE(stdd)

      ! Read attributes from input radiance
      WRITE(name, '("S", I1, "_radiance_a", A1)') iband, view
      status = nf90_open(Path_Join(path, TRIM(name)//'.nc'), NF90_NOWRITE, ncid)
      status = nf90_inq_varid(ncid, name, v_rad)
      status = nf90_get_att(ncid, v_rad, '_FillValue', fill)
      status = nf90_get_att(ncid, v_rad, 'add_offset', offset)
      status = nf90_get_att(ncid, v_rad, 'scale_factor', scale)
      status = nf90_get_att(ncid, v_rad, 'long_name', attr_long)
      status = nf90_get_att(ncid, v_rad, 'standard_name', attr_name)
      status = nf90_get_att(ncid, v_rad, 'units', attr_unit)
      status = nf90_close(ncid)

      WRITE(name, '("S", I1, "_radiance_i", A1)') iband, view
      IF (use_bstripe .AND. iband > 3) THEN
        CALL process_scene_band(path, view, iband, rads, nnmap2, nnear*2, std=stdd)
      ELSE
        CALL process_scene_band(path, view, iband, rads, nnmap, nnear, std=stdd)
      END IF
      status = nf90_create(Path_Join(path, TRIM(name)//'.nc'), NF90_NETCDF4, ncid)
      status = nf90_def_dim(ncid, "columns", irsize(1), dim1)
      status = nf90_def_dim(ncid, "rows", irsize(2), dim2)

      ! Create radiance variable
      status = nf90_def_var(ncid, name, nf90_short, [dim1, dim2], v_rad, &
                            chunksizes=irsize, deflate_level=2, shuffle=.True.)
      status = nf90_put_att(ncid, v_rad, '_FillValue', fill)
      status = nf90_put_att(ncid, v_rad, 'add_offset', offset)
      status = nf90_put_att(ncid, v_rad, 'scale_factor', scale)
      status = nf90_put_att(ncid, v_rad, 'valid_min', -3000_GbcsInt2)
      status = nf90_put_att(ncid, v_rad, 'valid_max', 32767_GbcsInt2)
      status = nf90_put_att(ncid, v_rad, 'long_name', attr_long)
      status = nf90_put_att(ncid, v_rad, 'standard_name', attr_name)
      status = nf90_put_att(ncid, v_rad, 'units', attr_unit)

      ! Create radiance standard deviation variable
      IF (calc_std) THEN
        WRITE(name, '("S", I1, "_stddev_i", A1)') iband, view
        status = nf90_def_var(ncid, name, nf90_short, [dim1, dim2], v_std, &
                              chunksizes=irsize, deflate_level=2, shuffle=.True.)
        status = nf90_put_att(ncid, v_std, '_FillValue', fill)
        status = nf90_put_att(ncid, v_std, 'add_offset', 0.0)
        status = nf90_put_att(ncid, v_std, 'scale_factor', scale)
        status = nf90_put_att(ncid, v_std, 'valid_min', 0_GbcsInt2)
        status = nf90_put_att(ncid, v_std, 'valid_max', 32767_GbcsInt2)
        status = nf90_put_att(ncid, v_std, 'long_name', "Standard deviation of radiance observations")
        status = nf90_put_att(ncid, v_std, 'cell_methods', 'area: standard_deviation')
        status = nf90_put_att(ncid, v_std, 'units', attr_unit)
      END IF

      WRITE(attr_long, '(A,I2,A)') "Regrid to i-stripe using ", nnear, " nearest neighbours."
      IF (use_bstripe .AND. iband > 3) THEN
        attr_long = TRIM(attr_long) // " Including b-stripe."
      END IF
      status = nf90_put_att(ncid, nf90_global, 's3regrid', attr_long)
      status = nf90_enddef(ncid)

      ! Output
      CALL Pack_Data(rads, offset, scale, -3000, 32767, INT(fill))
      status = nf90_put_var(ncid, v_rad, rads)

      IF (calc_std) THEN
        CALL Pack_Data(stdd, 0.0, scale, 0, 32767, INT(fill))
        status = nf90_put_var(ncid, v_std, stdd)
      END IF

      status = nf90_close(ncid)
    END DO

  END SUBROUTINE process_view

  PURE SUBROUTINE Pack_Data(values, offset, scale, vmin, vmax, fill)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    REAL, DIMENSION(:,:), INTENT(inout) :: values
    REAL,                 INTENT(in)    :: offset
    REAL,                 INTENT(in)    :: scale
    INTEGER,              INTENT(in)    :: vmin
    INTEGER,              INTENT(in)    :: vmax
    INTEGER,              INTENT(in)    :: fill

    values = ANINT((values - offset) / scale)
    WHERE(values < vmin .OR. values > vmax)
      values = fill
    END WHERE

  END SUBROUTINE Pack_Data

END PROGRAM s3regrid
