!--------------------------------------------------------------------------------------------
!P+ Preprocess_SLSTR
! NAME:
!     Preprocess_SLSTR
!
! PURPOSE:
!
! Test driver program for SLSTR_Preprocessor
!
! Preprocess radiance band files 1-6 in nadir and oblique views from scene folder provided as first command line argument
! and write resulting reprocessed band files (retaining the same name) to the output folder provided as the second
! command line argument.
!
! By default the preprocessor calculates neighbourhoods based on the across- and along-track distances. If a third
! command line argument "simple" is supplied, fall back to an existing, simple but fast method.
!
! CATEGORY:
!     Data I/O
!
! LANGUAGE:
!     Fortran-95
!
! PUBLIC DATA:
!     None
!
! MODULES:
!
! CONTAINS:
!    process_view
!
! NOTES:
!
! CREATION HISTORY:
!     Written by:  Niall McCarroll 12/01/21
!                University of Reading
!
!  Copyright 2020 Chris Merchant, Owen Embury, and Claire Bulgin
!  The Institute of Atmospheric and Environmental Science
!  The University of Edinburgh, UK
!
!  This file is part of the GBCS software package.
!
!  The GBCS software package is free software: you can redistribute it
!  and/or modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation, either version 3 of
!  the License, or (at your option) any later version.
!
!  The GBCS software package is distributed in the hope that it will be
!  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with the GBCS software package.
!  If not, see <http://www.gnu.org/licenses/>.
!
!P-
!--------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  !S+
  ! NAME:
  !     safe_write_out
  !
  ! PURPOSE:
  !>    Write out an array to a netcdf4 file
  !
  ! CATEGORY:
  !
  ! LANGUAGE:
  !     Fortran-95
  !
  ! CALLING SEQUENCE:
  !     safe_write_out('/path/to/file.nc','field1',1000,1000,arr,-99999)
  !
  ! ARGUMENTS:
  !>@ARG{path, in, CHARACTER(LEN=*)} path of the file to be written / overwritten
  !>@ARG{field_name, in, CHARACTER(LEN=*)} the name of the variable to write into the file
  !>@ARG{w, in, INTEGER} the width of the array
  !>@ARG{h, in, INTEGER} the height of the array
  !>@ARG{v, in, INTEGER} the number of variables in te array
  !>@ARG{values, in, REAL\, DIMENSION(w\,h\,v)} array of real values to write
  !>@ARG{fill_value, in, REAL} a value to represent missing values in the array
  !>@ARG{suffixes, in, REAL\, DIMENSION(v)} array of suffixes to use for each variable written
  !
  ! CALLS:
  !     nf90_create, nf90_def_dim, nf90_def_var, nf90_put_att, nf90_enddef, nf90_put_var, nf90_close
  !
  ! SIDE EFFECTS:
  !     Stops program if code indicates a netcdf module operation has failed
  !
  ! RESTRICTIONS:
  !
  ! PROCEDURE:
  !
  ! CREATION HISTORY:
  !     11/11/20  NM  Creation
  !
  !S-
  !------------------------------------------------------------------------------
  SUBROUTINE safe_write_out(path,field_name,w,h,v,values,fill_value,suffixes)
    USE netcdf
    USE SLSTR_Preprocessor
    IMPLICIT NONE

    ! -----------
    ! Arguments
    ! -----------
    CHARACTER(LEN=*),  INTENT(IN)  :: path, field_name
    INTEGER, INTENT(IN) :: w, h, v
    REAL, DIMENSION(w,h,v), INTENT(IN) :: values
    REAL, INTENT(IN) :: fill_value
    CHARACTER(LEN=*), DIMENSION(v), INTENT(IN) :: suffixes

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: ncid
    INTEGER :: status
    INTEGER :: x_dimid, y_dimid, varid
    INTEGER, DIMENSION(2) :: dimids
    INTEGER :: var_index
    REAL, DIMENSION(w,h) :: value_slice
    INTEGER, DIMENSION(v) :: varids

    status = nf90_create(path, NF90_CLOBBER, ncid)
    CALL handle_err(status)

    status = nf90_def_dim(ncid, "columns", w, x_dimid)
    CALL handle_err(status)

    status = nf90_def_dim(ncid, "rows", h, y_dimid)
    CALL handle_err(status)

    dimids =  (/ x_dimid, y_dimid /)

    DO var_index=1,v
      status = nf90_def_var(ncid, TRIM(field_name)//'_'//TRIM(suffixes(var_index)), NF90_REAL4, dimids, varid)
      CALL handle_err(status)

      status = nf90_put_att(ncid, varid, "_FillValue", fill_value)
      CALL handle_err(status)
      varids(var_index) = varid
    END DO

    status = nf90_enddef(ncid)
    CALL handle_err(status)

    DO var_index=1,v
      value_slice = values(:,:,var_index)
      status = nf90_put_var(ncid, varids(var_index), value_slice)
      CALL handle_err(status)
    END DO

    status = nf90_close(ncid)
    CALL handle_err(status)

  END SUBROUTINE safe_write_out


!------------------------------------------------------------------------------
!S+
! NAME:
!     show_neighbourhood_stats
!
! PURPOSE:
!     dump statistics on a computed neighbourhood, and perform some integrity checks on the neighbourhood.
!
! CATEGORY:
!
! LANGUAGE:
!     Fortran-95
!
! CALLING SEQUENCE:
!
! ARGUMENTS:
!     neighbourhood   - a computed array of computed neighbourhoods, organised by (COL,ROW)
!
! CALLS:
!
! SIDE EFFECTS:
!     Statistics are written to stdout
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     15/01/21  NM  Creation
!
!S-
!------------------------------------------------------------------------------
SUBROUTINE show_neighbourhood_stats(neighbourhood, title)
  USE SLSTR_Preprocessor
  IMPLICIT NONE

  ! -----------
  ! Arguments
  ! -----------
  TYPE(NEIGHBOURHOOD_MAP), INTENT(IN) :: neighbourhood
  CHARACTER(LEN=*), INTENT(IN) :: title

  ! ---------------
  ! Local Variables
  ! ---------------
  REAL dist_from_expected, x, y, sqdist
  INTEGER :: distance_bucket, distance_gt_10, distance_gt_10_orphan
  INTEGER :: row, col, k, rows, cols, x_index, y_index
  INTEGER, DIMENSION(10) :: distance_distribution
  INTEGER, DIMENSION(10) :: distance_distribution_orphan
  INTEGER, DIMENSION(4) :: window_shrink_stats_x, window_shrink_stats_y
  INTEGER :: i, shrink
  INTEGER :: main_count_a, orphan_count_a, main_count_b, orphan_count_b, null_count
  INTEGER :: source
  INTEGER :: y_min_index, y_max_index, x_min_index, x_max_index

  main_count_a = 0
  orphan_count_a = 0
  main_count_b = 0
  orphan_count_b = 0
  null_count = 0
  distance_distribution(:) = 0
  distance_distribution_orphan(:) = 0
  distance_gt_10 = 0
  distance_gt_10_orphan = 0
  window_shrink_stats_x(:) = 0
  window_shrink_stats_y(:) = 0

  cols = SIZE(neighbourhood%entries,1)
  rows = SIZE(neighbourhood%entries,2)

  ! print some handy statistics on the sources of pixels in the neighbourhood
  DO row=1,rows
    DO col=1,cols
      sqdist=0
      DO k=1,MAX_K_NEAREST_NEIGHBOURS
        source = neighbourhood%entries(col,row)%source(k)
        if (source == ORPHAN_PIXEL_SOURCE_A .or. source == ORPHAN_PIXEL_SOURCE_B) THEN
          IF (source == ORPHAN_PIXEL_SOURCE_A) THEN
            orphan_count_a = orphan_count_a + 1
          ELSE
            orphan_count_b = orphan_count_b + 1
          END IF
          y = neighbourhood%entries(col,row)%y(k)
          ! track in 2D how far (in pixels) the neighbourhood pixel is from the idealized "2x" location
          distance_bucket = ABS(y-2*row)+1
          IF (distance_bucket < 11) THEN
            distance_distribution_orphan(distance_bucket) = 1+distance_distribution_orphan(distance_bucket)
          ELSE
            distance_gt_10_orphan = distance_gt_10_orphan + 1
          END IF
        END IF
        if (source == MAIN_PIXEL_SOURCE_A .or. source == MAIN_PIXEL_SOURCE_B) THEN
          IF (source == MAIN_PIXEL_SOURCE_A) THEN
            main_count_a = main_count_a + 1
          ELSE
            main_count_b = main_count_b + 1
          END IF
          x = neighbourhood%entries(col,row)%x(k)
          y = neighbourhood%entries(col,row)%y(k)
          ! track in 1D how far (in pixels) the neighbourhood pixel is from the idealized "2x" location
          dist_from_expected = SQRT((x-(2*col))**2 + (y-(2*row))**2)
          distance_bucket = NINT(dist_from_expected)+1
          IF (distance_bucket < 11) THEN
            distance_distribution(distance_bucket) = 1+distance_distribution(distance_bucket)
          ELSE
            distance_gt_10 = distance_gt_10 + 1
          END IF
        END IF
        if (source == NULL_PIXEL_SOURCE) THEN
          null_count = null_count + 1
        ELSE
          ! run integrity checks.  do the neighbourhhod pixels lie in the expected search window?
          IF (.false.) THEN
            x_index = neighbourhood%entries(col,row)%x(k)
            y_index = neighbourhood%entries(col,row)%y(k)

            y_min_index = MAX(2*row - search_height,1)
            y_max_index = MIN(2*row + search_height, 2*rows)
            x_min_index = 1
            x_max_index = cols*2

            IF (source == MAIN_PIXEL_SOURCE_A .or. source == MAIN_PIXEL_SOURCE_B) THEN
              x_min_index = MAX(2*col - search_width, 1)
              x_max_index = MIN(2*col + search_width, 2*cols)
            END IF

            IF (x_index < x_min_index .or. x_index > x_max_index) THEN
              PRINT *, 'ERROR - x index in neighbourhood is out of bounds'
              PRINT *, x_index
              STOP
            END IF
            IF (y_index < y_min_index .or. y_index > y_max_index) THEN
              PRINT *, 'ERROR - y index in neighbourhood is out of bounds'
              PRINT *, y_index
              STOP
            END IF

            DO shrink = 1,4
              IF (x_index < (x_min_index+shrink) .or. x_index > (x_max_index-shrink)) THEN
                window_shrink_stats_x(shrink) = window_shrink_stats_x(shrink)+1
              END IF
              IF (y_index < (y_min_index+shrink) .or. y_index > (y_max_index-shrink)) THEN
                window_shrink_stats_y(shrink) = window_shrink_stats_y(shrink)+1
              END IF
            END DO

            ! values at each entry should be ordered by increasing squared distance, check this
            IF (sqdist > neighbourhood%entries(col,row)%squared_distances(k)) THEN
              PRINT *, 'ERROR - distances in neighbourhood are out of order'
              PRINT *, neighbourhood%entries(col,row)%squared_distances
              STOP
            END IF
            sqdist = neighbourhood%entries(col,row)%squared_distances(k)

            IF (sqdist > MAX_NEIGHBOUR_DISTANCE*MAX_NEIGHBOUR_DISTANCE) THEN
              PRINT *, 'ERROR - distances in neighbourhood exceeds maximum'
              STOP
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO
  PRINT *, ''
  PRINT *, '***', title, '***'
  PRINT *, ''
  PRINT *, '  MAIN PIXELS USED(A)=', main_count_a
  PRINT *, '  MAIN PIXELS USED(B)=', main_count_b
  PRINT *, 'ORPHAN PIXELS USED(A)=', orphan_count_a
  PRINT *, 'ORPHAN PIXELS USED(B)=', orphan_count_b
  PRINT *, '  NULL PIXELS USED   =', null_count
  PRINT *, ''
  PRINT *, 'SHRINKAGE STATS'
  DO i=1,4
    PRINT *, i, window_shrink_stats_x(i), window_shrink_stats_y(i)
  END DO
  PRINT *, ''
  PRINT *, 'MAIN PIXEL LOCATION DISTRIBUTION OF DISTANCES FROM EXPECTED'
  DO i=1,10
    PRINT *, i, distance_distribution(i)
  END DO
  PRINT *, '        >10', distance_gt_10
  PRINT *, ''
  PRINT *, 'ORPHAN PIXEL LOCATION DISTRIBUTION OF DISTANCES FROM EXPECTED'
  DO i=1,10
    PRINT *, i, distance_distribution_orphan(i)
  END DO
  PRINT *, '        >10', distance_gt_10_orphan
END SUBROUTINE show_neighbourhood_stats

!------------------------------------------------------------------------------
!S+
! NAME:
!     process_view
!
! PURPOSE:
!     given a folder containing an SLSTR scene and a view (nadir or oblique), loads the cartesian
!     locations for visible and ir pixel grids and computes the neighborhood of visible pixels closest to
!     each IR pixel.  The neigbourhood can then be used in calls to process_scene_band to map visible
!     bands from the same scene and view to a new IR based grid.
!
!     Note that the neighbourhood is defined according to the parameters and source code in the
!     slstr_preprocessor module
!
! CATEGORY:
!
! LANGUAGE:
!     Fortran-95
!
! CALLING SEQUENCE:
!   process_view('o','/path/to/scene','/path/to/output/folder',.false.,.false.,(\'mean'\))
!
! ARGUMENTS:
!     view_type       - 'o' for oblique view or 'n' for nadir view
!     scene_folder    - path to a folder storing the various input files from an SLSTR scene
!     output_folder   - path to a folder to store the various output files produced by the preprocessor
!     simple          - .true. iff the processing should use the simpler original method
!     stats           - .true. iff the processing should print neighbourhood stats
!
! CALLS:
!     Invokes compute_scene_neighbourhood and process_scene_band in module slstr_preprocessor
!
! SIDE EFFECTS:
!     Output files are written to the output_folder
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     11/11/20  NM  Creation
!
!S-
!------------------------------------------------------------------------------
SUBROUTINE process_view(view_type, scene_folder, output_folder, simple, stats, compare, effective_k)
  USE SLSTR_Preprocessor
  USE GbcsPath
  USE netcdf
  IMPLICIT NONE

  ! -----------
  ! Arguments
  ! -----------
  CHARACTER(1), INTENT(IN) :: view_type
  CHARACTER(LEN=*), INTENT(IN) :: scene_folder
  CHARACTER(LEN=*), INTENT(IN) :: output_folder
  LOGICAL, INTENT(IN) :: simple, stats, compare
  INTEGER, INTENT(IN) :: effective_k

  ! ---------------
  ! Local Variables
  ! ---------------
  INTEGER :: ir_width, ir_height, band, rad_ncid, status

  CHARACTER(1) :: band_str
  CHARACTER(256) :: output_field_name
  TYPE(NEIGHBOURHOOD_MAP) :: neighbourhood_a, neighbourhood_ab
  REAL, ALLOCATABLE, DIMENSION(:,:) :: vis_output_radiance
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: vis_output_radiances
  REAL :: fill_value
  INTEGER :: function_index
  INTEGER, DIMENSION(4) :: functions
  CHARACTER(16), DIMENSION(4) :: function_names

  REAL :: start_time, end_time, elapsed_time

  CALL cpu_time(start_time)

  functions(1) = FUNCTION_MEAN
  function_names(1) = 'mean'

  functions(2) = FUNCTION_SD
  function_names(2) = 'sd'

  functions(3) = FUNCTION_MAX
  function_names(3) = 'max'

  functions(4) = FUNCTION_MIN_MAX_DIFF
  function_names(4) = 'min_max_diff'

  ir_height = 1200
  IF (view_type == 'n') THEN
    ! nadir view
    ir_width = 1500
  ELSE
    ! oblique view
    ir_width = 900
  END IF

  ALLOCATE(vis_output_radiance(ir_width,ir_height))
  ALLOCATE(vis_output_radiances(ir_width,ir_height,4))

  IF (.not. simple) THEN
    IF (stats) THEN
      IF (view_type == 'n') THEN
        PRINT *, 'Nadir View'
      ELSE
        PRINT *, 'Oblique View'
      END IF
      PRINT *, 'max neighbourhood size', MAX_K_NEAREST_NEIGHBOURS
      PRINT *, 'effective neighbourhood size', effective_k
      PRINT *, 'max distance (m)', MAX_NEIGHBOUR_DISTANCE
    END IF
    CALL compute_scene_neighbourhood(view_type,scene_folder,neighbourhood_a,'a')
    IF (stats) THEN
      CALL show_neighbourhood_stats(neighbourhood_a, &
              'Neighbourhood A stats - view='//view_type)
    END IF
    CALL compute_scene_neighbourhood(view_type,scene_folder,neighbourhood_ab,'b')
    IF (stats) THEN
      CALL show_neighbourhood_stats(neighbourhood_ab, &
              'Neighbourhood B stats - view='//view_type)
    END IF
    CALL merge_neighbourhoods(neighbourhood_a,neighbourhood_ab)
    IF (stats) THEN
      CALL show_neighbourhood_stats(neighbourhood_ab, &
              'Neighbourhood A+B stats - view='//view_type)
    END IF
  END IF

  rad_ncid = safe_open(TRIM(scene_folder)//'/'//'S1_radiance_a'//view_type//'.nc')
  fill_value = safe_get_real_attribute(rad_ncid,'S1_radiance_a'//view_type,'_FillValue')
  status = nf90_close(rad_ncid)

  DO band = 1,6
    WRITE(band_str,'(I1)') band
    WRITE(output_field_name, '("S",I1,"_radiance_i",A1)') band, view_type

    DO function_index = 1,4
      IF (band > 3) THEN
        ! bands 4,5,6 should have data avaialble for both a and b stripes
        CALL process_scene_band(view_type,scene_folder,band,vis_output_radiance, &
                neighbourhood_ab,effective_k,functions(function_index))
      ELSE
        ! bands 3,4,5 may only have the a stripe
        CALL process_scene_band(view_type,scene_folder,band,vis_output_radiance, &
                neighbourhood_a,effective_k,functions(function_index))
      END IF
      vis_output_radiances(:,:,function_index) = vis_output_radiance
    END DO

    WHERE (vis_output_radiances == MISSING_R)
      vis_output_radiances = fill_value
    END WHERE

    CALL safe_write_out(Path_Join(output_folder,'S'//band_str//'_radiance_i'//view_type//'.nc'), &
      output_field_name,ir_width,ir_height,4,vis_output_radiances,fill_value,function_names)

  END DO
  CALL cpu_time(end_time)
    elapsed_time = end_time - start_time

    IF (view_type == 'n') THEN
      PRINT *, 'Nadir View processing CPU time (secs) ', elapsed_time
    ELSE
      PRINT *, 'Oblique View processing CPU time (secs) ', elapsed_time
    END IF
END SUBROUTINE process_view

PROGRAM Preprocess_SLSTR
  USE SLSTR_Preprocessor
  IMPLICIT NONE

  ! ---------------
  ! Local Variables
  ! ---------------
  CHARACTER(256) :: scene_folder
  CHARACTER(256) :: output_folder

  CHARACTER(256) :: option
  CHARACTER(256) :: option_value

  INTEGER :: effective_k, i, max_distance, window_height, window_width
  LOGICAL :: simple, stats, compare

  simple = .false.
  stats = .false.
  compare = .false.
  effective_k = MAX_K_NEAREST_NEIGHBOURS
  max_distance = 10000
  window_height = 0
  window_width = 0

  CALL GET_COMMAND_ARGUMENT(1,scene_folder)
  CALL GET_COMMAND_ARGUMENT(2,output_folder)

  i = 3
  DO WHILE (i <= COMMAND_ARGUMENT_COUNT())
    CALL GET_COMMAND_ARGUMENT(i,option)
    i = i + 1
    IF (option == '--simple') THEN
      simple = .true.
    ELSE IF (option == '--stats') THEN
      stats = .true.
    ELSE IF (option == '--compare') THEN
      compare = .true.
    ELSE IF (option == '--effective_k') THEN
      CALL GET_COMMAND_ARGUMENT(i,option_value)
      IF (option_Value == '') THEN
        EXIT
      END IF
      i = i + 1
      READ(option_value,'(I3.1)') effective_k
    ELSE IF (option == '--window_height') THEN
      CALL GET_COMMAND_ARGUMENT(i,option_value)
      IF (option_Value == '') THEN
        EXIT
      END IF
      i = i + 1
      READ(option_value,'(I3.1)') window_height
    ELSE IF (option == '--window_width') THEN
      CALL GET_COMMAND_ARGUMENT(i,option_value)
      IF (option_Value == '') THEN
        EXIT
      END IF
      i = i + 1
      READ(option_value,'(I3.1)') window_width
    ELSE IF (option == '--max_distance') THEN
      CALL GET_COMMAND_ARGUMENT(i,option_value)
      IF (option_Value == '') THEN
        EXIT
      END IF
      i = i + 1
      READ(option_value,'(I6.1)') max_distance
    ELSE
      PRINT *, 'Unknown option ', option
      STOP
    END IF
  END DO

  MISSING_R = -1.0e+30
  MAX_NEIGHBOUR_DISTANCE = max_distance
  IF (window_height > 0) THEN
    search_height = window_height
  END IF
  IF (window_width > 0) THEN
    search_width = window_width
  END IF

  CALL process_view('n',scene_folder,output_folder,simple,stats,compare,effective_k)
  CALL process_view('o',scene_folder,output_folder,simple,stats,compare,effective_k)
END PROGRAM Preprocess_SLSTR
