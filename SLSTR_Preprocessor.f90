!------------------------------------------------------------------------------
!M+
! NAME:
!       SLSTR_Preprocessor
!
! PURPOSE:
!>      Regrid Vis and NIR observations at 1km using nearest neighbourhood and
!>      including orphan pixels
!
! CATEGORY:
!       Data I/O
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       USE SLSTR_Preprocessor
!
! PUBLIC DATA:
!       None
!
! MODULES:
!       GbcsPath
!
! CONTAINS:
!       handle_err
!       safe_open
!       safe_get_real_data
!       safe_get_real_attribute
!       safe_get_int_data
!       safe_get_dimlen
!       insert_neighbour
!       find_neighbours
!       build_neighbourhood_map
!       apply_function
!       apply_neighbours
!       apply_simple_aggregation
!       load_radiance_data
!       process_scene_band
!       align_cosmetics
!       merge_neighbourhoods
!       compute_scene_neighbourhood
!
! DERIVED TYPES:
!       NEIGHBOURHOOD_ENTRY
!       NEIGHBOURHOOD_MAP
!       Locations
!
! NOTES:
!
! CREATION HISTORY:
!     Written by:  Niall McCarroll 12/01/21
!                University of Reading
!
!  Copyright 2021 Owen Embury, and Claire Bulgin
!  Department of Meteorology, University of Reading, UK
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
!M-
!--------------------------------------------------------------------------------------------
MODULE SLSTR_Preprocessor
  IMPLICIT NONE

  ! ----------------
  ! Module variables
  ! ----------------

  !> Real number to represent missing values in I/O. 
  REAL :: MISSING_R = -1.0e+30

  !> Integer value to represent missing values in I/O.
  INTEGER :: MISSING_I = -32768
  
  !> Exclude neighbours further away than this (metres).
  REAL :: MAX_NEIGHBOUR_DISTANCE = 10000

  !> Define a search window around the area (2*x,2*y) in which to look for closest visible pixels
  INTEGER :: search_width = 8, search_height = 6     ! Search window size, in pixels


  ! -----------------
  ! Module parameters
  ! -----------------

  !> Whether to perform extra computation to process cosmetically filled IR pixels the same as the original pixels
  LOGICAL :: align_cosmetic_pixels = .false.

  !> Maximum neighbourhood size.
  INTEGER, PARAMETER :: MAX_K_NEAREST_NEIGHBOURS = 15

  !> Mask to use to check if a pixel is cosmetically filled, in flags_ao.nc/flags_an.nc
  INTEGER, PARAMETER :: COSMETIC_PIXEL_MASK = 256      ! Identify COSMETIC pixels with these bits flagged

  !> internal: code for an invalid pixel index
  INTEGER, PARAMETER :: INVALID_PIXEL = -1

  !> internal: code for the source of a neighbourhood pixel from main dataset (stripe A)
  INTEGER, PARAMETER :: MAIN_PIXEL_SOURCE_A = 0

  !> internal: code for the source of a neighbourhood pixel from orphan dataset (stripe A)
  INTEGER, PARAMETER :: ORPHAN_PIXEL_SOURCE_A = 1

  !> internal: code for the source of a neighbourhood pixel from main dataset (stripe B)
  INTEGER, PARAMETER :: MAIN_PIXEL_SOURCE_B = 2

  !> internal: code for the source of a neighbourhood pixel from orphan dataset (stripe B)
  INTEGER, PARAMETER :: ORPHAN_PIXEL_SOURCE_B = 3

  !> Allowed confidence code for visible pixels
  INTEGER, PARAMETER :: EXCEPTION_SATURATION = 16

  !> Constant representing the mean function
  INTEGER, PARAMETER :: FUNCTION_MEAN = 1

  !> Constant representing the std deviation function
  INTEGER, PARAMETER :: FUNCTION_SD = 2

  !> Constant representing the max function
  INTEGER, PARAMETER :: FUNCTION_MAX = 3

  !> Constant representing the max-min function
  INTEGER, PARAMETER :: FUNCTION_MIN_MAX_DIFF = 4

  ! -------------
  ! Derived Types
  ! -------------

  !> Define a structure grouping three arrays, organised by K into a neighbourhood
  TYPE NEIGHBOURHOOD_ENTRY
    INTEGER, DIMENSION(MAX_K_NEAREST_NEIGHBOURS) :: x = INVALID_PIXEL   !< Image coordinate
    INTEGER, DIMENSION(MAX_K_NEAREST_NEIGHBOURS) :: y = INVALID_PIXEL   !< Image coordinate
    INTEGER, DIMENSION(MAX_K_NEAREST_NEIGHBOURS) :: source = INVALID_PIXEL  !< Pixel source (a, b, orphan etc.)
    REAL,    DIMENSION(MAX_K_NEAREST_NEIGHBOURS) :: d = 0   !< Distance
    INTEGER                                      :: n = 0   !< Number of pixels in neighbourhood
  END TYPE NEIGHBOURHOOD_ENTRY

  !> Define a neighbourhood as an array of entries and metadata on which stripe(s) the neighbourhood is built on
  TYPE NEIGHBOURHOOD_MAP
    LOGICAL :: include_a_stripe
    LOGICAL :: include_b_stripe
    TYPE(NEIGHBOURHOOD_ENTRY), ALLOCATABLE, DIMENSION(:,:) :: entries
  END TYPE NEIGHBOURHOOD_MAP

  !> Define a structure to hold the pixel across (x) and along (y) track distances for each (col,row)
  TYPE Locations
    INTEGER :: width = 0
    INTEGER :: height = 0
    REAL,    ALLOCATABLE, DIMENSION(:,:) :: x     !< Across track pixel coordinate
    REAL,    ALLOCATABLE, DIMENSION(:,:) :: y     !< Along track pixel coordinate
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: conf  !< Confidence flags
  END TYPE Locations

CONTAINS
!------------------------------------------------------------------------------
!S+
! NAME:
!       handle_err
!
! PURPOSE:
!>      Passed a status code returned from a netcdf module operation, if the code
!>      indicates an error, print more information and stop the program
!
! CATEGORY:
!       netCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL handle_err(code)
!
! ARGUMENTS:
!>@ARG{status, in, INTEGER} netCDF return status
!
! CALLS:
!       nf90_strerror
!
! SIDE EFFECTS:
!       Stops program if code indicates a netcdf module operation has failed
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!       11/11/20  NM  Creation
!S-
!------------------------------------------------------------------------------
  SUBROUTINE handle_err(status)
    USE netcdf
    IMPLICIT NONE
    ! -----------
    ! Arguments
    ! -----------
    INTEGER, INTENT (IN) :: status

    IF (status /= nf90_noerr) THEN
      PRINT *, nf90_strerror(status)
      STOP "Stopped"
    END IF
  END SUBROUTINE handle_err


!------------------------------------------------------------------------------
!F+
! NAME:
!       safe_open
!
! PURPOSE:
!>      Open a netcdf4 file and return its fileid
!
! CATEGORY:
!       netCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       ncid = safe_open('/path/to/file.nc')
!
! ARGUMENTS:
!>@ARG{path, in, CHARACTER(LEN=*)} path of the file to be opened
!
! FUNCTION RESULT:
!>@RES{id, INTEGER} ID of opened file
!
! CALLS:
!       nf90_open
!
! SIDE EFFECTS:
!       Stops program if code indicates a netcdf module operation has failed
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!       11/11/20  NM  Creation
!F-
!------------------------------------------------------------------------------
  FUNCTION safe_open(path)
    USE netcdf
    IMPLICIT NONE
    ! -----------
    ! Arguments
    ! -----------
    CHARACTER(LEN=*),  INTENT(IN)  :: path

    ! ----------------
    ! Function Result
    ! ----------------
    INTEGER :: safe_open

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: ncid
    INTEGER :: status

    status = nf90_open(path, NF90_NOWRITE, ncid)
    CALL handle_err(status)
    safe_open = ncid
  END FUNCTION safe_open


!------------------------------------------------------------------------------
!S+
! NAME:
!       safe_get_real_data
!
! PURPOSE:
!>      Reads real values from an opened netcdf4 file into an allocated array
!
! CATEGORY:
!       netCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL safe_get_real_data(ncid, name, values, fill_value)
!
! ARGUMENTS:
!>@ARG{ncid, in, INTEGER} ID of the opened file (returned from safe_open)
!>@ARG{name, in, CHARACTER (LEN = *)} name of the variable
!>@ARG{values, out, REAL\, DIMENSION(:\,:)} array of real values which will be read, organised by (COL,ROW)
!>@ARG{fill_value, in, REAL} a real value to represent missing values in the array
!
! CALLS:
!       nf90_inq_varid
!       nf90_get_var
!       nf90_get_att
!
! SIDE EFFECTS:
!       Stops program if code indicates a netcdf module operation has failed
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!       11/11/20  NM  Creation
!S-
!------------------------------------------------------------------------------
  SUBROUTINE safe_get_real_data(ncid, name, values, fill_value)
    USE netcdf
    IMPLICIT NONE
    ! -----------
    ! Arguments
    ! -----------
    INTEGER,              INTENT(in)  :: ncid
    CHARACTER(LEN=*),     INTENT(in)  :: name
    REAL, DIMENSION(:,:), INTENT(out) :: values
    REAL,                 INTENT(in)  :: fill_value

    ! ---------------
    ! Local Variables
    ! ---------------
    REAL    :: file_fill_value
    INTEGER :: varid
    INTEGER :: status
    REAL    :: scale_factor

    status = nf90_inq_varid(ncid,name,varid)
    CALL handle_err(status)
    status = nf90_get_var(ncid,varid,values)
    CALL handle_err(status)
    status = nf90_get_att(ncid,varid,"scale_factor",scale_factor)
    CALL handle_err(status)
    status = nf90_get_att(ncid,varid,"_FillValue",file_fill_value)
    CALL handle_err(status)
    WHERE (values /= file_fill_value)
      values = values * scale_factor
    ELSE WHERE
      values = fill_value
    END WHERE

  END SUBROUTINE safe_get_real_data


!------------------------------------------------------------------------------
!F+
! NAME:
!       safe_get_real_attribute
!
! PURPOSE:
!>      Reads real value for named attribute from field in an opened netcdf4 file
!
! CATEGORY:
!       netCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       value = safe_get_real_attribute(ncid, field_name, attr_name)
!
! ARGUMENTS:
!>@ARG{ncid, in, INTEGER} ID of the opened file (returned from safe_open)
!>@ARG{field_name, in, CHARACTER(LEN=*)} name of the variable
!>@ARG{attr_name, in, CHARACTER(LEN=*)} name of the attribute
!
! RETURNS:
!@RES{value, REAL} value of attribute
!
! CALLS:
!       nf90_inq_varid
!       nf90_get_att
!
! SIDE EFFECTS:
!       Stops program if code indicates a netcdf module operation has failed
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     11/11/20  NM  Creation
!F-
!------------------------------------------------------------------------------
  FUNCTION safe_get_real_attribute(ncid, field_name, attr_name)
    USE netcdf
    IMPLICIT NONE
    ! -----------
    ! Arguments
    ! -----------
    INTEGER,          INTENT(in) :: ncid
    CHARACTER(LEN=*), INTENT(in) :: field_name
    CHARACTER(LEN=*), INTENT(in) :: attr_name

    ! ----------------
    ! Function Result
    ! ----------------
    REAL :: safe_get_real_attribute

    ! ---------------
    ! Local Variables
    ! ---------------
    REAL    :: value
    INTEGER :: varid
    INTEGER :: status

    status = nf90_inq_varid(ncid, field_name, varid)
    CALL handle_err(status)
    status = nf90_get_att(ncid, varid, attr_name, value)
    CALL handle_err(status)
    safe_get_real_attribute = value

  END FUNCTION safe_get_real_attribute


!------------------------------------------------------------------------------
!S+
! NAME:
!       safe_get_int_data
!
! PURPOSE:
!>      Reads integer values from an opened netcdf4 file into an allocated array
!
! CATEGORY:
!       netCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL safe_get_int_data(ncid, name, values, fill_value)
!
! ARGUMENTS:
!>@ARG{ncid, in, INTEGER} ID of the opened file (returned from safe_open)
!>@ARG{name, in, CHARACTER(LEN=*)} name of the variable
!>@ARG{values, out, INTEGER\, DIMENSION(:\,:)} array of integer values which will be read, organised by (COL,ROW)
!>@ARG{fill_value, in, INTEGER} an integer value to represent missing values in the array
!
! CALLS:
!       nf90_inq_varid
!       nf90_get_var
!       nf90_get_att
!
! SIDE EFFECTS:
!       Stops program if code indicates a netcdf module operation has failed
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     11/11/20  NM  Creation
!S-
!------------------------------------------------------------------------------
  SUBROUTINE safe_get_int_data(ncid, name, values, fill_value)
    USE netcdf
    IMPLICIT NONE
    ! -----------
    ! Arguments
    ! -----------
    INTEGER,                  INTENT(in)  :: ncid
    CHARACTER(LEN=*),         INTENT(in)  :: name
    INTEGER, DIMENSION(:,:),  INTENT(out) :: values
    INTEGER,                  INTENT(in)  :: fill_value

    ! ---------------
    ! Local Variables
    ! ---------------

    INTEGER :: varid
    INTEGER :: status
    INTEGER :: file_fill_value

    status = nf90_inq_varid(ncid,name,varid)
    CALL handle_err(status)
    status = nf90_get_var(ncid,varid,values)
    CALL handle_err(status)
    status = nf90_get_att(ncid, varid, "_FillValue", file_fill_value)
    IF (status == nf90_noerr) THEN
      WHERE (values == file_fill_value)
        values = fill_value
      END WHERE
    END IF

  END SUBROUTINE safe_get_int_data


!------------------------------------------------------------------------------
!F+
! NAME:
!       safe_get_dimlen
!
! PURPOSE:
!>      Get the length of the specified dimension
! CATEGORY:
!       netCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       dimlen = safe_get_dimlen(ncid, name)
!
! ARGUMENTS:
!>@ARG{ncid, in, INTEGER} ID of the opened file (returned from safe_open)
!>@ARG{name, in, CHARACTER(LEN=*)} name of the variable
!
! CALLS:
!       nf90_inq_dimid
!       nf90_inquire_dimension
!F-
!------------------------------------------------------------------------------
  FUNCTION safe_get_dimlen(ncid, name)
    USE netcdf
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    INTEGER,          INTENT(in) :: ncid
    CHARACTER(LEN=*), INTENT(in) :: name

    ! ---------
    ! Function result
    ! ---------
    INTEGER                      :: safe_get_dimlen

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: status, dimid

    safe_get_dimlen = -1

    status = nf90_inq_dimid(ncid, name, dimid)
    IF (status /= nf90_noerr) RETURN

    status = nf90_inquire_dimension(ncid, dimid, LEN=safe_get_dimlen)

  END FUNCTION safe_get_dimlen


!------------------------------------------------------------------------------
!S+
! NAME:
!       insert_neighbour
!
! PURPOSE:
!>      Insert a neighbour (squared-distance, source, x_index, y_index) into a
!>      neighbourhood that is already sorted by increasing distance. The
!>      neighbourhood may be modified by adding a neighbour and potentially
!>      shuffling/removing existing ones the assigned may be increased by 1 if
!>      the number of neighbours stored in the neighbourhood increases
!
! CATEGORY:
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL insert_neighbour(neighbourhood,distance,source,x_index,y_index)
!
! ARGUMENTS:
!>@ARG{neighbourhood, inout, NEIGHBOURHOOD_ENTRY} the neighbourhood into which the neigbour is to be inserted
!>@ARG{distance, in, REAL} the squared distance of the neighbour
!>@ARG{source, in, INTEGER} the source of the neighbour (orphan|main, a|b)
!>@ARG{x_index, in, INTEGER} the x-index (column) of the pixel on the IR grid
!>@ARG{y_index, in, INTEGER} the y-index (row) of the pixel on the IR grid
!
! CALLS:
!       None
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     29/01/21  NM  Creation
!S-
!------------------------------------------------------------------------------
  PURE SUBROUTINE insert_neighbour(neighbourhood,distance,source,x_index,y_index)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    TYPE(NEIGHBOURHOOD_ENTRY), INTENT(inout) :: neighbourhood
    REAL,                      INTENT(in)    :: distance
    INTEGER,                   INTENT(in)    :: source
    INTEGER,                   INTENT(in)    :: x_index
    INTEGER,                   INTENT(in)    :: y_index

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: insert_after_index

    ! work out the position to insert the neighbour.  zero means insert at the first position.
    insert_after_index = 0
    IF (neighbourhood%n > 0) THEN
      ! Some neighbours already assigned, so work out where to insert the neighbour (if anywhere)
      DO insert_after_index = neighbourhood%n, 0, -1
        IF (insert_after_index == 0) EXIT
        IF (distance >= neighbourhood%d(insert_after_index)) EXIT
        IF (insert_after_index < MAX_K_NEAREST_NEIGHBOURS) THEN
          ! Shuffle neighbours up to make space for the neighbour to be added
          neighbourhood%d(insert_after_index+1) = neighbourhood%d(insert_after_index)
          neighbourhood%source(insert_after_index+1) = neighbourhood%source(insert_after_index)
          neighbourhood%x(insert_after_index+1) = neighbourhood%x(insert_after_index)
          neighbourhood%y(insert_after_index+1) = neighbourhood%y(insert_after_index)
        END IF
      END DO
    END IF

    ! insert the neighbour if it belongs in the neighbourhood
    IF (insert_after_index < MAX_K_NEAREST_NEIGHBOURS) THEN
      neighbourhood%d(insert_after_index+1) = distance
      neighbourhood%source(insert_after_index+1) = source
      neighbourhood%x(insert_after_index+1) = x_index
      neighbourhood%y(insert_after_index+1) = y_index
      IF (neighbourhood%n < MAX_K_NEAREST_NEIGHBOURS) THEN
        neighbourhood%n = neighbourhood%n + 1
      END IF
    END IF

  END SUBROUTINE insert_neighbour


!------------------------------------------------------------------------------
!S+
! NAME:
!       find_neighbours
!
! PURPOSE:
!>      Find the closest radiance band pixels to a specific IR band pixel
!
! CATEGORY:
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL find_neighbours(ir_x_index, ir_y_index, ir_x, ir_y, &
!                            main_cartesian, orphan_cartesian, neighborhood)
!
! ARGUMENTS:
!>@ARG{ir_x_index, in, INTEGER} the column of the IR pixel
!>@ARG{ir_y_index, in, INTEGER} the row of the IR pixel
!>@ARG{ir_x, in, REAL} across track distance of the IR pixel
!>@ARG{ir_y, in, REAL} along track distance of the IR pixel
!>@ARG{main_cartesian, in, Locations} visible main pixel locations and confidences
!>@ARG{orphan_cartesian, in, Locations} viisble orphan pixel locations
!>@ARG{stripe, in, CHARACTER(1)} the stripe which the cartesian data represents, 'a' or 'b'
!>@ARG{neighbourhood, inout, NEIGHBOURHOOD_ENTRY} the neighbourhood map that is being constructed
!
! CALLS:
!       insert_neighbour
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     23/10/20  NM  Creation
!     06/11/20  NM  Refactor after code review with CB,AW,OE
!     13/05/21  OE  Rewrite to avoid array temporaries, improve performace
!S-
!------------------------------------------------------------------------------
  PURE SUBROUTINE find_neighbours(ir_x_index, ir_y_index, ir_x, ir_y, &
                                  main_cartesian, orphan_cartesian, stripe, neighbourhood)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    INTEGER,                   INTENT(in)  :: ir_x_index
    INTEGER,                   INTENT(in)  :: ir_y_index
    REAL,                      INTENT(in)  :: ir_x
    REAL,                      INTENT(in)  :: ir_y
    TYPE(Locations),           INTENT(in)  :: main_cartesian
    TYPE(Locations),           INTENT(in)  :: orphan_cartesian
    CHARACTER(1),              INTENT(in)  :: stripe
    TYPE(NEIGHBOURHOOD_ENTRY), INTENT(out) :: neighbourhood

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: vis_x_min_index, vis_x_max_index, vis_y_min_index, vis_y_max_index
    REAL :: max_dist_sq
    INTEGER :: visible_width, visible_height
    INTEGER :: main_source, orphan_source

    INTEGER :: elem, line
    REAL :: v

    IF (stripe == 'a') THEN
      main_source = MAIN_PIXEL_SOURCE_A
      orphan_source = ORPHAN_PIXEL_SOURCE_A
    ELSE
      main_source = MAIN_PIXEL_SOURCE_B
      orphan_source = ORPHAN_PIXEL_SOURCE_B
    END IF

    visible_width = main_cartesian%width
    visible_height = main_cartesian%height

    ! if the across or along track distances of the IR pixel are missing, set all neighbours to missing and
    ! return immediately
    IF (ir_x == MISSING_R .or. ir_y == MISSING_R) THEN
      RETURN
    ENDIF

    ! First work out the search window in terms of radiance indices
    vis_x_min_index = MAX(2*ir_x_index - search_width, 1)
    vis_x_max_index = MIN(2*ir_x_index + search_width, visible_width)
    vis_y_min_index = MAX(2*ir_y_index - search_height,1)
    vis_y_max_index = MIN(2*ir_y_index + search_height, visible_height)

    ! Compute the threshold for squared distances - neighbours must be closer than this
    max_dist_sq = MAX_NEIGHBOUR_DISTANCE**2


    ! iterate over the search window and build an ordered list of the closest K
    ! values in the neighborhood structure
    DO line = vis_y_min_index, vis_y_max_index
      ! main pixel array
      DO elem = vis_x_min_index, vis_x_max_index
        ! ignore missing and cosmetically filled pixels
        IF (main_cartesian%x(elem, line) == MISSING_R) CYCLE
        IF (IAND(main_cartesian%conf(elem, line), COSMETIC_PIXEL_MASK) /= 0) CYCLE

        v = (main_cartesian%x(elem, line) - ir_x) ** 2 + (main_cartesian%y(elem, line) - ir_y) ** 2
        IF (v >= max_dist_sq) CYCLE   ! ignore if the squared distance is too far away
        CALL insert_neighbour(neighbourhood, v, main_source, elem, line)
      END DO

      ! orphan pixel array
      DO elem = 1, orphan_cartesian%width
        IF (orphan_cartesian%x(elem, line) == MISSING_R) EXIT

        v = (orphan_cartesian%x(elem, line) - ir_x) ** 2 + (orphan_cartesian%y(elem, line) - ir_y) ** 2
        IF (v >= max_dist_sq) CYCLE
        CALL insert_neighbour(neighbourhood, v, orphan_source, elem, line)
      END DO
    END DO

  END SUBROUTINE find_neighbours


!------------------------------------------------------------------------------
!S+
! NAME:
!       build_neighbourhood_map
!
! PURPOSE:
!>      Find the K closest radiance band pixels to a every IR band pixel and store
!>      the resulting mapping in a neighbourhood map.  Considers main and orphan
!>      band pixels.
!
! CATEGORY:
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL build_neighbourhood_map(ir, visnir, orphan, neighbourhood)
!
! ARGUMENTS:
!>@ARG{ir, in, Locations} ir main pixel locations
!>@ARG{visnir, in, Locations} visible main pixel locations and confidences
!>@ARG{orphan, in, Locations} visible orphan pixel locations
!>@ARG{neighbourhood, inout, NEIGHBOURHOOD_ENTRY} the neighbourhood map that is being constructed
!
! CALLS:
!       find_neighbours
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     23/10/20  NM  Creation
!     06/11/20  NM  Refactor after code review with CB,AW,OE
!S-
!------------------------------------------------------------------------------
  PURE SUBROUTINE build_neighbourhood_map(ir, visnir, orphan, &
                                          neighbourhood)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    TYPE(Locations),         INTENT(in)    :: ir
    TYPE(Locations),         INTENT(in)    :: visnir
    TYPE(Locations),         INTENT(in)    :: orphan
    TYPE(NEIGHBOURHOOD_MAP), INTENT(inout) :: neighbourhood

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: elem, line
    REAL :: x, y
    CHARACTER(1) :: stripe

    IF (neighbourhood%include_a_stripe) THEN
      stripe = 'a'
    ELSE
      stripe = 'b'
    END IF

    ! Now loop over each IR pixel location and find its neighbours
    DO line = 1,ir%height
      DO elem = 1,ir%width
        x = ir%x(elem,line)
        y = ir%y(elem,line)
        CALL find_neighbours(elem, line, x, y, visnir, orphan, &
                             stripe, neighbourhood%entries(elem, line))
      END DO
    END DO
  END SUBROUTINE build_neighbourhood_map


!------------------------------------------------------------------------------
!F+
! NAME:
!       apply_function
!
! PURPOSE:
!>      Aggregate a neighbourhood of radiance pixel values using the specified function.
!>      Customise according to taste to add new functions.
!
! CATEGORY:
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       aggregated_value = apply_function(function_code, pixel_values, pixel_value_count)
!
! ARGUMENTS:
!>@ARG{function_code, in, INTEGER} the code of the function to apply to neighbourhood, see FUNCTION_* parameters
!>@ARG{values, in, REAL\, DIMENSION(:)} array containing values to aggregate
!>@ARG{value_count, in, INTEGER} the number of values in the values array
!
! FUNCTION RESULT:
!>@RES{aggregated_value, REAL} result of applying the function to the array of pixel values
!
! CALLS:
!       None
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     17/02/21  NM  Creation
!F-
!------------------------------------------------------------------------------
  PURE FUNCTION apply_function(function_code, values, value_count)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    INTEGER,            INTENT(in) :: function_code
    REAL, DIMENSION(:), INTENT(in) :: values
    INTEGER,            INTENT(in) :: value_count

    ! ----------------
    ! Function Result
    ! ----------------
    REAL :: apply_function

    ! ---------------
    ! Local Variables
    ! ---------------
    REAL :: vsum, vsumsq

    apply_function = MISSING_R
    IF (value_count <= 0) RETURN

    vsum = 0.0
    vsumsq = 0.0

    SELECT CASE (function_code)
      CASE (FUNCTION_MEAN)
        vsum = SUM(values(1:value_count))
        apply_function = vsum/value_count
      CASE (FUNCTION_SD)
        IF (value_count > 1) THEN
          vsum = SUM(values(1:value_count))
          vsumsq = SUM(values(1:value_count)**2)
          apply_function = SQRT((vsumsq/value_count) - (vsum/value_count)**2)
        ELSE
          apply_function = 0
        END IF
      CASE (FUNCTION_MAX)
        apply_function = MAXVAL(values(1:value_count))
      CASE (FUNCTION_MIN_MAX_DIFF)
        apply_function = MAXVAL(values(1:value_count)) - MINVAL(values(1:value_count))
    END SELECT

  END FUNCTION apply_function



!------------------------------------------------------------------------------
!S+
! NAME:
!       apply_neighbours
!
! PURPOSE:
!>      Aggregate radiance pixels to regrid onto the IR band.  Currently, this
!>      computes the mean value of the neighbours.  Customise according to taste.
!
! CATEGORY:
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL apply_neighbours(vis_radiance_a, vis_orphans_a, vis_exception_a, &
!                             vis_radiance_b, vis_orphans_b, vis_exception_b, &
!                             neighbourhood, vis_output_radiance, function_code)
!
! ARGUMENTS:
!>@ARG{vis_radiance_a, in, REAL\, DIMENSION(:\,:)} matrix of (col,row) radiance band pixels, stripe a
!>@ARG{vis_orphans_a, in, REAL\, DIMENSION(:\,:)} matrix of (col,row) radiance orphan band pixels, stripe a
!>@ARG{vis_exception_a, in, INTEGER\, DIMENSION(:\,:)} matrix of (col,row) exception values for radiance band pixels, stripe a
!>@ARG{vis_radiance_b, in, REAL\, DIMENSION(:\,:)} matrix of (col,row) radiance band pixels, stripe b
!>@ARG{vis_orphans_b, in, REAL\, DIMENSION(:\,:)} matrix of (col,row) radiance orphan band pixels, stripe b
!>@ARG{vis_exception_b, in, INTEGER\, DIMENSION(:\,:)} matrix of (col,row) exception values for radiance band pixels, stripe b
!>@ARG{neighbourhood, in, NEIGHBOURHOOD_ENTRY\, DIMENSION(:\,:)} the neighbourhood map that was populated by build_neighbourhood_map
!>@ARG{effective_k, in, INTEGER} use the closest n neighbours, must be less than or equal to MAX_K_NEAREST_NEIGHBOURS
!>@ARG{vis_output_radiance, inout, REAL\, DIMENSION(:\,:)} matrix of (col,row) to store output aggregated/regridded pixels
!>@ARG{function_code, in, INTEGER} the name of the function to apply to neighbourhood, see FUNCTION_* parameters
!
! CALLS:
!       apply_function
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     23/10/20  NM  Creation
!     06/11/20  NM  Refactor after code review with CB,AW,OE
!     17/02/21  NM  Refactor out apply_function
!S-
!------------------------------------------------------------------------------
  PURE SUBROUTINE apply_neighbours(vis_radiance_a, vis_orphans_a, vis_exception_a, &
                                   vis_radiance_b, vis_orphans_b, vis_exception_b, &
                                   neighbourhood, effective_k, vis_output_radiance, function_code)

    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------

    REAL,    DIMENSION(:,:), INTENT(in)    :: vis_radiance_a
    REAL,    DIMENSION(:,:), INTENT(in)    :: vis_orphans_a
    INTEGER, DIMENSION(:,:), INTENT(in)    :: vis_exception_a
    REAL,    DIMENSION(:,:), INTENT(in)    :: vis_radiance_b
    REAL,    DIMENSION(:,:), INTENT(in)    :: vis_orphans_b
    INTEGER, DIMENSION(:,:), INTENT(in)    :: vis_exception_b
    TYPE(NEIGHBOURHOOD_MAP), INTENT(in)    :: neighbourhood
    INTEGER,                 INTENT(IN)    :: effective_k
    REAL,    DIMENSION(:,:), INTENT(inout) :: vis_output_radiance
    INTEGER,                 INTENT(in)    :: function_code

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: ir_x, ir_y, ik, n_x, n_y, n_s
    REAL :: lookup
    INTEGER :: ir_width, ir_height

    REAL, DIMENSION(MAX_K_NEAREST_NEIGHBOURS) :: values
    INTEGER :: vcount

    ir_width = SIZE(vis_output_radiance,1)
    ir_height = SIZE(vis_output_radiance,2)

    ! Loop over all IR pixels
    DO ir_y = 1, ir_height
      DO ir_x = 1, ir_width
        ! Track the count and sum of the neighbours of the IR pixel
        vcount = 0
        ! Loop over the neighbour indices
        DO ik = 1, effective_k
          n_x = neighbourhood%entries(ir_x,ir_y)%x(ik)
          n_y = neighbourhood%entries(ir_x,ir_y)%y(ik)
          n_s = neighbourhood%entries(ir_x,ir_y)%source(ik)
          lookup = MISSING_R
          ! We check the main pixel exception flags and reject pixels which have any bits
          ! other than the saturation bit set. There is no need to check the orphan pixel
          ! flags as missing / invalid pixels will not be present.
          SELECT CASE(n_s)
            CASE(MAIN_PIXEL_SOURCE_A)
              IF (IOR(vis_exception_a(n_x,n_y),EXCEPTION_SATURATION) == EXCEPTION_SATURATION) THEN
                lookup = vis_radiance_a(n_x,n_y)
              END IF
            CASE(ORPHAN_PIXEL_SOURCE_A)
              lookup = vis_orphans_a(n_x,n_y)
            CASE(MAIN_PIXEL_SOURCE_B)
              IF (IOR(vis_exception_b(n_x,n_y),EXCEPTION_SATURATION) == EXCEPTION_SATURATION) THEN
                lookup = vis_radiance_b(n_x,n_y)
              END IF
            CASE(ORPHAN_PIXEL_SOURCE_B)
              lookup = vis_orphans_b(n_x,n_y)
          END SELECT

          IF (lookup /= MISSING_R) THEN
            vcount = vcount + 1
            values(vcount) = lookup
          END IF
        END DO
        vis_output_radiance(ir_x,ir_y) = apply_function(function_code,values,vcount)
      END DO
    END DO
  END SUBROUTINE apply_neighbours


!------------------------------------------------------------------------------
!S+
! NAME:
!       apply_simple_aggregation
!
! PURPOSE:
!>      Aggregate radiance pixels to regrid onto the IR band.  Currently, this
!>      computes the mean value of the visible pixels at locations:
!>      (2x,2y),(2x+1,2y),(2x+1,2y+1),(2x,2y+1)
!
! CATEGORY:
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL apply_simple_aggregation(vis_radiance, vis_output_radiance, function_code)
!
! ARGUMENTS:
!>@ARG{vis_radiance, in, REAL\, DIMENSION(:\,:)} matrix of (col,row) radiance band pixels
!>@ARG{vis_output_radiance, inout, REAL\, DIMENSION(:\,:)} matrix of (col,row) to store output aggregated/regridded pixels
!>@ARG{function_code, in, INTEGER} the name of the function to apply, see FUNCTION_* parameters
!
! CALLS:
!       apply_function
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     11/01/21  NM  Creation
!     17/02/21  NM  Refactor out apply_function
!     28/06/21  OE  Fix for granules with odd number of scanlines
!S-
!------------------------------------------------------------------------------
  PURE SUBROUTINE apply_simple_aggregation(vis_radiance, vis_output_radiance, function_code)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    REAL, DIMENSION(:,:), INTENT(in)    :: vis_radiance
    REAL, DIMENSION(:,:), INTENT(inout) :: vis_output_radiance
    INTEGER,              INTENT(in)    :: function_code

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: ir_height, ir_width, vis_height
    INTEGER :: ir_x, ir_y, x_off, y_off, vcount
    REAL :: vis_value
    REAL, DIMENSION(4) :: values

    ir_width = SIZE(vis_output_radiance,1)
    ir_height = SIZE(vis_output_radiance,2)
    vis_height = SIZE(vis_radiance, 2)

    ! -- Some granules do not have exactly 2x as many visible scanlines as IR
    !    so check heights of both bands
    DO ir_y=1, MIN(ir_height, vis_height/2)
      DO ir_x=1, ir_width
        vcount = 0
        DO y_off = -1, 0
          DO x_off = -1, 0
            vis_value = vis_radiance(2*ir_x+x_off, 2*ir_y+y_off)
            IF (vis_value /= MISSING_R) THEN
              vcount = vcount + 1
              values(vcount) = vis_value
            END IF
          END DO
        END DO
        vis_output_radiance(ir_x,ir_y) = apply_function(function_code,values,vcount)
      END DO
    END DO

    ! -- Process the final row if we have an odd number of visible scanlines
    IF (ir_y <= ir_height .AND. MODULO(vis_height, 2) == 1) THEN
      DO ir_x=1, ir_width
        vcount = 0
        DO x_off = -1, 0
          vis_value = vis_radiance(2*ir_x+x_off, 2*ir_y-1)
          IF (vis_value /= MISSING_R) THEN
            vcount = vcount + 1
            values(vcount) = vis_value
          END IF
        END DO
        vis_output_radiance(ir_x,ir_y) = apply_function(function_code,values,vcount)
      END DO
    END IF

  END SUBROUTINE apply_simple_aggregation


!------------------------------------------------------------------------------
!S+
! NAME:
!       load_radiance_data
!
! PURPOSE:
!>      Loads data from scene file into radiance, orphan and exception arrays
!
! CATEGORY:
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL load_radiance_data(path, view, band, stripe, &
!                               vis_radiance, vis_orphans, vis_exception)
!
! ARGUMENTS:
!>@ARG{path, in, CHARACTER(LEN=*)} path to a folder storing the various files from an SLSTR scene
!>@ARG{view, in, CHARACTER(LEN=*)} 'o' for oblique view or 'n' for nadir view
!>@ARG{band, in, INTEGER} the radiance band to process, in the range 1 to 6
!>@ARG{stripe, in, CHARACTER(1)} the stripe which the cartesian data represents, 'a' or 'b'
!>@ARG{vis_radiance, out, REAL\, DIMENSION(:\,:)\, ALLOCATABLE} (output) 2D array to hold the output visible pixels
!>@ARG{vis_orphans, out, REAL\, DIMENSION(:\,:)\, ALLOCATABLE} (output) 2D array to hold the output orphan pixels
!>@ARG{vis_exception, out, REAL\, DIMENSION(:\,:)\, ALLOCATABLE} (output) 2D array to hold the output exception values
!
! CALLS:
!       safe_open
!       safe_get_real_data
!       safe_get_int_data
!       nf90_close
!       Path_Join
!
! SIDE EFFECTS:
!       File i/o
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     4/02/21  NM  Creation
!S-
!------------------------------------------------------------------------------
  SUBROUTINE load_radiance_data(path, view, band, stripe, &
                                vis_radiance, vis_orphans, vis_exception)
    USE GbcsPath
    USE netcdf
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    CHARACTER(LEN=*),                     INTENT(in) :: path
    CHARACTER(LEN=*),                     INTENT(in) :: view
    INTEGER,                              INTENT(in) :: band
    CHARACTER(LEN=1),                     INTENT(in) :: stripe
    REAL,    DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: vis_radiance
    REAL,    DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: vis_orphans
    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: vis_exception

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: ncid, dimid, status
    INTEGER :: width, height, orphans
    CHARACTER(40) :: orphan_name, radiance_name, exception_name


    WRITE(radiance_name, '("S",I1,"_radiance_",A1,A1)') band, stripe, view
    WRITE(orphan_name, '("S",I1,"_radiance_orphan_",A1,A1)') band, stripe, view
    WRITE(exception_name, '("S",I1,"_exception_",A1,A1)') band, stripe, view

    ncid = safe_open(Path_Join(path, TRIM(radiance_name)//'.nc'))

    ! Work out the sizes of the arrays
    width = safe_get_dimlen(ncid, 'columns')
    height = safe_get_dimlen(ncid, 'rows')
    orphans = safe_get_dimlen(ncid, 'orphan_pixels')

    ALLOCATE(vis_radiance(width, height), &
             vis_orphans(orphans, height), &
             vis_exception(width, height))

    CALL safe_get_real_data(ncid, TRIM(radiance_name), vis_radiance, MISSING_R)
    CALL safe_get_real_data(ncid, TRIM(orphan_name), vis_orphans, MISSING_R)
    CALL safe_get_int_data(ncid, TRIM(exception_name), vis_exception, MISSING_I)
    status = nf90_close(ncid)

  END SUBROUTINE load_radiance_data


!------------------------------------------------------------------------------
!S+
! NAME:
!       process_scene_band
!
! PURPOSE:
!>      Process visible pixel values from a radiance band to regrid onto an IR band
!>      and write the resulting dataset to an allocated output array.
!
! CATEGORY:
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL process_scene_band(path, view, band, vis_ir_radiance, &
!                               neighbourhood, effective_k, function_code)
!
!
!
! ARGUMENTS:
!>@ARG{path, in, CHARACTER(LEN=*)} path to a folder storing the various files from an SLSTR scene
!>@ARG{view, in, CHARACTER(LEN=*)} 'o' for oblique view or 'n' for nadir view
!>@ARG{band, in, INTEGER} the radiance band to process, in the range 1 to 6
!>@ARG{vis_ir_radiance, inout, REAL\, DIMENSION(:\,:)} (output) 2D array to hold the output visible pixels arranged
!                (COL,ROW) on the IR grid.  This will be exactly half the width and height of the vis_radiance array
!>@ARG{neighborhood, in, NEIGHBOURHOOD_MAP} the neighbourhood map that has been constructed
!>@ARG{effective_k, in, INTEGER} use the closest n neighbours, must be less than or equal to MAX_K_NEAREST_NEIGHBOURS
!>@ARG{function_code, in, INTEGER} the name of the function to apply, see FUNCTION_* parameters
!
! CALLS:
!       load_radiance_data
!       apply_neighbours (if neighbourhood allocated)
!       apply_simple_aggregation
!
! SIDE EFFECTS:
!     All cells in vis_ir_radiance are assigned
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! EXAMPLE:
!> @EX
!>      USE SLSTR_Preprocessor
!>      ! Declare variable to hold the neighbourhood
!>      TYPE(NEIGHBOURHOOD_MAP) :: neighbourhood
!>      ! Build neighbourhood mapping for nadir, opening files as needed from scene folder
!>      ! This should eventually work with a folder or path to a compressed file, right now has to be a folder
!>      CALL compute_scene_neighbourhood('n','/path/to/scene',neighbourhood,'a')
!>      ! Process a visual scene from radiance band 4
!>      CALL process_scene_band('n','/path/to/scene',4,vis_output_radiance,neighbourhood,5,FUNCTION_MEAN)
!> @endcode
!
! CREATION HISTORY:
!     11/11/20  NM  Creation
!     29/01/21  NM  Add stripe support
!
!S-
!------------------------------------------------------------------------------
  SUBROUTINE process_scene_band(path, view, band, vis_ir_radiance, &
                                neighbourhood, effective_k, function_code)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    CHARACTER(LEN=*),        INTENT(in)    :: path
    CHARACTER(LEN=*),        INTENT(in)    :: view
    INTEGER,                 INTENT(in)    :: band
    REAL,    DIMENSION(:,:), INTENT(inout) :: vis_ir_radiance
    TYPE(NEIGHBOURHOOD_MAP), INTENT(in)    :: neighbourhood
    INTEGER,                 INTENT(IN)    :: effective_k
    INTEGER,                 INTENT(in)    :: function_code

    ! ---------------
    ! Local Variables
    ! ---------------
    REAL, ALLOCATABLE, DIMENSION(:,:) :: vis_radiance_a, vis_radiance_b
    REAL, ALLOCATABLE, DIMENSION(:,:) :: vis_orphans_a, vis_orphans_b
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: vis_exception_a, vis_exception_b

    IF (function_code < FUNCTION_MEAN .or. function_code > FUNCTION_MIN_MAX_DIFF) THEN
        PRINT *, 'Unrecognized function code', function_code
        RETURN
    END IF

    IF (effective_k <= 0 .or. effective_k > MAX_K_NEAREST_NEIGHBOURS) THEN
        PRINT *, 'effective_k must be greater than 0 and less than or equal to MAX_K_NEAREST_NEIGHBOURS', effective_k
        RETURN
    END IF

    IF (ALLOCATED(neighbourhood%entries)) THEN
      ! if the neighbourhood contains a-stripe neighbours, then load radiance data from that stripe
      IF (neighbourhood%include_a_stripe) THEN
        CALL load_radiance_data(path,view,band,'a',vis_radiance_a,vis_orphans_a,vis_exception_a)
      END IF

      ! if the neighbourhood contains b-stripe neighbours, then load radiance data from that stripe
      IF (neighbourhood%include_b_stripe) THEN
        CALL load_radiance_data(path,view,band,'b',vis_radiance_b,vis_orphans_b,vis_exception_b)
      END IF

      CALL apply_neighbours(vis_radiance_a,vis_orphans_a,vis_exception_a, &
              vis_radiance_b,vis_orphans_b,vis_exception_b, &
              neighbourhood,effective_k,vis_ir_radiance,function_code)
    ELSE
      CALL load_radiance_data(path,view,band,'a',vis_radiance_a,vis_orphans_a,vis_exception_a)
      CALL apply_simple_aggregation(vis_radiance_a,vis_ir_radiance,function_code)
    END IF

  END SUBROUTINE process_scene_band


!------------------------------------------------------------------------------
!S+
! NAME:
!       align_cosmetics
!
! PURPOSE:
!>      Try to adjust the cartesian coordinates of cosmetic pixels so that they
!>      exactly match the coordinates of the pixels which contributed their values
!
! CATEGORY:
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL align_cosmetics(ir_cartesian, s8_bt)
!
! ARGUMENTS:
!>@ARG{ir_cartesian, inout, Locations}  ir main pixel locations
!>@ARG{s8_bt, in, REAL\, DIMENSION(:\,:)}  values from he S8_BT band
!
! CALLS:
!       None
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     25/11/20  NM  Creation
!S-
!------------------------------------------------------------------------------
  PURE SUBROUTINE align_cosmetics(ir_cartesian, s8_bt)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    TYPE(Locations),         INTENT(inout) :: ir_cartesian
    REAL,    DIMENSION(:,:), INTENT(in)    :: s8_bt

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: ir_x, ir_y, x, y, off_x, off_y, n_x, n_y
    INTEGER :: ir_width, ir_height
    INTEGER :: miss_count, value_match_count
    REAL :: cosmetic_value, xc, yc, xc2, yc2, dist, n_dist

    ir_width = ir_cartesian%width
    ir_height = ir_cartesian%height

    miss_count = 0   ! count the number of cosmetic pixels that could not be aligned

    ! Now loop over each IR pixel location
    DO ir_y = 1,ir_height
      DO ir_x = 1,ir_width
        ! check to see if this pixel is cosmetically filled
        IF (IAND(ir_cartesian%conf(ir_x,ir_y),COSMETIC_PIXEL_MASK) == COSMETIC_PIXEL_MASK) THEN
          ! search the 8 immediate neighbours for pixels with the same value in band 8
          ! and find the closest one
          cosmetic_value = s8_bt(ir_x,ir_y)
          value_match_count = 0
          xc = ir_cartesian%x(ir_x,ir_y)
          yc = ir_cartesian%y(ir_x,ir_y)
          n_dist = 1000000000.0
          n_x = -1
          n_y = -1
          DO off_y = -1,1
            DO off_x = -1,1
              IF (off_y /= 0 .or. off_x /= 0) THEN
                ! calculate the x and y indices of the neighbour
                x = ir_x + off_x
                y = ir_y + off_y
                ! check the indices are valid and the neighbour is not itself cosmetically filled
                IF (x >= 1 .and. x <= ir_width .and. y >= 1 .and. y <= ir_height .and. &
                    IAND(ir_cartesian%conf(x,y),COSMETIC_PIXEL_MASK) == COSMETIC_PIXEL_MASK) THEN
                  IF (cosmetic_value == s8_bt(x,y)) THEN
                    ! neighbour has the same value in band 8, see if it the closest one found?
                    value_match_count = value_match_count + 1
                    xc2 = ir_cartesian%x(x,y)
                    yc2 = ir_cartesian%y(x,y)
                    dist = SQRT((xc-xc2)**2 + (yc-yc2)**2)
                    IF (dist < n_dist) THEN
                      n_x = x
                      n_y = y
                      n_dist = dist
                    END IF
                  END IF
                END IF
              END IF
            END DO
          END DO
          IF (value_match_count == 0) THEN
            ! we did not find the neighbouring pixel with the same value
            miss_count = miss_count + 1
          ELSE
            ! align the x and y coordinates of the cosmetically filled pixel to those of the
            ! neighbour
            ir_cartesian%x(ir_x,ir_y) = ir_cartesian%x(n_x,n_y)
            ir_cartesian%y(ir_x,ir_y) = ir_cartesian%y(n_x,n_y)
          END IF
        END IF
      END DO
    END DO
  END SUBROUTINE align_cosmetics


!------------------------------------------------------------------------------
!S+
! NAME:
!       merge_neighbourhoods
!
! PURPOSE:
!>      merge 2 neighbourhoods (presumably computed from different stripes)
!>      retaining the nearest K neighbours. The second of the neighbourhoods
!>      will hold the merged neighbourhood
!
! CATEGORY:
!
! LANGUAGE:
!      Fortran-95
!
! CALLING SEQUENCE:
!       CALL merge_neighbourhoods(neighbourhood_1_in,neighbourhood_2_inout)
! ARGUMENTS:
!>@ARG{neighbourhood_1_in, in, NEIGHBOURHOOD_ENTRY\, DIMENSION(:\,:)} first neighbourhood to merge
!>@ARG{neighbourhood_2_inout, inout, NEIGHBOURHOOD_ENTRY\, DIMENSION(:\,:)} second neighbourhood to be merged/modified
!
! CALLS:
!       insert_neighbour
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     29/01/21  NM  Creation
!S-
!------------------------------------------------------------------------------
  PURE SUBROUTINE merge_neighbourhoods(neighbourhood_1_in,neighbourhood_2_inout)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    TYPE(NEIGHBOURHOOD_MAP), INTENT(in)    :: neighbourhood_1_in
    TYPE(NEIGHBOURHOOD_MAP), INTENT(inout) :: neighbourhood_2_inout

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: width, height, x, y, i, source

    width = SIZE(neighbourhood_1_in%entries,1)
    height = SIZE(neighbourhood_1_in%entries,2)

    neighbourhood_2_inout%include_a_stripe = neighbourhood_1_in%include_a_stripe .or. neighbourhood_2_inout%include_a_stripe
    neighbourhood_2_inout%include_b_stripe = neighbourhood_1_in%include_b_stripe .or. neighbourhood_2_inout%include_b_stripe

    ! Merge sort would be more efficient here as the input neighbourhoods are already sorted
    ! However this would only save a second or so at most, so about 1-2% of the overall execution time
    DO x = 1,width
      DO y = 1, height

        ! now insert cells from the first neighbourhood into the second
        DO i = 1, neighbourhood_1_in%entries(x,y)%n
          IF (neighbourhood_1_in%entries(x,y)%source(i) == INVALID_PIXEL) EXIT

          CALL insert_neighbour(neighbourhood_2_inout%entries(x,y), &
                  neighbourhood_1_in%entries(x,y)%d(i), &
                  neighbourhood_1_in%entries(x,y)%source(i), &
                  neighbourhood_1_in%entries(x,y)%x(i), &
                  neighbourhood_1_in%entries(x,y)%y(i))
        END DO
      END DO
    END DO

  END SUBROUTINE merge_neighbourhoods


!------------------------------------------------------------------------------
!S+
! NAME:
!       Allocate_Locations
!
! PURPOSE:
!>      Allocate a locations structure
!
! CATEGORY:
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL Allocate_Locations(width, height, flags, array)
!
! ARGUMENTS:
!>@ARG{width, in, INTEGER} Width of image
!>@ARG{height, in, INTEGER} Height of image
!>@ARG{flags, in, LOGICAL} If true we will also allocate the conf flags array
!>@ARG{array, out, TYPE\(Locations\)} Pixel locations array
!
! CALLS:
!       None
!
! CREATION HISTORY:
!     27/11/24  OE  Creation
!S-
!------------------------------------------------------------------------------
  SUBROUTINE Allocate_Locations(width, height, flags, array)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    INTEGER, INTENT(in) :: width
    INTEGER, INTENT(in) :: height
    LOGICAL, INTENT(in) :: flags
    TYPE(Locations), INTENT(out) :: array

    array%width = width
    array%height = height
    ALLOCATE(array%x(width, height), array%y(width, height))
    IF (flags) ALLOCATE(array%conf(width, height))

  END SUBROUTINE Allocate_Locations


!------------------------------------------------------------------------------
!S+
! NAME:
!       compute_scene_neighbourhood
!
! PURPOSE:
!>      given a folder containing an SLSTR scene, a view (nadir or oblique) and
!>      stripe (a or b), loads the cartesian locations for visible and ir pixel
!>      grids and computes the neighborhood of visible pixels closest to each IR
!>      pixel.  The neigbourhood can then be used in calls to process_scene_band
!>      to map visible bands from the same scene and view to a new IR based grid.
!>      Neighbourhoods from a and b stripes can also be merged by calling
!>      merge_neighbourhoods
!
! CATEGORY:
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL compute_scene_neighbourhood(path, view, stripe, neighbourhood)
!
! ARGUMENTS:
!>@ARG{path, in, CHARACTER(LEN=*)} path to a folder storing the various files from an SLSTR scene
!>@ARG{view, in, CHARACTER(1)} 'o' for oblique view or 'n' for nadir view
!>@ARG{stripe, in, CHARACTER(1)} pass 'a' to compute the neighbourhood for the a-stripe, or 'b' to use the b-stripe
!>@ARG{neighbourhood, out, NEIGHBOURHOOD_MAP} an array of neighbourhood maps, one for each IR pixel
!
! CALLS:
!       safe_open
!       safe_get_real_data
!       safe_get_int_data
!       align_cosmetics
!       build_neighbourhood_map
!
! SIDE EFFECTS:
!       File i/o
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! EXAMPLE:
!> @EX
!>      USE slstr_scene_preprocessor
!>      ! Declare variable to hold the neighbourhood
!>      TYPE(NEIGHBOURHOOD_ENTRY) :: neighbourhood
!>      ! Build neighbourhood mapping for nadir, opening files as needed from scene folder
!>      ! This should eventually work with a folder or path to a compressed file, right now has to be a folder
!>      CALL compute_scene_neighbourhood('n','/path/to/scene',neighbourhood,'a')
!
! CREATION HISTORY:
!     11/11/20  NM  Creation
!     29/01/21  NM  Add stripe support
!S-
!------------------------------------------------------------------------------
  SUBROUTINE compute_scene_neighbourhood(path, view, stripe, neighbourhood)
    USE netcdf
    USE GbcsPath
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    CHARACTER(LEN=*),        INTENT(in)  :: path
    CHARACTER(LEN=1),        INTENT(in)  :: view
    CHARACTER(LEN=1),        INTENT(in)  :: stripe
    TYPE(NEIGHBOURHOOD_MAP), INTENT(out) :: neighbourhood

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER :: status
    INTEGER :: vis_ncid, flags_ncid, ir_ncid, iflags_ncid, s8_bt_ncid
    INTEGER :: visible_width, visible_height, ir_width, ir_height, orphan_width
    INTEGER :: x, y

    ! Create data structures to hold the x- and y- 2D arrays for the IR and main / orphan pixels
    TYPE(Locations) :: ir_cartesian
    TYPE(Locations) :: main_cartesian
    TYPE(Locations) :: orphan_cartesian
    REAL, ALLOCATABLE, DIMENSION(:,:) :: s8_bt


    IF (stripe == 'a') THEN
      neighbourhood%include_a_stripe = .true.
      neighbourhood%include_b_stripe = .false.
    ELSE
      neighbourhood%include_a_stripe = .false.
      neighbourhood%include_b_stripe = .true.
    END IF

    vis_ncid = safe_open(Path_Join(path, 'cartesian_'//stripe//view//'.nc'))
    ir_ncid = safe_open(Path_Join(path, 'cartesian_i'//view//'.nc'))
    flags_ncid = safe_open(Path_Join(path, 'flags_'//stripe//view//'.nc'))

    ! Work out the sizes of the arrays
    visible_width = safe_get_dimlen(vis_ncid, 'columns')
    visible_height = safe_get_dimlen(vis_ncid, 'rows')
    orphan_width = safe_get_dimlen(vis_ncid, 'orphan_pixels')
    ir_width = safe_get_dimlen(ir_ncid, 'columns')
    ir_height = safe_get_dimlen(ir_ncid, 'rows')

    ! Allocate the arrays for each data structure
    ALLOCATE(neighbourhood%entries(ir_width,ir_height))
    CALL Allocate_Locations(ir_width, ir_height, align_cosmetic_pixels, ir_cartesian)
    CALL Allocate_Locations(visible_width,visible_height, .True., main_cartesian)
    CALL Allocate_Locations(orphan_width,visible_height, .False., orphan_cartesian)

   ! Load the cartesian coordinates arrays and flags from the scene
    CALL safe_get_real_data(vis_ncid, 'x_'//stripe//view, main_cartesian%x, MISSING_R)
    CALL safe_get_real_data(vis_ncid, 'y_'//stripe//view, main_cartesian%y, MISSING_R)
    CALL safe_get_int_data(flags_ncid, 'confidence_'//stripe//view, main_cartesian%conf, MISSING_I)

    CALL safe_get_real_data(ir_ncid, 'x_i'//view, ir_cartesian%x, MISSING_R)
    CALL safe_get_real_data(ir_ncid, 'y_i'//view, ir_cartesian%y, MISSING_R)

    status = nf90_close(flags_ncid)
    status = nf90_close(ir_ncid)

    IF (align_cosmetic_pixels) THEN
      ! tweak cartesian of cosmetically filled pixels to be exactly the same as the original (donor) pixels
      s8_bt_ncid = safe_open(Path_Join(path, 'S8_BT_i'//view//'.nc'))
      iflags_ncid = safe_open(Path_Join(path, 'flags_i'//view//'.nc'))
      ALLOCATE(s8_bt(ir_width,ir_height))
      CALL safe_get_int_data(iflags_ncid, 'confidence_i'//view, ir_cartesian%conf, MISSING_I)
      CALL safe_get_real_data(s8_bt_ncid, 'S8_BT_i'//view, s8_bt, MISSING_R)
      status = nf90_close(iflags_ncid)
      status = nf90_close(s8_bt_ncid)

      CALL align_cosmetics(ir_cartesian, s8_bt)
    END IF

    CALL safe_get_real_data(vis_ncid, 'x_orphan_'//stripe//view, orphan_cartesian%x, MISSING_R)
    CALL safe_get_real_data(vis_ncid, 'y_orphan_'//stripe//view, orphan_cartesian%y, MISSING_R)

    status = nf90_close(vis_ncid)

    ! Much of the effort in searching for neighbours involves searching the orphan pixels
    ! First determine how many columns within the orphan matrix we need to search
    ! All valid values are clumped together with low column indices in orphan data, it seems
    orphan_cartesian%width = 0

    DO y = 1,visible_height
      DO x = 1,orphan_width
        IF (orphan_cartesian%x(x,y) == MISSING_R) THEN
          orphan_cartesian%width = MAX(orphan_cartesian%width, x-1)
          EXIT
        END IF
      END DO
    END DO

    ! Now build the neighbourhood, using the loaded cartesian data
    CALL build_neighbourhood_map(ir_cartesian,main_cartesian,orphan_cartesian,neighbourhood)

  END SUBROUTINE compute_scene_neighbourhood

END MODULE SLSTR_Preprocessor
