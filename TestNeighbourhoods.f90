!--------------------------------------------------------------------------------------------
!P+ Preprocess_SLSTR
! NAME:
!     TestNeighbourhoods
!
! PURPOSE:
!
! Test driver program for SLSTR_Preprocessor
!
! Test the similarity/equality of neighbourhoods calcluated using old (search window) and new algorithms
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
!     Written by:  Niall McCarroll 28/06/21
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

SUBROUTINE dump_entry(entry)
  USE SLSTR_Preprocessor
  IMPLICIT NONE

  TYPE(NEIGHBOURHOOD_ENTRY), INTENT(IN) :: entry
  INTEGER :: k, source
  REAL :: x, y, sqdist

  DO k=1,MAX_K_NEAREST_NEIGHBOURS
      source = entry%source(k)
      x = entry%x(k)
      y = entry%y(k)
      sqdist = entry%squared_distances(k)
      IF (source /= NULL_PIXEL_SOURCE) THEN
        PRINT *, source, x, y, sqdist
      END IF
  END DO
END SUBROUTINE dump_entry

!------------------------------------------------------------------------------
!S+
! NAME:
!     compare_neighbourhood_entries
!
! PURPOSE:
!
! CATEGORY:
!
! LANGUAGE:
!     Fortran-95
!
! CALLING SEQUENCE:
!
! ARGUMENTS:
!     old_neighbourhood   - old algorithm computed array of computed neighbourhoods, organised by (COL,ROW)
!     new_neighbourhood   - new algorithm computed array of computed neighbourhoods, organised by (COL,ROW)
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
!     28/06/21  NM  Creation
!
!S-
!------------------------------------------------------------------------------
SUBROUTINE compare_neighbourhood_entries(old_neighbourhood, new_neighbourhood)
  USE SLSTR_Preprocessor
  IMPLICIT NONE

  ! -----------
  ! Arguments
  ! -----------
  TYPE(NEIGHBOURHOOD_MAP), INTENT(IN) :: old_neighbourhood
  TYPE(NEIGHBOURHOOD_MAP), INTENT(IN) :: new_neighbourhood

  ! ---------------
  ! Local Variables
  ! ---------------
  REAL old_x, old_y, new_x, new_y, old_sqdist, new_sqdist
  INTEGER :: row, col, k, rows, cols
  INTEGER :: old_main_count, old_orphan_count, new_main_count, new_orphan_count
  INTEGER :: old_source, new_source
  INTEGER :: old_entries_count, new_entries_count

  old_main_count = 0
  old_orphan_count = 0
  new_main_count = 0
  new_orphan_count = 0

  cols = SIZE(old_neighbourhood%entries,1)
  rows = SIZE(old_neighbourhood%entries,2)

  ! print some handy statistics on the sources of pixels in the neighbourhood
  DO row=1,rows
    DO col=1,cols
      old_entries_count = 0
      new_entries_count = 0
      DO k=1,MAX_K_NEAREST_NEIGHBOURS
        old_source = old_neighbourhood%entries(col,row)%source(k)
        if (old_source == MAIN_PIXEL_SOURCE_A .or. old_source == MAIN_PIXEL_SOURCE_B) then
            old_main_count = old_main_count + 1
            old_entries_count = old_entries_count + 1
        else if (old_source == ORPHAN_PIXEL_SOURCE_A .or. old_source == ORPHAN_PIXEL_SOURCE_B) then
            old_orphan_count = old_orphan_count + 1
            old_entries_count = old_entries_count + 1
        end if
        old_x = old_neighbourhood%entries(col,row)%x(k)
        old_y = old_neighbourhood%entries(col,row)%y(k)
        old_sqdist = old_neighbourhood%entries(col,row)%squared_distances(k)

        new_source = new_neighbourhood%entries(col,row)%source(k)
        new_x = new_neighbourhood%entries(col,row)%x(k)
        new_y = new_neighbourhood%entries(col,row)%y(k)
        new_sqdist = new_neighbourhood%entries(col,row)%squared_distances(k)
        if (new_source == MAIN_PIXEL_SOURCE_A .or. new_source == MAIN_PIXEL_SOURCE_B) then
            new_main_count = new_main_count + 1
            new_entries_count = new_entries_count + 1
        else if (new_source == ORPHAN_PIXEL_SOURCE_A .or. new_source == ORPHAN_PIXEL_SOURCE_B) then
            new_orphan_count = new_orphan_count + 1
            new_entries_count = new_entries_count + 1
        end if
        if (old_sqdist /= new_sqdist) THEN
          PRINT *, 'Sq dists differ at', row, col
          PRINT *, 'OLD'
          CALL dump_entry(old_neighbourhood%entries(col,row))
          PRINT *, 'NEW'
          CALL dump_entry(new_neighbourhood%entries(col,row))
          STOP
        end if
      END DO
      IF (old_entries_count /= new_entries_count) THEN
        PRINT *, 'Entries differ at', row, col
        PRINT *, 'OLD'
        CALL dump_entry(old_neighbourhood%entries(col,row))
        PRINT *, 'NEW'
        CALL dump_entry(new_neighbourhood%entries(col,row))
        STOP
      END IF
    END DO
  END DO
  PRINT *, ''
  PRINT *, '  MAIN PIXELS USED(old)=', old_main_count
  PRINT *, '  MAIN PIXELS USED(new)=', new_main_count
  PRINT *, 'ORPHAN PIXELS USED(old)=', old_orphan_count
  PRINT *, 'ORPHAN PIXELS USED(new)=', new_orphan_count
  PRINT *, ' TOTAL PIXELS USED(old)=', old_main_count+old_orphan_count
  PRINT *, ' TOTAL PIXELS USED(new)=', new_main_count+new_orphan_count

END SUBROUTINE compare_neighbourhood_entries

!------------------------------------------------------------------------------
!S+
! NAME:
!     compare_neighbourhoods
!
! PURPOSE:
!
! CATEGORY:
!
! LANGUAGE:
!     Fortran-95
!
! CALLING SEQUENCE:
!   compare_neighbourhoods('o','/path/to/scene',10)
!
! ARGUMENTS:
!     view_type       - 'o' for oblique view or 'n' for nadir view
!     scene_folder    - path to a folder storing the various input files from an SLSTR scene
!     effective_k     - the value of k
!
! CALLS:
!     Invokes compute_scene_neighbourhood in module slstr_preprocessor
!
! SIDE EFFECTS:
!
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! CREATION HISTORY:
!     28/06/21  NM  Creation
!
!S-
!------------------------------------------------------------------------------
SUBROUTINE compare_neighbourhoods(view_type, scene_folder, effective_k)
  USE SLSTR_Preprocessor
  USE GbcsPath
  USE netcdf
  IMPLICIT NONE

  ! -----------
  ! Arguments
  ! -----------
  CHARACTER(1), INTENT(IN) :: view_type
  CHARACTER(LEN=*), INTENT(IN) :: scene_folder
  INTEGER, INTENT(IN) :: effective_k

  ! ---------------
  ! Local Variables
  ! ---------------

  TYPE(NEIGHBOURHOOD_MAP) :: old_neighbourhood_a, new_neighbourhood_a

  REAL :: old_start_time, old_end_time, old_elapsed_time
  REAL :: new_start_time, new_end_time, new_elapsed_time

  use_new_neighbourhood_algorithm = .false.
  CALL cpu_time(old_start_time)
  CALL compute_scene_neighbourhood(view_type,scene_folder,old_neighbourhood_a,'a')
  CALL cpu_time(old_end_time)

  old_elapsed_time = old_end_time - old_start_time

  use_new_neighbourhood_algorithm = .true.
  CALL cpu_time(new_start_time)
  CALL compute_scene_neighbourhood(view_type,scene_folder,new_neighbourhood_a,'a')
  CALL cpu_time(new_end_time)

  new_elapsed_time = new_end_time - new_start_time

  PRINT *, 'old neighbourhood build time: ', old_elapsed_time
  PRINT *, 'new neighbourhood build time', new_elapsed_time

  CALL compare_neighbourhood_entries(old_neighbourhood_a,new_neighbourhood_a)
END SUBROUTINE compare_neighbourhoods

PROGRAM TestNeighbourhoods
  USE SLSTR_Preprocessor
  IMPLICIT NONE

  ! ---------------
  ! Local Variables
  ! ---------------
  CHARACTER(256) :: scene_folder

  CHARACTER(256) :: option
  CHARACTER(256) :: option_value

  INTEGER :: effective_k, i, max_distance, window_height, window_width
  LOGICAL :: simple, stats, compare

  simple = .false.
  stats = .false.
  compare = .false.
  effective_k = MAX_K_NEAREST_NEIGHBOURS
  max_distance = 4000
  window_height = 0
  window_width = 0

  CALL GET_COMMAND_ARGUMENT(1,scene_folder)

  i = 2
  DO WHILE (i <= COMMAND_ARGUMENT_COUNT())
    CALL GET_COMMAND_ARGUMENT(i,option)
    i = i + 1
    IF (option == '--effective_k') THEN
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

  CALL compare_neighbourhoods('n',scene_folder,effective_k)

END PROGRAM TestNeighbourhoods
