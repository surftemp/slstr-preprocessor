!------------------------------------------------------------------------------
!M+
! NAME:
!       GbcsNetCDF
!
! PURPOSE:
!>      Routines for interfacing with the NetCDF library and applying
!>      automatic scaling where appropriate
!
! CATEGORY:
!       NetCDF
!
! LANGUAGE:
!       Fortran-2003
!
! CALLING SEQUENCE:
!       USE GbcsNetCDF
!
! PUBLIC DATA:
!       None
!
! MODULES:
!       GbcsKinds
!       netcdf
!
! CONTAINS:
!       Get_DimID
!       Dim_Length
!       Get_VarID
!       Get_VarInfo
!       Get_Var_Dimensions
!       Write_Var
!       Quantize
!
! DERIVED TYPES:
!       VarInfo (private)
!
! NOTES:
!
! CREATION HISTORY:
!       Written by:   Owen Embury 08/03/2011
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
MODULE GbcsNetCDF
  ! ------------
  ! Modules used
  ! ------------
  USE GbcsKinds
  USE netcdf
  USE iso_c_binding

  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Get_DimID, Dim_Length
  PUBLIC :: Find_VarID, Get_VarID, Get_VarInfo, Get_Var_Dimensions
  PUBLIC :: Write_Var
  PUBLIC :: Quantize
  PUBLIC :: nf90_strerror

  ! --------------------
  ! Function overloading
  ! --------------------
  INTERFACE Write_Var
    MODULE PROCEDURE Write_Var1F, Write_Var1I, &
                     Write_Var2D, Write_Var2F, &
                     Write_Var2I, Write_Var2S, Write_Var2B
  END INTERFACE Write_Var

  INTERFACE Pack_Data
    MODULE PROCEDURE Pack_Data_Single, Pack_Data_Double
  END INTERFACE Pack_Data

  INTERFACE Unpack_Data
    MODULE PROCEDURE Unpack_Data_Single, Unpack_Data_Double
  END INTERFACE Unpack_Data

  INTERFACE Quantize
    MODULE PROCEDURE Quantize_Scalar, Quantize_1F, Quantize_2F, &
                     Quantize_2D
  END INTERFACE Quantize


  ! -------------
  ! Derived Types
  ! -------------
  TYPE VarInfo
    LOGICAL :: Is_Packed = .FALSE.
    LOGICAL :: Has_Range = .FALSE.
    LOGICAL :: Has_Fill  = .FALSE.
    LOGICAL :: Use_lsd   = .FALSE.
    REAL    :: offset    = 0.0
    REAL    :: scale     = 1.0
    REAL    :: precision = -1.0
    REAL    :: valid_min =-HUGE(0.0)
    REAL    :: valid_max = HUGE(0.0)
    REAL    :: value_min =-HUGE(0.0)
    REAL    :: value_max = HUGE(0.0)
    REAL    :: fillvalue
  END TYPE VarInfo

  ! -----------------
  ! Module parameters
  ! -----------------
  INTEGER, PRIVATE, PARAMETER :: max_nc_dims = 8

  ! -- Module name for error messages
  CHARACTER( * ), PRIVATE, PARAMETER :: MODULE_NAME = &
    'GbcsMod_DataReaders/GbcsNetCDF.f90'

CONTAINS
!------------------------------------------------------------------------------
!F+
! NAME:
!       Get_DimID
!
! PURPOSE:
!>      Find the NetCDF dimid for the given dimension name
!
! CATEGORY:
!       NetCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       dimid = Get_DimID(ncid, name)
!
! ARGUMENTS:
!>@ARG{ncid, in, INTEGER} netCDF file id
!>@ARG{name, in, CHARACTER(*)} dimension name
!
! FUNCTION RESULT:
!>@RES{dimid, INTEGER} netCDF dimension ID or -1 on error
!
! CALLS:
!       nf90_inq_dimid
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!       None
!
! PROCEDURE:
!>@PRO  This is just a wrapper around nf90_inq_dimid which sets the returned
!>      dimid to -1 on error.
!
! CREATION HISTORY:
!       11/03/11  OE  Creation
!F-
!------------------------------------------------------------------------------
  FUNCTION Get_DimID(ncid, name)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    INTEGER,          INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: name

    ! ---------
    ! Function result
    ! ---------
    INTEGER                      :: Get_DimID

    ! ----------------
    ! Local parameters
    ! ----------------
    CHARACTER( * ),  PARAMETER :: ROUTINE_NAME = 'Get_DimID'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: status

    status = nf90_inq_dimid(ncid, name, Get_DimID)
    IF (status /= nf90_noerr) Get_DimID = -1

  END FUNCTION Get_DimID


!------------------------------------------------------------------------------
!F+
! NAME:
!       Dim_Length
!
! PURPOSE:
!>      Get the length of the specified dimension
!
! CATEGORY:
!       NetCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       len = Dim_Length(ncid, name)
!
! ARGUMENTS:
!>@ARG{ncid, in, INTEGER} netCDF file id
!>@ARG{name, in, CHARACTER(*)} dimension name
!
! FUNCTION RESULT:
!>@RES{length, INTEGER} Length of dimension or -1 on error
!
! CALLS:
!       nf90_inq_dimid
!       nf90_inquire_dimension
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!       None
!
! PROCEDURE:
!>@PRO  This is just a wrapper around nf90_inq_dimid / nf90_inquire_dimension
!>      which sets the returned length to -1 on error
!
! CREATION HISTORY:
!       10/03/11  OE  Creation
!F-
!------------------------------------------------------------------------------
  FUNCTION Dim_Length(ncid, name)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    INTEGER,          INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: name

    ! ---------
    ! Function result
    ! ---------
    INTEGER                      :: Dim_Length

    ! ----------------
    ! Local parameters
    ! ----------------
    CHARACTER( * ),  PARAMETER :: ROUTINE_NAME = 'Dim_Length'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: status, dimid

    Dim_Length = -1

    status = nf90_inq_dimid(ncid, name, dimid)
    IF (status /= nf90_noerr) RETURN

    status = nf90_inquire_dimension(ncid, dimid, LEN=Dim_Length)

  END FUNCTION Dim_Length


!------------------------------------------------------------------------------
!F+
! NAME:
!       Get_VarID
!
! PURPOSE:
!>      Find the NetCDF varid for the given variable name
!
! CATEGORY:
!       NetCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       varid = Get_VarID(ncid, name)
!
! ARGUMENTS:
!>@ARG{ncid, in, INTEGER} netCDF file id
!>@ARG{name, in, CHARACTER(*)} dimension name
!>@ARG{errmsg, in, CHARACTER(*)\, OPTIONAL} Display this message on error
!
! FUNCTION RESULT:
!>@RES{varid, INTEGER} netCDF variable ID or -1 on error
!
! CALLS:
!       nf90_inq_varid
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!       None
!
! PROCEDURE:
!>@PRO  This is just a wrapper around nf90_inq_varid which sets the returned
!>      varid to -1 on error.
!
! CREATION HISTORY:
!       11/03/11  OE  Creation
!F-
!------------------------------------------------------------------------------
  FUNCTION Get_VarID(ncid, name, errmsg)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    INTEGER,                    INTENT(IN) :: ncid
    CHARACTER(LEN=*),           INTENT(IN) :: name
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: errmsg

    ! ---------
    ! Function result
    ! ---------
    INTEGER                      :: Get_VarID

    ! ----------------
    ! Local parameters
    ! ----------------
    CHARACTER( * ),  PARAMETER :: ROUTINE_NAME = 'Get_VarID'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: status

    IF (name == '') THEN
      Get_VarID = -1
      RETURN
    END IF
    status = nf90_inq_varid(ncid, name, Get_VarID)
    IF (status /= nf90_noerr) Get_VarID = -1
    IF (PRESENT(errmsg) .AND. status /= nf90_noerr) THEN
      WRITE(*,*) errmsg//" Unable to read: "//name
    END IF

  END FUNCTION Get_VarID


!------------------------------------------------------------------------------
!F+
! NAME:
!       Find_VarID
!
! PURPOSE:
!>      Search for a netCDF variable with the specified metadata
!
! CATEGORY:
!       NetCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       varid = Find_VarID(ncid, standard_name, axis)
!
! ARGUMENTS:
!>@ARG{ncid, in, INTEGER} netCDF file id
!>@ARG{standard_name, in, CHARACTER(*)\, OPTIONAL} standard_name to match
!>@ARG{standard_name, in, CHARACTER(*)\, OPTIONAL} axis to match
!
! FUNCTION RESULT:
!>@RES{varid, INTEGER} netCDF variable ID or -1 on error
!
! CALLS:
!       nf90_inq_varid
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!       None
!
! PROCEDURE:
!
! CREATION HISTORY:
!       26/10/2018  OE  Creation
!F-
!------------------------------------------------------------------------------
  FUNCTION Find_VarID(ncid, standard_name, axis)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    INTEGER,                    INTENT(IN) :: ncid
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: standard_name
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: axis

    ! ---------
    ! Function result
    ! ---------
    INTEGER                      :: Find_VarID

    ! ----------------
    ! Local parameters
    ! ----------------
    CHARACTER( * ),  PARAMETER :: ROUTINE_NAME = 'Find_VarID'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: status, nvar, varid
    CHARACTER(LEN=60) :: name

    Find_VarID = -1

    status = nf90_inquire(ncid, nVariables=nvar)
    IF (status /= nf90_noerr) RETURN

    DO varid=1, nvar
      IF (PRESENT(standard_name)) THEN
        status = nf90_get_att(ncid, varid, "standard_name", name)
        IF (status /= 0) CYCLE
        IF (TRIM(name) /= standard_name) CYCLE
      END IF
      IF (PRESENT(axis)) THEN
        status = nf90_get_att(ncid, varid, "axis", name)
        IF (status /= 0) CYCLE
        IF (TRIM(name) /= axis) CYCLE
      END IF

      Find_VarID = varid
    END DO

  END FUNCTION Find_VarID


!------------------------------------------------------------------------------
!F+
! NAME:
!       Get_VarInfo
!
! PURPOSE:
!>      Get special CF attributes for variable packing, range etc.
!
! CATEGORY:
!       NetCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       inf = Get_VarInfo(ncid, varid)
!
! ARGUMENTS:
!>@ARG{ncid, in, INTEGER} netCDF file id
!>@ARG{varid, in, INTEGER} netCDF variable id
!
! FUNCTION RESULT:
!>@RES{varinfo, gbcsnetcdf::varinfo} VarInfo structure
!
! CALLS:
!       nf90_get_att
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!       None
!
! PROCEDURE:
!
!
! CREATION HISTORY:
!       17/05/12  OE  Creation
!F-
!------------------------------------------------------------------------------
  FUNCTION Get_VarInfo(ncid, varid) RESULT(inf)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    INTEGER,          INTENT(IN) :: ncid
    INTEGER,          INTENT(IN) :: varid

    ! ---------
    ! Function result
    ! ---------
    TYPE(VarInfo)                :: inf

    ! ----------------
    ! Local parameters
    ! ----------------
    CHARACTER( * ),  PARAMETER :: ROUTINE_NAME = 'Get_VarInfo'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: status
    REAL               :: rvalue
    REAL, DIMENSION(2) :: range

    status = nf90_get_att(ncid, varid, "add_offset", rvalue)
    IF (status == nf90_noerr) THEN
      inf%Is_Packed = .TRUE.
      inf%offset = rvalue
    END IF
    status = nf90_get_att(ncid, varid, "scale_factor", rvalue)
    IF (status == nf90_noerr) THEN
      inf%Is_Packed = .TRUE.
      inf%scale = rvalue
    END IF
    status = nf90_get_att(ncid, varid, "least_significant_digit", rvalue)
    IF (status == nf90_noerr) THEN
      inf%Use_lsd = .TRUE.
      inf%precision = rvalue
    END IF

    status = nf90_get_att(ncid, varid, "valid_range", range)
    IF (status == nf90_noerr) THEN
      inf%Has_Range = .TRUE.
      inf%valid_min = range(1)
      inf%valid_max = range(2)
      inf%value_min = inf%offset + inf%scale*inf%valid_min
      inf%value_max = inf%offset + inf%scale*inf%valid_max
    ELSE
      ! Try checking valid_min / valid_max attributes
      status = nf90_get_att(ncid, varid, "valid_min", rvalue)
      IF (status == nf90_noerr) THEN
        inf%Has_Range = .TRUE.
        inf%valid_min = rvalue
        inf%value_min = inf%offset + inf%scale*inf%valid_min
      END IF
      status = nf90_get_att(ncid, varid, "valid_max", rvalue)
      IF (status == nf90_noerr) THEN
        inf%Has_Range = .TRUE.
        inf%valid_max = rvalue
        inf%value_max = inf%offset + inf%scale*inf%valid_max
      END IF
    END IF

    status = nf90_get_att(ncid, varid, "_FillValue", inf%fillvalue)
    inf%Has_Fill = status == nf90_noerr

  END FUNCTION Get_VarInfo


!------------------------------------------------------------------------------
!F+
! NAME:
!       Get_Var_Dimensions
!
! PURPOSE:
!>      Get the variable dimensions
!
! CATEGORY:
!       NetCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       ndims = Get_Var_Dimensions(ncid, varid, dimids)
!
! ARGUMENTS:
!>@ARG{ncid, in, INTEGER} netCDF file id
!>@ARG{varid, in, INTEGER} netCDF variable id
!>@ARG{dimids, out, INTEGER(:)} dimension ids
!
! FUNCTION RESULT:
!>@RES{ndims, INTEGER} If positive: Number of dimensions; otherwise NetCDF error status
!
! CALLS:
!       nf90_inquire_dimension
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!       None
!
! PROCEDURE:
!
!
! CREATION HISTORY:
!       11/03/11  OE  Creation
!F-
!------------------------------------------------------------------------------
  FUNCTION Get_Var_Dimensions(ncid, varid, dimids)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    INTEGER,                         INTENT(IN)  :: ncid
    INTEGER,                         INTENT(IN)  :: varid
    INTEGER, DIMENSION(max_nc_dims), INTENT(OUT) :: dimids

    ! ---------
    ! Function result
    ! ---------
    INTEGER                             :: Get_Var_Dimensions

    ! ----------------
    ! Local parameters
    ! ----------------
    CHARACTER( * ),  PARAMETER :: ROUTINE_NAME = 'Get_Var_Dimensions'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: ndims

    Get_Var_Dimensions = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
    IF (Get_Var_Dimensions == nf90_noerr) Get_Var_Dimensions = ndims

  END FUNCTION Get_Var_Dimensions


!------------------------------------------------------------------------------
!F+
! NAME:
!       Pack_Data
!
! PURPOSE:
!>      Applys scale/offset transformation for packing data into short or byte
!>      integers. The function rounds to the nearest whole number but does not
!>      convert the type to integer. That should be done by the NetCDF library
!
! CATEGORY:
!       NetCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL Pack_Data(value, inf)
!
! ARGUMENTS:
!>@ARG{value, inout, REAL} value to pack
!>@ARG{inf, in, gbcsnetcdf::varinfo} VarInfo structure
!
! CALLS:
!       None
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!       None
!
! PROCEDURE:
!
! CREATION HISTORY:
!       Written by:   Owen Embury 09/03/2011
!                     IAES, University of Edinburgh
!F-
!------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE Pack_Data_Single(value, inf)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    REAL,          INTENT(INOUT) :: value
    TYPE(VarInfo), INTENT(IN)    :: inf

    IF (inf%Has_Range .AND. inf%Has_Fill) THEN
      IF ((value < inf%value_min) .OR. (value > inf%value_max)) THEN
        value = inf%fillvalue
        RETURN
      END IF
    END IF

    value = ANINT((value - inf%offset) / inf%scale)

  END SUBROUTINE Pack_Data_Single

  ELEMENTAL SUBROUTINE Pack_Data_Double(value, inf)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    REAL(KIND=GbcsDble), INTENT(INOUT) :: value
    TYPE(VarInfo),       INTENT(IN)    :: inf

    IF (inf%Has_Range .AND. inf%Has_Fill) THEN
      IF ((value < inf%value_min) .OR. (value > inf%value_max)) THEN
        value = inf%fillvalue
        RETURN
      END IF
    END IF

    value = ANINT((value - inf%offset) / inf%scale)

  END SUBROUTINE Pack_Data_Double


!------------------------------------------------------------------------------
!F+
! NAME:
!       Unpack_Data
!
! PURPOSE:
!>      Applys scale/offset transformation for unpacking data from short or
!>      byte integers.
!
! CATEGORY:
!       NetCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL = Unpack_Data(value, inf)
!
! ARGUMENTS:
!>@ARG{value, inout, REAL} packed data
!>@ARG{inf, in, gbcsnetcdf::varinfo} VarInfo structure
!
! CALLS:
!       None
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!       None
!
! PROCEDURE:
!
! CREATION HISTORY:
!       Written by:   Owen Embury 09/03/2011
!                     IAES, University of Edinburgh
!F-
!------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE Unpack_Data_Single(value, inf)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    REAL,          INTENT(INOUT) :: value
    TYPE(VarInfo), INTENT(IN)    :: inf

    IF (inf%Has_Range .AND. &
        ((value < inf%valid_min) .OR. (value > inf%valid_max))) THEN
      RETURN
    ELSE
      value = inf%offset + (value * inf%scale)
    END IF

  END SUBROUTINE Unpack_Data_Single

  ELEMENTAL SUBROUTINE Unpack_Data_Double(value, inf)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    REAL(KIND=GbcsDble), INTENT(INOUT) :: value
    TYPE(VarInfo),       INTENT(IN)    :: inf

    IF (inf%Has_Range .AND. &
        ((value < inf%valid_min) .OR. (value > inf%valid_max))) THEN
      RETURN
    ELSE
      value = inf%offset + (value * inf%scale)
    END IF

  END SUBROUTINE Unpack_Data_Double


!------------------------------------------------------------------------------
!S+
! NAME:
!       Get_Var_Ranges
!
! PURPOSE:
!>      Internal subroutine for calculating the nc_start, nc_count arguments
!>      when for read/writting a subset of a netCDF variable to an array which
!>      may have fewer dimensions.
!>
!>      start and count refer to the in memory array which may have fewer
!>      dimensions than the netCDF variable.
!>      start is assumed to correspond to the slowest varying dimensions
!>      count is assumed to correspond to the fastest varying dimensions
!
! CATEGORY:
!       NetCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL Get_Var_Ranges(start, count, nc_ndims, nc_start, nc_count, nc_map, remap)
!
! ARGUMENTS:
!>@ARG{start, in, INTEGER(:)}
!>@ARG{count, in, INTEGER(:)}
!>@ARG{nc_ndims, in, INTEGER}
!>@ARG{nc_start, out, INTEGER(max_nc_dims)}
!>@ARG{nc_count, out, INTEGER(max_nc_dims)}
!>@ARG{nc_map, out, INTEGER(max_nc_dims)}
!>@ARG{remap, in, INTEGER(:)}
!
! CALLS:
!       None
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!       None
!
! PROCEDURE:
!
! CREATION HISTORY:
!       Written by:   Owen Embury 14/01/2013
!                     IAES, University of Edinburgh
!S-
!------------------------------------------------------------------------------
  PURE SUBROUTINE Get_Var_Ranges(start, count, nc_ndims, nc_start, nc_count, nc_map, remap)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    INTEGER, DIMENSION(:),           INTENT(IN)  :: start
    INTEGER, DIMENSION(:),           INTENT(IN)  :: count
    INTEGER,                         INTENT(IN)  :: nc_ndims
    INTEGER, DIMENSION(max_nc_dims), INTENT(OUT) :: nc_start
    INTEGER, DIMENSION(max_nc_dims), INTENT(OUT) :: nc_count
    INTEGER, DIMENSION(max_nc_dims), INTENT(OUT) :: nc_map
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)  :: remap

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: i, ndims

    nc_start        = 1
    IF (SIZE(start) > nc_ndims) THEN
      ! Number of dimensions of memory array exceeds that of netCDF variable
      nc_start(1:nc_ndims) = start(1:nc_ndims)
    ELSE
      nc_start(nc_ndims+1-SIZE(start):nc_ndims) = start
    END IF

    ndims = SIZE(count)
    nc_count        = 1
    nc_count(1:ndims)   = count

    nc_map          = 1
    nc_map(1:ndims) = (/1, (PRODUCT(count(:i)), i=1, ndims-1) /)

    IF (PRESENT(remap)) THEN
      ndims = SIZE(remap)
      nc_start(1:ndims) = nc_start(remap)
      nc_count(1:ndims) = nc_count(remap)
      nc_map(1:ndims)   = nc_map(remap)
    END IF

  END SUBROUTINE Get_Var_Ranges


!------------------------------------------------------------------------------
!F+
! NAME:
!       Write_Var
!
! PURPOSE:
!>      Writes data out to NetCDF file. Applying data packing if required.
!>      Intended to write "chunks" of data out to a file with one "unlimited"
!>      dimension and takes a single positional argument which is the starting
!>      index along the "unlimited" dimension.
!
! CATEGORY:
!
! LANGUAGE:
!       Fortran-2003
!
! CALLING SEQUENCE:
!       status = Write_Var(ncid, varid, values, start)
!
! ARGUMENTS:
!>@ARG{ncid, in, INTEGER} NetCDF file id
!>@ARG{varid, in, INTEGER} NetCDF variable id
!>@ARG{values, in, } data to store
!>@ARG{start, in, INTEGER\, DIMENSION(:)} index along final dimensions to start write
!
! RESULT:
!>@RES{status, INTEGER} NetCDF error status
!
! CALLS:
!       Get_Var_Ranges
!       Get_VarInfo
!       Quantize
!       Pack_Data
!       nf90_inquire_variable
!       nf90_put_var
!
! SIDE EFFECTS:
!>@SIDE Calls NetCDF library to output data
!
! RESTRICTIONS:
!       Size of values must not exceed size of NetCDF variable
!
! PROCEDURE:
!
! CREATION HISTORY:
!       Written by:   Owen Embury 26/09/2008
!                     IAES, University of Edinburgh
!       09/03/2011:   Moved to GbcsNetCDF module and added called to Is_Packed
!                     and Pack_Data functions.
!       26/10/2018:   Add inplace quantization (requires Fortran-2003)
!F-
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!F+
! NAME:
!       Write_Var1F
!
! PURPOSE:
!       Write a 1d array of Single precision REALs to NetCDF variable
  FUNCTION Write_Var1F(ncid, varid, values, start ) RESULT(status)
    IMPLICIT NONE

    ! ---------
    ! Arguments
    ! ---------
    INTEGER,                               INTENT(IN) :: ncid
    INTEGER,                               INTENT(IN) :: varid
    REAL(KIND=GbcsReal), DIMENSION(:),     INTENT(IN) :: values
    INTEGER,             DIMENSION(:),     INTENT(IN) :: start

    ! ---------
    ! Function result
    ! ---------
    INTEGER               :: status

    ! ----------------
    ! Local parameters
    ! ----------------
    CHARACTER( * ),        PARAMETER :: ROUTINE_NAME = 'Write_Var1F'

    ! ---------------
    ! Local variables
    ! ---------------
    TYPE(VarInfo)         :: inf
    INTEGER, DIMENSION(max_nc_dims) :: nc_start, nc_count, nc_map
    INTEGER               :: ndims
    REAL, ALLOCATABLE, DIMENSION(:) :: buffer


    status = nf90_inquire_variable(ncid, varid, ndims=ndims)
    IF (status /= nf90_noerr) RETURN

    CALL Get_Var_Ranges(start, SHAPE(values), ndims, nc_start, nc_count, nc_map)
    inf = Get_VarInfo(ncid, varid)

    IF (inf%Is_Packed .OR. inf%Use_lsd) THEN
      buffer = values
      IF (inf%Use_lsd)   CALL Quantize(buffer, inf%precision)
      IF (inf%Is_Packed) CALL Pack_Data(buffer, inf)
      status = nf90_put_var(ncid, varid, buffer, nc_start(1:ndims), nc_count(1:ndims))
    ELSE
      status = nf90_put_var(ncid, varid, values, nc_start(1:ndims), nc_count(1:ndims))
    END IF

  END FUNCTION Write_Var1F

!------------------------------------------------------------------------------
!F+
! NAME:
!       Write_Var1I
!
! PURPOSE:
!       Write a 1d array of INTEGERs to NetCDF variable
  FUNCTION Write_Var1I(ncid, varid, values, start ) RESULT(status)
    IMPLICIT NONE

    ! ---------
    ! Arguments
    ! ---------
    INTEGER,                               INTENT(IN) :: ncid
    INTEGER,                               INTENT(IN) :: varid
    INTEGER,             DIMENSION(:),     INTENT(IN) :: values
    INTEGER,             DIMENSION(:),     INTENT(IN) :: start

    ! ---------
    ! Function result
    ! ---------
    INTEGER               :: status

    ! ----------------
    ! Local parameters
    ! ----------------
    CHARACTER( * ),        PARAMETER :: ROUTINE_NAME = 'Write_Var1I'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER, DIMENSION(max_nc_dims) :: nc_start, nc_count, nc_map
    INTEGER               :: ndims


    status = nf90_inquire_variable(ncid, varid, ndims=ndims)
    IF (status /= nf90_noerr) RETURN

    CALL Get_Var_Ranges(start, SHAPE(values), ndims, nc_start, nc_count, nc_map)

    status = nf90_put_var(ncid, varid, values, nc_start(1:ndims), nc_count(1:ndims))

  END FUNCTION Write_Var1I

!------------------------------------------------------------------------------
!F+
! NAME:
!       Write_Var2D
!
! PURPOSE:
!       Write a 2d array of Double precision REALs to NetCDF variable
  FUNCTION Write_Var2D(ncid, varid, values, start ) RESULT(status)
    IMPLICIT NONE

    ! ---------
    ! Arguments
    ! ---------
    INTEGER,                               INTENT(IN) :: ncid
    INTEGER,                               INTENT(IN) :: varid
    REAL(KIND=GbcsDble), DIMENSION(:,:),   INTENT(IN) :: values
    INTEGER,             DIMENSION(:),     INTENT(IN) :: start

    ! ---------
    ! Function result
    ! ---------
    INTEGER               :: status

    ! ----------------
    ! Local parameters
    ! ----------------
    CHARACTER( * ),        PARAMETER :: ROUTINE_NAME = 'Write_Var2D'

    ! ---------------
    ! Local variables
    ! ---------------
    TYPE(VarInfo)         :: inf
    INTEGER, DIMENSION(max_nc_dims) :: nc_start, nc_count, nc_map
    INTEGER               :: ndims
    REAL(KIND=GbcsDble), ALLOCATABLE, DIMENSION(:,:) :: buffer


    status = nf90_inquire_variable(ncid, varid, ndims=ndims)
    IF (status /= nf90_noerr) RETURN

    CALL Get_Var_Ranges(start, SHAPE(values), ndims, nc_start, nc_count, nc_map)
    inf = Get_VarInfo(ncid, varid)

    IF (inf%Is_Packed .OR. inf%Use_lsd) THEN
      buffer = values
      IF (inf%Use_lsd)   CALL Quantize(buffer, inf%precision)
      IF (inf%Is_Packed) CALL Pack_Data(buffer, inf)
      status = nf90_put_var(ncid, varid, buffer, nc_start(1:ndims), nc_count(1:ndims))
    ELSE
      status = nf90_put_var(ncid, varid, values, nc_start(1:ndims), nc_count(1:ndims))
    END IF

  END FUNCTION Write_Var2D

!------------------------------------------------------------------------------
!F+
! NAME:
!       Write_Var2F
!
! PURPOSE:
!       Write a 2d array of Single precision REALs to NetCDF variable
  FUNCTION Write_Var2F(ncid, varid, values, start ) RESULT(status)
    IMPLICIT NONE

    ! ---------
    ! Arguments
    ! ---------
    INTEGER,                               INTENT(IN) :: ncid
    INTEGER,                               INTENT(IN) :: varid
    REAL(KIND=GbcsReal), DIMENSION(:,:),   INTENT(IN) :: values
    INTEGER,             DIMENSION(:),     INTENT(IN) :: start

    ! ---------
    ! Function result
    ! ---------
    INTEGER               :: status

    ! ----------------
    ! Local parameters
    ! ----------------
    CHARACTER( * ),        PARAMETER :: ROUTINE_NAME = 'Write_Var2F'

    ! ---------------
    ! Local variables
    ! ---------------
    TYPE(VarInfo)         :: inf
    INTEGER, DIMENSION(max_nc_dims) :: nc_start, nc_count, nc_map
    INTEGER               :: ndims
    REAL, ALLOCATABLE, DIMENSION(:,:) :: buffer


    status = nf90_inquire_variable(ncid, varid, ndims=ndims)
    IF (status /= nf90_noerr) RETURN

    CALL Get_Var_Ranges(start, SHAPE(values), ndims, nc_start, nc_count, nc_map)
    inf = Get_VarInfo(ncid, varid)

    IF (inf%Is_Packed .OR. inf%Use_lsd) THEN
      buffer = values
      IF (inf%Use_lsd)   CALL Quantize(buffer, inf%precision)
      IF (inf%Is_Packed) CALL Pack_Data(buffer, inf)
      status = nf90_put_var(ncid, varid, buffer, nc_start(1:ndims), nc_count(1:ndims))
    ELSE
      status = nf90_put_var(ncid, varid, values, nc_start(1:ndims), nc_count(1:ndims))
    END IF

  END FUNCTION Write_Var2F

!------------------------------------------------------------------------------
!F+
! NAME:
!       Write_Var2I
!
! PURPOSE:
!       Write a 2d array of INTEGERs to NetCDF variable
  FUNCTION Write_Var2I(ncid, varid, values, start ) RESULT(status)
    IMPLICIT NONE

    ! ---------
    ! Arguments
    ! ---------
    INTEGER,                               INTENT(IN) :: ncid
    INTEGER,                               INTENT(IN) :: varid
    INTEGER,             DIMENSION(:,:),   INTENT(IN) :: values
    INTEGER,             DIMENSION(:),     INTENT(IN) :: start

    ! ---------
    ! Function result
    ! ---------
    INTEGER               :: status

    ! ----------------
    ! Local parameters
    ! ----------------
    CHARACTER( * ),        PARAMETER :: ROUTINE_NAME = 'Write_Var2I'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER, DIMENSION(max_nc_dims) :: nc_start, nc_count, nc_map
    INTEGER               :: ndims


    status = nf90_inquire_variable(ncid, varid, ndims=ndims)
    IF (status /= nf90_noerr) RETURN

    CALL Get_Var_Ranges(start, SHAPE(values), ndims, nc_start, nc_count, nc_map)

    status = nf90_put_var(ncid, varid, values, nc_start(1:ndims), nc_count(1:ndims))

  END FUNCTION Write_Var2I

!------------------------------------------------------------------------------
!F+
! NAME:
!       Write_Var2S
!
! PURPOSE:
!       Write a 2d array of 2-byte INTEGERs to NetCDF variable
  FUNCTION Write_Var2S(ncid, varid, values, start ) RESULT(status)
    IMPLICIT NONE

    ! ---------
    ! Arguments
    ! ---------
    INTEGER,                               INTENT(IN) :: ncid
    INTEGER,                               INTENT(IN) :: varid
    INTEGER(KIND=GbcsInt2), DIMENSION(:,:),INTENT(IN) :: values
    INTEGER,             DIMENSION(:),     INTENT(IN) :: start

    ! ---------
    ! Function result
    ! ---------
    INTEGER               :: status

    ! ----------------
    ! Local parameters
    ! ----------------
    CHARACTER( * ),        PARAMETER :: ROUTINE_NAME = 'Write_Var2S'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER, DIMENSION(max_nc_dims) :: nc_start, nc_count, nc_map
    INTEGER               :: ndims


    status = nf90_inquire_variable(ncid, varid, ndims=ndims)
    IF (status /= nf90_noerr) RETURN

    CALL Get_Var_Ranges(start, SHAPE(values), ndims, nc_start, nc_count, nc_map)

    status = nf90_put_var(ncid, varid, values, nc_start(1:ndims), nc_count(1:ndims))

  END FUNCTION Write_Var2S

!------------------------------------------------------------------------------
!F+
! NAME:
!       Write_Var2B
!
! PURPOSE:
!       Write a 2d array of 1-byte INTEGERs to NetCDF variable
  FUNCTION Write_Var2B(ncid, varid, values, start ) RESULT(status)
    IMPLICIT NONE

    ! ---------
    ! Arguments
    ! ---------
    INTEGER,                               INTENT(IN) :: ncid
    INTEGER,                               INTENT(IN) :: varid
    INTEGER(KIND=GbcsInt1), DIMENSION(:,:),INTENT(IN) :: values
    INTEGER,             DIMENSION(:),     INTENT(IN) :: start

    ! ---------
    ! Function result
    ! ---------
    INTEGER               :: status

    ! ----------------
    ! Local parameters
    ! ----------------
    CHARACTER( * ),        PARAMETER :: ROUTINE_NAME = 'Write_Var2B'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER, DIMENSION(max_nc_dims) :: nc_start, nc_count, nc_map
    INTEGER               :: ndims


    status = nf90_inquire_variable(ncid, varid, ndims=ndims)
    IF (status /= nf90_noerr) RETURN

    CALL Get_Var_Ranges(start, SHAPE(values), ndims, nc_start, nc_count, nc_map)

    status = nf90_put_var(ncid, varid, values, nc_start(1:ndims), nc_count(1:ndims))

  END FUNCTION Write_Var2B


!------------------------------------------------------------------------------
!F+
! NAME:
!       Quantize
!
! PURPOSE:
!>      Round a REAL number to the nearest power of two which provides the
!>      specified precision.
!
! CATEGORY:
!       NetCDF
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL Quantize(value, precision)
!
! ARGUMENTS:
!>@ARG{value, inout, REAL} Number to quantize
!>@ARG{precision, in, REAL} Desired precision (e.g. 0.1)
!
! CALLS:
!       None
!
! SIDE EFFECTS:
!       None
!
! RESTRICTIONS:
!       None
!
! PROCEDURE:
!
! CREATION HISTORY:
!       19/12/2017  OE  Creation
!F-
!------------------------------------------------------------------------------
  PURE SUBROUTINE Quantize_Scalar(value, precision)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    REAL, INTENT(INOUT) :: value
    REAL, INTENT(IN)    :: precision

    ! ----------------
    ! Local parameters
    ! ----------------
    REAL, PARAMETER :: ilog2 = -1.0 / LOG(2.0)

    ! ---------------
    ! Local variables
    ! ---------------
    REAL :: scale

    scale = 2.0 ** CEILING(ilog2 * LOG(precision))

    IF (1/scale > SPACING(value)) THEN
      value = ANINT(value*scale)/scale
    END IF

  END SUBROUTINE Quantize_Scalar

  PURE SUBROUTINE Quantize_1F(array, precision)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    REAL, DIMENSION(:), INTENT(INOUT) :: array
    REAL,               INTENT(IN)    :: precision

    ! ----------------
    ! Local parameters
    ! ----------------
    REAL, PARAMETER :: ilog2 = -1.0 / LOG(2.0)

    ! ---------------
    ! Local variables
    ! ---------------
    REAL :: scale
    INTEGER :: i

    scale = 2.0 ** CEILING(ilog2 * LOG(precision))

    DO i=1, SIZE(array)
      IF (1/scale > SPACING(array(i))) THEN
        array(i) = ANINT(array(i)*scale)/scale
      END IF
    END DO

  END SUBROUTINE Quantize_1F

  PURE SUBROUTINE Quantize_2F(array, precision)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    REAL, DIMENSION(:,:), INTENT(INOUT) :: array
    REAL,                 INTENT(IN)    :: precision

    ! ----------------
    ! Local parameters
    ! ----------------
    REAL, PARAMETER :: ilog2 = -1.0 / LOG(2.0)

    ! ---------------
    ! Local variables
    ! ---------------
    REAL :: scale
    INTEGER :: i, j

    scale = 2.0 ** CEILING(ilog2 * LOG(precision))

    DO j=1, SIZE(array,2)
      DO i=1, SIZE(array,1)
        IF (1/scale > SPACING(array(i,j))) THEN
          array(i,j) = ANINT(array(i,j)*scale)/scale
        END IF
      END DO
    END DO

  END SUBROUTINE Quantize_2F

  PURE SUBROUTINE Quantize_2D(array, precision)
    IMPLICIT NONE
    ! ---------
    ! Arguments
    ! ---------
    REAL(KIND=GbcsDble), DIMENSION(:,:), INTENT(INOUT) :: array
    REAL,                                INTENT(IN)    :: precision

    ! ----------------
    ! Local parameters
    ! ----------------
    REAL(KIND=GbcsDble), PARAMETER :: ilog2 = -1.0 / LOG(2.0)

    ! ---------------
    ! Local variables
    ! ---------------
    REAL(KIND=GbcsDble) :: scale
    INTEGER :: i, j

    scale = 2.0 ** CEILING(ilog2 * LOG(precision))

    DO j=1, SIZE(array,2)
      DO i=1, SIZE(array,1)
        IF (1/scale > SPACING(array(i,j))) THEN
          array(i,j) = ANINT(array(i,j)*scale)/scale
        END IF
      END DO
    END DO

  END SUBROUTINE Quantize_2D

END MODULE GbcsNetCDF
