# slstr-preprocessor

Fortran library code and test driver program for aggregating visual channels in SLSTR scenes to the IR grid

* The preprocessor works by first computing a neighbourhood of visual pixels in a scene that are closest to each IR pixel
  
  * Neighbourhoods can be defined using the N nearest neighbours that lie within a maximum radius
  * Cartesian distances are used, based on along-track and across-track offsets of each pixel
  * Cosmetically filled pixels are ignored to avoid making a double contribution to the output
    
* For each visual band, IR compatible rasters are then computed from all pixels which are not nan or cosmetically filled within
  the neighbourhood of each IR pixel, using the specified aggregation function(s).

  * If no suitable pixels appear in the neighbourhood, a nan value is used
  * MEAN, MAX, SD and MAX-MIN aggregations are supported and other ones can easily be added
    
* For comparison purposes, a `simple` mode is provided which ignores the along and across track distances and aggregates the visible
  band value for an IR pixel at (x,y) by simply taking the mean of the 4 visible bands at (2x,2y),(2x+1,2y),(2x,2y+1),(2x+1,2y+1)
  
This software processes L1 scenes from the Sentinel 3 SLSTR satellites as described here:

https://sentinel.esa.int/documents/247904/1872792/Sentinel-3-SLSTR-Product-Data-Format-Specification-Level-1

## Fortran module and test program

[fortran integration module](SLSTR_Preprocessor.f90) should be linked with a client program to integrate this functionality.

Usage of this module can be described by the following pattern:

```
USE slstr_preprocessor

! Declare variable to hold the neighbourhood, the array will be allocated inside the compute_scene_neighbourhood
TYPE(NEIGHBOURHOOD_MAP), ALLOCATABLE, DIMENSION(:,:) :: neighbourhood
! Build neighbourhood mapping for nadir, opening files as needed from scene folder
! This should eventually work with a folder or path to a compressed file, right now has to be a folder
CALL compute_scene_neighbourhood('n','/path/to/scene',neighbourhood)

! Process visual band 4 from the nadir view
! vis_output_radiance is a 2D array pre-allocated to receive the output
REAL, ALLOCATABLE, DIMENSION(:,:) :: vis_output_radiance
ALLOCATE(vis_output_radiance(ir_width,ir_height))
CALL process_scene_band('n','/path/to/scene',4,vis_output_radiance,neighbourhood)

! Process further bands from the same scene and nadir view, by calling process_scene_band further times

! To process the scene band without using neighbourhood map (using old behaviour) omit the optional neighbourhood
REAL, ALLOCATABLE, DIMENSION(:,:) :: vis_output_radiance
ALLOCATE(vis_output_radiance(ir_width,ir_height))
CALL process_scene_band('n','/path/to/scene',4,vis_output_radiance)
```

A [standalone test program "scene_tester"](Preprocess_SLSTR.f90) which can be linked with this module to test it is also provided.

To build the test program clone this repo and then use:

```
make clean
make
```

This should build an executable `Preprocess_SLSTR` in the current folder.  To run the test program, supply the input folder containing an SLSTR scene 
and the output folder into which the output files S*_radiance_in.nc and S*radiance_io.nc are written:

* run with neighbourhood calculation:

```
./Preprocess_SLSTR path/to/data/S3A_SL_1_RBT____20200101T235811_20200102T000111_20200103T033745_0179_053_187_3420_LN2_O_NT_003.SEN3 /tmp
```

* run without neighbourhood calculation (old behaviour):

```
./Preprocess_SLSTR path/to/data/S3A_SL_1_RBT____20200101T235811_20200102T000111_20200103T033745_0179_053_187_3420_LN2_O_NT_003.SEN3 /tmp simple
```

