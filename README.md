# slstr-preprocessor

Fortran library code and test driver program for aggregating visual channels in SLSTR scenes to the IR grid

* The preprocessor works by first computing a neighbourhood of visual pixels in a scene that are closest to each IR pixel
  
  * Neighbourhoods are defined using the K nearest neighbours that lie within a maximum radius
  * Cartesian distances are used, based on along-track and across-track offsets of each pixel
  * Cosmetically filled pixels are ignored to avoid making a double contribution to the output
  * Choose from using the a-stripe or (where available) the b-stripe or a combination of a-stripe and b-stripe
    
* For each visual band, IR compatible rasters are then computed from all pixels which are not missing or cosmetically filled within
  the neighbourhood of each IR pixel, using the specified aggregation function(s).

  * If no suitable pixels appear in the neighbourhood, a missing value is used
  * MEAN, MAX, SD (population standard deviation) and MAX-MIN aggregations are supported and other ones can easily be added
    
* For comparison purposes, a `simple` mode is provided which ignores the along and across track distances and aggregates the visible
  band value for an IR pixel at indexes (x,y) by agggregating the 4 visible pixels at indexes (2x,2y),(2x+1,2y),(2x,2y+1),(2x+1,2y+1)
  
This software processes L1 scenes from the Sentinel 3A/3B SLSTR satellites as described here:

https://sentinel.esa.int/documents/247904/1872792/Sentinel-3-SLSTR-Product-Data-Format-Specification-Level-1

## Fortran module and test program

A [fortran integration module](SLSTR_Preprocessor.f90) should be linked with a client program to integrate this functionality.

Usage of this module can be described by the following code snippet:

```f90
USE slstr_preprocessor

! The module variable MISSING_R must first be set to represent missing values in array data returned by the module
! For example:
MISSING_R = -1.0e+30

! The module variable MAX_NEIGHBOUR_DISTANCE (defaults to 10000) - sets the maximum distance in metres a VIS pixel 
! can lie from the IR pixel center to be included in the neighbourhood calculation
! For example:
MAX_NEIGHBOUR_DISTANCE = 5000 

! Declare variable to hold the neighbourhood, the array will be allocated inside the compute_scene_neighbourhood
TYPE(NEIGHBOURHOOD_MAP), ALLOCATABLE, DIMENSION(:,:) :: neighbourhood
! Build neighbourhood mapping for the a-stripe, nadir view, opening files as needed from scene folder
CALL compute_scene_neighbourhood('n','/path/to/scene',neighbourhood,'a')

! Process visual band 4 from the nadir view, using the constructed neighbourhood mapping, and applying the mean function
! Possible functions are represented by module constants FUNCTION_MEAN, FUNCTION_MAX, FUNCTION_MIN_MAX_DIFF, FUNCTION_SD
! vis_output_radiance is a 2D array pre-allocated to receive the output.  Use up to the 5 nearest neighbours.
! use width=1500,height=1200 for nadir view
! use width=900,height=1200 for oblique view
REAL, ALLOCATABLE, DIMENSION(:,:) :: vis_output_radiance
ALLOCATE(vis_output_radiance(1500,1200))
CALL process_scene_band('n','/path/to/scene',4,vis_output_radiance,neighbourhood,5,FUNCTION_MEAN)

! Process further bands from the same scene and nadir view, by calling process_scene_band further times

! To process the scene band without using neighbourhood map (using the simple method) do not build the neighbourhood 
! The simple method agggregates pixels at location (2x,2y),(2x+1,2y),(2x,2y+1),(2x+1,2y+1)
TYPE(NEIGHBOURHOOD_MAP) :: empty_neighbourhood
REAL, ALLOCATABLE, DIMENSION(:,:) :: vis_output_radiance
ALLOCATE(vis_output_radiance(1500,1200))
CALL process_scene_band('n','/path/to/scene',4,vis_output_radiance,empty_neighbourhood,0,FUNCTION_MEAN)
```

In the example above, a neighbourhood was constructed on the pixels in the a-stripe.  A neighbourhood can similarly be constructed on the b-stripe:

```f90
CALL compute_scene_neighbourhood('n','/path/to/scene',neighbourhood,'b')
```

It is also possible to create a combined neighbourhood that includes the closest pixels drawn from both the a- and b-stripes, 
using the `merge_neighborhoods` subroutine:

```f90
TYPE(NEIGHBOURHOOD_MAP), ALLOCATABLE, DIMENSION(:,:) :: neighbourhood_a
TYPE(NEIGHBOURHOOD_MAP), ALLOCATABLE, DIMENSION(:,:) :: neighbourhood_ab
CALL compute_scene_neighbourhood('n','/path/to/scene',neighbourhood_a,'a')
CALL compute_scene_neighbourhood('n','/path/to/scene',neighbourhood_ab,'b')
CALL merge_neighbourhoods(neighbourhood_a,neighbourhood_ab)  ! the combined a/b stripe neighborhood is then stored in neighbourhood_ab
```

Note that the behaviour of `process_scene_band` using a neighbourhood which includes b-stripe data will raise an error and cause the program to stop,
when called for a visible channel which does not include a b-stripe in the input data.

## Customising compile time options

Edit the module parameters of [fortran integration module](SLSTR_Preprocessor.f90) to customise its behaviour

The most commonly used parameter sets the maximum limit on the number of neighbours that can used (the `effective_k` parameter passed to `process_scene_band`)

* MAX_K_NEAREST_NEIGHBOURS (defaults to 10) - controls how many neighbours (K) at most can be aggegrated to compute an output pixel

```f90
!> Maximum neighbourhood size.
INTEGER, PARAMETER :: MAX_K_NEAREST_NEIGHBOURS = 10
```

## Building the code and test program

A [standalone test program Preprocess_SLSTR](Preprocess_SLSTR.f90) which can be linked with this module is also provided.

To build the test program clone this repo, cd into its root directory and then run:

```
make clean
make
```

This should build an executable `Preprocess_SLSTR` in the current folder.  The Makefile expects:

* gfortran to be installed
* netcdf fortran libraries (libnetcdff) installed to `/usr/local/lib` (libraries) and `/usr/local/include` (headers)

## Running the module via the test program Preprocess_SLSTR

To run the test program, supply the input folder containing an SLSTR scene 
and specify the output folder into which the output files S*_radiance_in.nc and S*radiance_io.nc are written.

* first, unzip the scene if zipped.  In the example below, the file is unzipped into the folder containing the Preprocess_SLSTR executable.

```
unzip S3A_SL_1_RBT____20200101T235811_20200102T000111_20200103T033745_0179_053_187_3420_LN2_O_NT_003.zip
```

* run with neighbourhood calculation, output files to /tmp:

```
./Preprocess_SLSTR S3A_SL_1_RBT____20200101T235811_20200102T000111_20200103T033745_0179_053_187_3420_LN2_O_NT_003.SEN3 /tmp
```

* run without neighbourhood calculation (old behaviour):

```
./Preprocess_SLSTR S3A_SL_1_RBT____20200101T235811_20200102T000111_20200103T033745_0179_053_187_3420_LN2_O_NT_003.SEN3 /tmp 
         --simple
```

* run with neighbourhood calculation and print basic neighbourhood statistics:

```
./Preprocess_SLSTR S3A_SL_1_RBT____20200101T235811_20200102T000111_20200103T033745_0179_053_187_3420_LN2_O_NT_003.SEN3 /tmp 
         --stats
```

* run with neighbourhood calculation, using K=6 and setting the maximum distance of a neighbour to 4000m:

```
./Preprocess_SLSTR S3A_SL_1_RBT____20200101T235811_20200102T000111_20200103T033745_0179_053_187_3420_LN2_O_NT_003.SEN3 /tmp 
          --effective_k 6 --max_distance 4000
```

## Output files and variables

For each visible channel C (1-6), a pair of output files are generated with the data regridded on the IR nadir and oblique angle grids, for example for channel 1 these are named:

* `S1_radiance_in.nc` - for nadir view
* `S1_radiance_io.nc` - for oblique view

Inside each file, with the default configuration, variables will be created for each aggregation function supported.

For example, in `S1_radiance_in.nc` the following variables will be created:

* `S1_radiance_max` - value is the maximum pixel value within the neighbourhood
* `S1_radiance_mean` - value is the mean pixel value within the neighbourhood
* `S1_radiance_min_max_diff` - value is the difference between the max and min values within the neighbourhood
* `S1_radiance_sd` - value is the population standard deviation of values in the neighbourhood

