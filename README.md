# RAVen
RAVen project data processing.
Data processing steps in short
1.	Prepare the merged HDF5 files
1.1.	The HDF5 files need to be merged by ohfa.py module to have all variables in one file and with the right mask (on- or off-shore, or without any mask for the full domain)
The module ohfa.py is an Odim Hdf5 file Aggregator, which aggregates two or more Odim HDF5 files into a single one. It was programmed by Christophe Ferauge (2014) and modified for RAVen project by Maryna Lukach (2017). 
The modified module ohfa.py does the merging and applies the land mask if it is needed. The -o and -c options are added and masking of on- or off- shore parts of the data is done based on the mask provided in -c option. In the processing the output was placed in onshore, offshore or full subdirectories of h5.merged directory of the main directory.

The call of the ohfa.py module:
ohfa.py
with correct parameters and input files.
Input:  volume HDF5 files (dBZ, V, W) and the directory for the output
-v - verbose
-d - debug
-o – offshore (store True)
-c – coastline (string) the full path and the name of the land-mask generated based on the coastlines.
Output: one HDF5 file with all data placed in separate datasets will be written to the provided output directory
An example of the call:
python ohfa.py -v -d -o -c ./bejab_landcoast_validmask.hdf ./bin/2016101401450500dBZ.vol.h5   ./bin/2016101401450500V.vol.h5 ./bin/2016101401450500W.vol.h5  ./h5.merged/offshore/


The run will combine the dBZ, V and W data from the input files in one HDF5 file, apply the mask and write the output file in the ./h5.merged/offshore/ directory.

2.	Run the Vol2Bird algorithm over merged hdf5 files.
The algorithm takes as input all the files from the directory provided by –i parameter. The output will be placed in the directory defined by –o parameter. Additional –t parameter gives a possibility to run the module only for the files in –i directory that start with a given YYYYmmDD (applicable only for the files having the timestamp at the beginning of their name).
In the processing the outputs of the Vol2Bird algorithm were placed in full_vp, onshore_vp or offshore_vp subdirectories of directory h5.merged of the main directory.
The processing for RAVen was done with the off- or on-shore part set to "nodata" and with the intern Vol2Bird parameter REQUIRE_VRAD set to True. (used options from the option.comfig file can be found in the apendix)
The attribute "nodata" in OPERA data model is defined as "The raw value used to denote areas void of data (never radiated)". The parameter REQUIRE_VRAD makes Vol2Bird require from a range gate to have a valid radial velocity to contribute to the calculation.
With these settings, the total number of points (on- + off- shore), used for calculation of the bird reflectivity for each altitude level, equals the number of points used for the calculation on the whole domain.
This number of points can be found in "/dataset1/ data6" group of the Vol2Bird output.
The call of the vol2bird.py module with correct parameters and input files:
python vol2bird.py
Input: all input files in the directory provided in –i parameter
–i - the input directory
–o – the output directory
–t - (YYYYmmDD format) -  gives a possibility to run the module only for the files that start with a given YYYYmmDD (applicable only for the files having the timestamp at the beginning of their name).
Output: bird’s vertical profiles
	An example of the call:
python vol2bird.py -i ./h5.merged/onshore -o ./h5.merged/onshore_vp -t 201609
The algorithm will proceed through the onshore data of September 2016.

3.	Calculate MTR’s for vertical profiles from vol2bird algorithm
The module mtr_vbird.py calculates “specific MTRs “ – MTRs per altitude level for 20 levels between 200m and 4000m a.m.s.l. The MTR per level is calculated as density*speed*3.6/5, where for onshore and offshore domains the speed is taken from the full domain and the density is calculated from the birds reflectivity values. These values are non-zero for all three datasets (full, on-shore, off-shore) and the density calculation is based on the same assumption of 11 cm2 bird cross-section.
The MTR [birds/ km/ hour] is calculated from the volume density [bids/km3], where the volume density is calculated the same way as before.
It is the bird reflectivity (eta [cm2/km3]) output of the vol2bird algorithm ("/dataset1/ data7") divided by the bird cross-section 11 cm2.
The speed used for all datasets (on-, off-shore, and full) is the speed from the full domain output.
For the altitude intervals between 200m and 800m, between 800m and 1400m, between 200m and 1400m the total MTR per interval and the mean speed are calculated. Based on these “specific MTRs” all mean MTRs are calculated. The output of the module is saved to the panda’s datasets and dumped to csv files.
In the processing the outputs of mtr_vbird algorithm were placed in MTR subdirectory of the main directory. The input for the module should be in onshore_vp, offshore_vp and full_vp subdirectories of the directory h5.merged of the main directory.
The hourly mean values can be found in bejab_vbird_avgmtralllevels_YYYYmm_*_vp.csv files.
python mtr_vbird.py
Input: all vertical profiles in the directory provided in –i parameter
–i - the input directory
–o – the output directory
–t - (YYYYmmDD format) -  gives a possibility to run the module only for the files with a given timestamp in the name of the file.
Output: .csv files with calculated means.
This module requires some additional python libraries:
Pandas
An example of the call:
python mtr_vbird.py -t 201609 -i  ./h5.merged/offshore_vp/ -o ./MTR/
