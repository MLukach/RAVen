import sys
import warnings
import os.path
from optparse import OptionParser
import numpy as np
import wradlib as wrl
import h5py
import re
from cartopy.io import shapereader
import cartopy.crs as ccrs
from shapely.ops import cascaded_union
import raven_helper as rh

"""
	The module generates land mask for an input file (InputPathToHDF5File) based on the polygon from the shapefile (-s option). 
	The output HDF5 file will be written to output directory (OutputDirectory).
		input: [options] InputPathToHDF5File OutputDirectory
		output: HDF5 file with the (0/1) mask and the same number of datasets as the input file. 
	__________________________________________________________________
	Module requires some additional python libraries:
		wradlib, cartopy and shapely
	__________________________________________________________________
	Maryna Lukach, 2017 (RMI). 
"""

#-----------------------------------------------------------------------------------------
# Main function:

def main():
	
	usage = """%prog [options] InputPathToHDF5File OutputDirectory
    Use the option -h or --help to get all possible options
    """
	#determine basename of the application
	appbasename = os.path.basename(os.path.splitext(sys.argv[0])[0])

    #determine the basedir  of the application
	basedir = os.path.dirname(sys.argv[0]) + os.sep

    #the default logfile is located in the same directory as the program and
    #has the program name with the extension '.log'
	dfltlogfile = basedir + appbasename + '.log'
    
	# parsing options and arguments
	parser = OptionParser(usage=usage)
	parser.add_option("-d", "--debug", action="store_true", dest="debug", default=False, help="output debug messages during program execution")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="print messages (debug,info,warning,...) about the program execution to the console (stderr)")
	parser.add_option("-s", "--shape_file",  dest="shapeFile", default="/mnt/netapp/home/malukach/repos/vol2bird/src/vol2bird/data/bejab/bejab_poly_valid", help="The path and the name of the shapefiel")
	parser.add_option("-c", "--site_code",  dest="code", default="bejab", help="The site code for the output filename (default is 'bejab' and the filename will be bejab_landcoast_validmask.hdf).")
	
	(options, files) = parser.parse_args()
	debug = options.debug
	verbose = options.verbose
	shapeFile = options.shapeFile
	code = options.code
	if(verbose):print len(files)
	
	# get module input values
	if(len(files) < 1):DfltInputPath = "/mnt/netapp/home/malukach/repos/vol2bird/src/vol2bird/data/bejab/2016091423550500dBZ.pvol.h5" 
	else:DfltInputPath = files[0]
	if(verbose):print "input path is %s" %DfltInputPath
	
	if(len(files) < 2):DfltOutputPath = "/mnt/netapp/home/malukach/repos/vol2bird/src/vol2bird/data/bejab/" 
	else:DfltOutputPath = files[1]
	if(verbose):print "Output path is %s" %DfltOutputPath
	
	#create a logging object
	logger = rh.create_logger(appbasename, dfltlogfile, options)
	logger.info('Starting script %s ...' % sys.argv[0])
    
	inputFiles = files[:-1]
	outputDir = files[-1]

	if options.debug:
		logger.debug("inputfiles : {0}".format(inputFiles))
		logger.debug("outputDir : {0}".format(outputDir))
	
	logger.info('Opening the input file ...')
	f = rh.open_hdf5_file(DfltInputPath)
	logger.info('SUCCEEDED opening the input file ...')
	
	# get the number of datasets in the input file
	ntilt = rh.get_elevations_number(f)
	if(verbose):print "ntilt is %d" %ntilt
	(nelangle,nrays,nbins,rscale,rstart) = rh.get_grid_parameters(f,ntilt)
	sitecoords = rh.get_site_latlon(f)
	rh.close_hdf5_file(f)
	
	# get radar projection for the bins lat lon calculation
	proj_radar = wrl.georef.create_osr("aeqd", lat_0=sitecoords[1],lon_0=sitecoords[0])
	# get Earth radius at the latitude of the radar for the bins lat lon calculation
	radius = wrl.georef.get_earth_radius(sitecoords[1], proj_radar)
	
	# read the shapefile
	shp = shapereader.Reader(shapeFile)
	geoms = shp.geometries()
	polygon = cascaded_union(list(geoms))
	projection = ccrs.PlateCarree()
	
	# make neew output file
	fileExists = os.path.isfile("%s%s_landcoast_validmask.hdf" %(DfltOutputPath,code))
	if(verbose):print "%s%slandcoast_validmask.hdf file exists is %s" %(DfltOutputPath,code,fileExists)
	
	if (fileExists):
		fout = h5py.File("%s%s_landcoast_validmask.hdf" %(DfltOutputPath,code),"r+")
	else:
		fout = h5py.File("%s%s_landcoast_validmask.hdf" %(DfltOutputPath,code),"w")
	
	# create a list of 3D coordinate arrays where each entity corresponds to a given elevation 
	coord = []
	radlon = [] 
	radlat = []
	height = []
	# fill the list in with the coordinates
	for t in range(ntilt):
		# add empty arrays per elevation
		coord.append(np.empty((nrays[t], nbins[t], 3)))
		radlon.append(np.empty((nrays[t], nbins[t])))
		radlat.append(np.empty((nrays[t], nbins[t])))
		height.append(np.empty((nrays[t], nbins[t])))
		# get radar bins coordinates
		coord[t][...] = wrl.georef.sweep_centroids(nrays[t],rscale[t],nbins[t],nelangle[t])
		if(verbose):print coord[t].shape
		# calculate radar bins lat, lon and height positions
		radlon[t], radlat[t], height[t] = wrl.georef.polar2lonlatalt_n(coord[t][:,:,0],np.degrees(coord[t][:,:,1]),coord[t][:,:,2],sitecoords,re=radius,ke=4./3.)
		# make mask based on the polygon
		thismask = rh.inpolygon(polygon, radlon[t][:,:].ravel(), radlat[t][:,:].ravel())
		radmask = np.empty((nrays[t],nbins[t]), dtype=int)
		if(verbose):print thismask.size
		if(verbose):print (nrays[t],nbins[t])
		radmask[:,:] = np.asarray(thismask, dtype=int).reshape(nrays[t],nbins[t])
		# update or write down the mask
		if(fileExists):
			# if file already exists, update existing mask with current one
			arr = np.zeros((nrays[t],nbins[t]), dtype=int)
			arr = fout["/dataset%d/data1/data" %(t+1)]
			radmask[arr[...]==1] = 1
			fout["/dataset%d/data1/data" %(t+1)] = radmask
		else:
			# create new group and new dataset in the new file
			group = fout.create_group("/dataset%d/data1" %(e+1))
			dset = fout.create_dataset("/dataset%d/data1/data" %(e+1),data=radmask)
		
	fout.close()
	return

	
#-----------------------------------------------------------------------------------------
if __name__ == "__main__":
        main()



