#! /usr/bin/env python
"""
ohfa.py  --  ohfa.py  --  Odim Hdf5 file Aggregator

Aggregates two or more Odim HDF5 files into a single one

Makes use of/and inspired from the threetoone.py script from Maryna Lukach

Christophe Ferauge , november 2014

Modifications for the RAVen project by Maryna Lukach, 2017
-o and -c options are added and masking of on- or off- shore parts of the data is done based on the mask provided in -c option
"""


import sys
import os
import os.path
import socket
import re
from optparse import OptionParser
import time
import ohfa_helper as oh
import h5py

def main():
    usage = """%prog [options] inputfile1, inputfile2, ..., inputfileN outputdir

    Use the option -h or --help to get all possible options
    """
    #determine basename of the application
    appbasename = os.path.basename(os.path.splitext(sys.argv[0])[0])

    #determine the basedir  of the application
    basedir = os.path.dirname(sys.argv[0]) + os.sep

    #the default logfile is located in the same directory as the program and
    #has the program name with the extension '.log'
    dfltlogfile = basedir + appbasename + '.log'

    #parsing options and arguments ./testdata/bejab_landcoast_validmask.hdf
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="print messages (debug,info,warning,...) about the program execution to the console (stderr)")
    parser.add_option("-d", "--debug", action="store_true", dest="debug", default=False, help="output debug messages during program execution")
    parser.add_option("-e", "--exclude", type="float", action="append", dest="excludeAngles", default=[], help="Exclude this angle from the aggregation. The angle must be given as a float. By repeating this option you can give multiple angles to exclude.")
    parser.add_option("-n", "--namebase", type="string", action="store", dest="nameBase", help="Part of the output filename put in front of the name.")
    parser.add_option("-t", "--normalizeTime", type="int", action="store", dest="normalizeTimePeriod", help="Normalize the time of the aggregated file to a specific period in minutes. This will only change the timestamp in /what/time in the top level group and not altering other time attributes, neither the filename")
    parser.add_option("-c", "--coastline", type="string", action="store", dest="OffShoreFile", default=None, help="coastline data filename for the offshore/ onshore data splitcing. If None - no data splitcing is done.")
    parser.add_option("-o", "--off-shore", action="store_true", dest="offShore", default=False, help="off-shore option (True for offshore data onle, False for onshore data only)")
    
    (options, files) = parser.parse_args()

    if len(files) < 3 :
        if options.verbose or options.debug:
            print "\nOops, at least two input files and one output path should be given as arguments, terminating ...\n"
            parser.print_usage()
        sys.exit(1)

    inputFiles = files[:-1]
    outputDir = files[-1]

    #create a logging object
    logger = oh.create_logger(appbasename, dfltlogfile, options)
    logger.info('Starting script %s ...' % sys.argv[0])

    if options.debug:
      logger.debug("inputfiles : {0}".format(inputFiles))
      logger.debug("outputDir : {0}".format(outputDir))

    #try to open the input files and put file objects in an array
    #If a failure occurs, terminate the program
    inpFilesObjs = [];
    for inpFile in inputFiles:
      try:
        inpFilesObjs.append( oh.open_hdf5_file(inpFile,logger) )
        if oh.getObjectType(inpFilesObjs[-1]) != 'PVOL':
          logger.error('File {0} does not represent a polar volume !'.format(inpFilesObjs[-1].filename))
          raise
      except:
        logger.error('FAILED opening all {0} input files ...\n'.format(len(inputFiles)))
        oh.terminate(logger, 1)
    logger.info('SUCCEEDED opening all {0} input files ...'.format(len(inputFiles)))
    
    if(options.OffShoreFile is not None):inputOffshoreFile = h5py.File(options.OffShoreFile,"r") 
    else: inputOffshoreFile = None
    
    offShore = options.offShore

    #sort files on their start timestamp
    logger.info('Sorting input files based on their start timestamp ...')
    inpFilesObjs.sort(key=oh.getStartTimestamp)
    firstTimeStamp = oh.getStartTimestamp(inpFilesObjs[0])
    lastTimeStamp = oh.getStartTimestamp(inpFilesObjs[-1])
    
    if (lastTimeStamp - firstTimeStamp==0):
        if options.debug:
	      #print "Same time stemp"
	      for ifo in inpFilesObjs:
	        logger.debug('File {0} - elevations {1}'.format(ifo.filename, oh.getElevations(ifo, logger)))
	        logger.debug('File {0} has moment {1} .'.format(ifo.filename, oh.getMoment(ifo)))

	    #create a basename for an output filename(s)
        if options.nameBase:
          outputBase = os.path.join(outputDir,options.nameBase)
        else:
          outputBase = os.path.join(outputDir, os.path.basename(os.path.splitext(inpFilesObjs[0].filename)[0].split(".vol")[0] + '.pvol'))    
        
        if options.debug:
	      logger.debug('Basename for the output filenames : {0}'.format(outputBase))
	    
		#do the aggregation
        oh.aggregate(inpFilesObjs, outputBase, logger, excludeAngles = options.excludeAngles, period = options.normalizeTimePeriod, oneFile=True, OffshoreFile=inputOffshoreFile, offShore=offShore)        
    else:
        
	    if options.debug:
	      for ifo in inpFilesObjs:
	        logger.debug('File {0} - elevations {1}'.format(ifo.filename, oh.getElevations(ifo, logger)))
	        logger.debug('File {0} has timestamp {1} .'.format(ifo.filename, time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(oh.getStartTimestamp(ifo)))))
	    
	    #create a basename for an output filename(s)
	    if options.nameBase:
	      outputBase = os.path.join(outputDir,options.nameBase)
	    else:
	      outputBase = os.path.join(outputDir, os.path.basename(os.path.splitext(inpFilesObjs[0].filename)[0]))
	
	    if options.debug:
	      logger.debug('Basename for the output filenames : {0}'.format(outputBase))
	
	    #do the aggregation
	    oh.aggregate(inpFilesObjs, outputBase, logger, excludeAngles = options.excludeAngles, period = options.normalizeTimePeriod)
	
    #close all input files
    for ifo in inpFilesObjs:
      try:
        oh.close_hdf5_file(ifo, logger)
      except Exception as error:
        logger.error('FAILED closing all {0} input files ({1}) ...\n'.format(len(inputFiles), str(error)))
        oh.terminate(logger, 1)

    logger.info('SUCCEEDED closing all {0} input files ...'.format(len(inputFiles)))
    oh.terminate(logger, 0)


if __name__ == "__main__":
        main()
