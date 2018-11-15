#! /usr/bin/env python
"""
ohfa_helper.py  --  Helper module for the Odim Hdf5 file Aggregator

Aggregates two ore more Odim HDF5 files into a single one

Makes use of/and inspired from the threetoone.py script from Maryna Lukach

Christophe Ferauge , november 2014

Modifications for RAVen project by Maryna Lukach, 2017
added two functions for the files satisfying 1.2 version of ODIM HDF5 specification 
"""

import sys
import os
import socket
import re
import logging.handlers
import h5py
import numpy
import time
import calendar
import traceback
import numpy.ma as ma
import numpy as np
from pprint import pprint

# module wide variable which contains the name given to a dataset, this
# can be either 'dataset' or 'scan'
DATASETNAME = None
# module wide variable which contains the name given to an elevation angle
# of a dataset, this can be either 'angle' or 'elangle'
ANGLENAME = None
# module wide variables which contains the present ODIM  ormat or
# information model version
ODIM_CONVENTIONS = 'ODIM_H5/V2_2'
ODIM_VERSION = 'H5rad 2.2'
# module wide dictionnary which contains the ODIM source string for a
# specific radar
ODIM_SOURCE = {
    6451: 'WMO:06451,RAD:BX40,PLC:Zaventem,NOD:bezav,CTY:605,CMT:merged_volume_scan'}


def create_logger(appbasename, logfile, options, maxLogFileSize=1024 * 1024, maxLogFileRotation=5):
    """
    Creates a logger with a rotating log file of approx. 'maxLogFileSize'
    and a rotation of 'maxLogFileRotation'.

    Parameters
        appbasename : string representing the application basename to be used in the logging output
        logfile : string representing the full path of the log file
        options : object representing the application command line options
        maxLogFileSize : number representing the maximum size the logfile may grow before rotating
        maxLogFileRotation : number terpresenting the maximum number of logfile rotations

    returns the created logger object.
    """

    # set up a logger
    logger = logging.getLogger(appbasename)
    # sets the basic loglevel
    if options.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # logging output is always directed to a rotating logfile

    # format of the logging output to logfiles
    logfileformat = '%%(asctime)s %s %s[%%(process)s] : %%(levelname)s : %%(message)s' % (
        socket.gethostname(), appbasename)
    # create a formatter for logfile output
    logfile_formatter = logging.Formatter(logfileformat)
    # create a handler for rotating logfiles
    lfh = logging.handlers.RotatingFileHandler(
        logfile, maxBytes=maxLogFileSize, backupCount=maxLogFileRotation)
    # attach corresponding formatter to the logfile handler
    lfh.setFormatter(logfile_formatter)
    # add logfile handler to logger
    logger.addHandler(lfh)

    # if the verbos option is set, logging output is directed to the console,
    # having a simpler output format
    if options.verbose:
        # format of the logging output to the console
        consoleformat = '%s : %%(levelname)s : %%(message)s' % appbasename
        # create a formatter for the console output
        console_formatter = logging.Formatter(consoleformat)
        # create a handler for the console
        ch = logging.StreamHandler()
        # attach corresponding formatter to the console handler
        ch.setFormatter(console_formatter)
        # add logfile handler to logger
        logger.addHandler(ch)

    return logger


def terminate(logger=None, exitval=0):
    """
    Terminates the execution.

    Parameters
        logger : logging object
        exitval :  exit value , default 0

    Returns 0.
    """
    # write the info to the logfile
    if logger:
        logger.info('Terminating script %s (exitcode %d) ...\n' %
                    (sys.argv[0], exitval))
    # exit the execution code
    sys.exit(exitval)


def init_datasetname(fo, logger=None):
    """
    Determines the name of the hdf5 group used for the dataset

    Parameters
        fo : h5py file object
        logger : logging object

    Raises an exception when no valid name is found

    """
    global DATASETNAME
    datasetNames = ['dataset', 'scan']

    for name in datasetNames:
        if ("/" + name + "1") in fo:
            DATASETNAME = name
            return
    # Oops found no valid name
    errorMessage = "Could not find a valid name for the hdf5 group used for the dataset !!!"
    if logger:
        logger.error(errorMessage)
    raise Exception(errorMessage)


def init_anglename(fo, logger=None):
    """
    Determines the name of the hdf5 group used for the dataset

    Parameters
        fo : h5py file object
        logger : logging object

    Raises an exception when no valid name is found

    """
    global DATASETNAME
    global ANGLENAME
    angleNames = ['angle', 'elangle']

    for name in angleNames:
        if name in fo['/' + DATASETNAME + '1' + '/where'].attrs.keys():
            ANGLENAME = name
            return
    # Oops found no valid name
    errorMessage = "Could not find a valid name for the angle attribute in the dataset !!!"
    if logger:
        logger.error(errorMessage)
        raise Exception(errorMessage)


def open_hdf5_file(InputPath, logger=None):
    """
    Opens the hdf5 file in the reading mode. Writes info to the log file.

    Parameters
        InputPath : string representing the full path of the input file
        logger : logging object

    Returns the Id of the file object
    """
    # opening the hdf5 FileObject,
    if logger:
        logger.debug('Opening the hdf5 file {0} ...'.format(InputPath))

    # set FileObject to None
    FileObject = None

    try:
        # open the file in the reading mode
        FileObject = h5py.File(InputPath, 'r')
    except Exception as error:
        # write error to the logfile in case of exception
        if logger:
            logger.error('Could not open hdf5 file {0} ({1}) !'.format(
                InputPath, str(error)))
        raise error
    finally:
        # check if the file object is created
        if logger and FileObject is not None:
            # write to the logfile opening information
            logger.debug(
                'The hdf5 file {0} is successfully opened ...'.format(InputPath))

    init_datasetname(FileObject, logger)
    init_anglename(FileObject, logger)

    return (FileObject)


def close_hdf5_file(fo, logger=None):
    """
    Close given hdf5 file Object 'fo'.

    Parameters
        fo : h5py file object
        logger : logging object

    Raises an exception in case of problems
    """
    try:
        filename = fo.filename
        fo.close()
        if logger:
            logger.debug(
                'The hdf5 file {0} is successfully closed ...'.format(filename))
    except Exception as error:
        if logger:
            logger.error('Could not close hdf5 file {0} ({1}) !'.format(
                filename, str(error)))
        raise error


def getObjectType(fo, logger=None):
    """
    gets the Object type  of the ODIM hdf5 file.
    This function will try to extract the 'object' attribute from the
    top level 'what' group

      Parameters
          fo : h5py file object
          logger : logging object

      Returns the object type as a string
    """
    try:
        # return the object type (is a numpy string) and convert
        # it to a python string
        return str(fo['/what'].attrs['object'])
    except Exception as error:
        if logger:
            logger.error('Could not extract object type from ODIM hdf5 file {0} ({1}) !'.format(
                fo.filename, str(error)))
        raise error


def getStartTimestamp(fo, logger=None):
    """
    Gets the start timestamp of the ODIM hdf5 file
    timestamp from the first dataset,
    assuming it is the first generated
    the timestamp will be extracted from '/dataset1/what/startdate' and
    '/dataset1/what/starttime'
    Raises an exception if the timestamp extraction fails
    Possibly dataset1 can be another group name,eg. 'scan1', this is
    covered by the global variable 'DATASETNAME'

    Parameters
        fo : h5py file object
        logger : logging object

    Returns the starttimestamp as an epoch
    """
    global DATASETNAME

    try:
        what = '/' + DATASETNAME + '1/what'
        timestampString = fo[what].attrs['startdate'] + \
            fo[what].attrs['starttime'] + 'UTC'
        startTimestamp = time.strptime(timestampString, "%Y%m%d%H%M%S%Z")
        if logger:
            logger.debug('{0} has timestamp {1} .'.format(
                fo.filename, time.strftime("%Y-%m-%d %H:%M:%S", startTimestamp)))
        return calendar.timegm(startTimestamp)
    except Exception as error:
        if logger:
            logger.error('Could not extract timestamp from ODIM hdf5 file {0} ({1}) !'.format(
                fo.filename, str(error)))
        raise error


def getElevations(fo, logger=None):
    """
    Gets the elevations of a Polar volume (PVOL)

    Parameters
        fo : h5py file object
        logger : logging object

    Returns the elevations as a dictionary eg. {'dataset1' : 0.1, 'dataset2' : 0.3, ..., 'datasetN' : 25.0)
    """
    global DATASETNAME
    global ANGLENAME

    elevations = {}
    try:
        # first check if we really have to do with a polar volume

        # 1 extract elevations now
        # 1.1 get a list of all datasets
        datasets = [group for group in list(
            fo['/']) if group.startswith(DATASETNAME)]
        # 1.2 itterate through all datasets and extract elevations
        for ds in datasets:
            elevations[ds] = fo['/' + ds + '/where'] .attrs[ANGLENAME]
            if logger:
                logger.debug(
                    '{0} - {1} has elevation angle: {2}'.format(fo.filename, ds, elevations[ds]))
    except Exception as error:
        if logger:
            logger.error('Could not extract elevations from ODIM hdf5 file {0} ({1}) !'.format(
                fo.filename, str(error)))
        raise error

    return elevations

def getMoment(fo, logger=None):
    """
    Gets the moment of a Polar volume (PVOL)

    Parameters
        fo : h5py file object
        logger : logging object

    Returns the elevations as a dictionary eg. {'dataset1' : 0.1, 'dataset2' : 0.3, ..., 'datasetN' : 25.0)
    """
    global DATASETNAME
    global ANGLENAME
    #print "in getMoments"
    #print fo
    moments = {}
    try:
        # first check if we really have to do with a polar volume

        # 1 extract elevations now
        # 1.1 get a list of all datasets
        datasets = [group for group in list(
            fo['/']) if group.startswith(DATASETNAME)]
        # 1.2 itterate through all datasets and extract elevations
        for ds in datasets:
            moments[ds] = fo['/' + ds + '/data1/what'] .attrs['quantity']
            #print moments[ds]
            if logger:
                logger.debug(
                    '{0} - {1} has moment: {2}'.format(fo.filename, ds, elevations[ds]))
    except Exception as error:
        if logger:
            logger.error('Could not extract moment from ODIM hdf5 file {0} ({1}) !'.format(
                fo.filename, str(error)))
        raise error

    return moments

def isExcludedElevation(elangle, excludeAngles):
    """
      Determines if an elevation is to be excluded according to a specific list
      Parameter
          elangle : floating point number representing the elevation to check for exclusion
          excludeAngles : list of elevation angles to be excluded

      Returns
          True when the elevation should be excluded, False otherwise
    """
    for ea in excludeAngles:
        if numpy.allclose(elangle, ea):
            return True
    return False


def getDatasetNumber(dataset):
    """
    Helper function which gives the number of a given dataset

    Parameters
      dataset : string representing the dataset from which you want to retrieve the number
                eg. 'dataset10'

    Returns the number of the dataset as an integer
    """
    return int(dataset.lstrip('dataset'))


def getNewDatasetName(datasets):
    """
    Helper function which gives the name of the next (new) dataset

      Parameters:
        datasets : an ascending ordered string array representing the already existing datasets
                    eg. ['dataset1', 'dataset2', 'dataset3', 'dataset4', 'dataset5']

      Returns the name of the new dataset eg. 'dataset6'
    """
    if(len(datasets)):
        return 'dataset' + str(getDatasetNumber(datasets[-1]) + 1)
    else:
        return 'dataset1'


def getNormalizedTimeString(timeString, period=0):
    """
        Helper function which normalizes a timeString (eg. '153723') according to a certain
        period (eg. 5 minutes): eg. '153723' becomes '153500'
        if period <= 0, does nothing

            Parameters:
              timeString : the incoming timeString to be normalized
              period : normalization period in minutes, default 0

            Returns the normalized timeString
    """
    if( period <= 0):
        return timeString
    else:
        hour,minute = [timeString[0:2], timeString[2:4]]
        normalizedMinute = (int(minute)/period)*period
        return hour + '{:02d}00'.format(normalizedMinute)


def initHelperDict(cd={}):
    """
    Initializes a helper dictionary which will contains all data necessary to build up
    the compound files

      Parameters:

      Returns an initialized helper dict object
    """
    cd = {'common': {'group': {}, 'attributes': {}}, 'moments': {}, 'aux': {}}
    return cd


def constructHelperDictCommon(fos, cd):
    """
    Function which builts up a helper dictionary used for creating the compound files.
    This function contributes to that part which is common for all datasets

    Parameters
        fos : h5py file object array, containing the ODIM hdf5 file objects to be aggregated
        cd : helper dictionary which contains all data necessary to build up the compound files
    """
    # get the top level where if it exists
    if '/where' in fos[0]:
        cd['common']['group']['where'] = (fos[0])['/where']
    # get the top level what if it exists
    if '/what' in fos[0]:
        cd['common']['group']['what'] = (fos[0])['/what']
    # get the top level how if it exists
    if '/how' in fos[0]:
        cd['common']['group']['how'] = (fos[0])['/how']
        # get 'endepochs' from /how in the chronologic last file
        if 'endepochs' in (fos[-1])['/how'].attrs:
            cd['aux']['endepochs'] = (fos[-1])['/how'].attrs['endepochs']
    # get the top level attributes
    cd['common']['attributes'] = (fos[0]).attrs
    # pprint(cd)

def constructHelperDictData(fos, cd, excludeAngles=[]):
    """
    Function which builts a helper dictionary used for creating the compound files.
    This function contributes to the Dataset part of the dictionary

    Parameters
        fos : h5py file object array, containing the ODIM hdf5 file objects to be aggregated
        cd : helper dictionary which contains all data necessary to build up the compound files
        excludeAngles : list of angles to be excluded from the aggregation

    """
    # dictionary to convert odim HDF5 quantity to well known radar moment
    # abbreviations
    moments = {'DBZH': 'dBZ', 'VRAD': 'V', 'TH': 'dBuZ', 'KDP': 'KDP', 'ZDR': 'ZDR', 'PHIDP': 'PhiDP', 'RHOHV': 'RhoHV', 'WRAD':'W'}

    # itterate through all hdf5 file objects to build up the helper dictionary
    for fo in fos:
        # print '\n' + fo.filename
        # itterate through all datasets
        datasets = sorted([ds for ds in list(fo['/'])
                           if ds.startswith('dataset')], key=getDatasetNumber)
        for ds in datasets:
            # print '\t--> {0}'.format(ds)
            # exclude angles if requested
            if isExcludedElevation(fo['/' + ds + '/where'].attrs['elangle'], excludeAngles):
                # skip this dataset
                continue
            # itterate through all data within a dataset
            data = [da for da in list(fo['/' + ds]) if da.startswith('data')]
            for da in data:
                quantity = fo['/' + ds + '/' + da + '/what'].attrs['quantity']
                moment = moments[quantity]
                if moment in cd['moments'].keys():
                    existingDatasets = sorted(
                        cd['moments'][moment].keys(), key=getDatasetNumber)
                else:
                    cd['moments'][moment] = {}
                    existingDatasets = []
                newDatasetName = getNewDatasetName(existingDatasets)
                cd['moments'][moment][newDatasetName] = {}
                cd['moments'][moment][newDatasetName][
                    'data1'] = fo['/' + ds + '/' + da]
                cd['moments'][moment][newDatasetName][
                    'what'] = fo['/' + ds + '/what']
                cd['moments'][moment][newDatasetName][
                    'where'] = fo['/' + ds + '/where']
              


def constructCompoundFiles(cd, base, period=0, logger=None):
    """
    Construct the compound files now from the helper dictionnary

    Parameters
        cd : helper dictionary which contains all data necessary to build up the compound files
        base :  string representing the basename for each compound filename
        period : Normalization period in minutes. If period is greater than zero the time in the
                 toplevel /what/time string will be normailzed according to 'period'. For example
                 if /what/time = '172156' and period is 5, /what/time will be changed to '172000'
        logger : logging object
    """
    # itterate through all the moments
    for moment in cd['moments'].keys():
        outputPath = base + '.' + moment + '.hdf'
        if logger:
            logger.info('Creating output file {0} for moment {1} ...'.format(
                outputPath, moment))
        fo = h5py.File(outputPath,'w')# 'a')
        # construct the common /where, /what and /how
        for commonGroup in ['where', 'what', 'how']:
            if commonGroup in cd['common']['group']:
                group = fo.create_group(commonGroup)
                for attrName, attrVal in cd['common']['group'][commonGroup].attrs.items():
                    group.attrs[attrName] = attrVal
                # Now treat the special case for the 'endepochs' attribute in
                # the commongroup 'how'
                if commonGroup == 'how' and ('endepochs' in group.attrs):
                    group.attrs['endepochs'] = cd['aux']['endepochs']
                # Now treat the special case where you want to normalize /what/time
                # in the top level group
                if commonGroup == 'what' and ('time' in cd['common']['group']['what'].attrs.keys()) and (period >= 0):
                    whatTime = cd['common']['group']['what'].attrs['time']
                    normalizedWhatTime = getNormalizedTimeString(whatTime, period)
                    group.attrs['time'] = normalizedWhatTime
                    logger.info('Normalized the /what/time of "{0}" to "{1}" ...'.format(whatTime, normalizedWhatTime))

        # copy the common attributes
        for attrName, attrVal in cd['common']['attributes'].items():
            fo['/'].attrs[attrName] = attrVal
        # copy the datasets for the specific moment
        for dataset in cd['moments'][moment].keys():
            # print "Creating group {0} ... ".format(dataset)
            fo.create_group(dataset)
            fo.copy(cd['moments'][moment][dataset]['data1'], fo['/' + dataset])
            # rename every data group to data1 (if it has not that name
            # already) ... as there is only one quantity per dataset
            dataName = (cd['moments'][moment][dataset]
                        ['data1']).name.split('/')[-1]
            if dataName != 'data1':
                fo['/' + dataset + '/data1'] = fo['/' + dataset + '/' + dataName]
                del fo['/' + dataset + '/' + dataName]
            fo.copy(cd['moments'][moment][dataset]['what'], fo['/' + dataset])
            fo.copy(cd['moments'][moment][dataset]['where'], fo['/' + dataset])
        fo.close()



def constructOneCompoundFile(cd, base, period=0, logger=None, OffshoreFile=None, offShore=False):
    """
    Construct the compound files now from the helper dictionnary

    Parameters
        cd : helper dictionary which contains all data necessary to build up the compound files
        base :  string representing the basename for each compound filename
        period : Normalization period in minutes. If period is greater than zero the time in the
                 toplevel /what/time string will be normailzed according to 'period'. For example
                 if /what/time = '172156' and period is 5, /what/time will be changed to '172000'
        logger : logging object
    """
    if(OffshoreFile is None):
		outputPath = base + '.h5'
    else:
		if(offShore):
			outputPath = base + '.offshore.h5'
		else:
			outputPath = base + '.onshore.h5'
    if logger:
        logger.info('Creating output file {0} ...'.format(outputPath))
	# print 'Creating output file {0} for moment {1}
	# ...'.format(outputPath, moment)
    fo = h5py.File(outputPath, 'w')#'a')
    # index of the /data group inside each /dataset group (depends on the moment) 
    dataIndex = 1
    # Used to calculate the number of the radar bins in ranges interval from 5 km to 35 km and below 4 000 km height. 
    #H = np.arange(200,4200,200)

	# itterate through all the moments
    for moment in cd['moments'].keys():
        if(dataIndex == 1):
	        # construct the common /where, /what and /how
	        for commonGroup in ['where', 'what', 'how']:
	            if commonGroup in cd['common']['group']:
	                group = fo.create_group(commonGroup)
	                for attrName, attrVal in cd['common']['group'][commonGroup].attrs.items():
	                    group.attrs[attrName] = attrVal
	                # Now treat the special case for the 'endepochs' attribute in
	                # the commongroup 'how'
	                if commonGroup == 'how' and ('endepochs' in group.attrs):
	                    group.attrs['endepochs'] = cd['aux']['endepochs']
	                # Now treat the special case where you want to normalize /what/time
	                # in the top level group
	                if commonGroup == 'what' and ('time' in cd['common']['group']['what'].attrs.keys()) and (period >= 0):
	                    whatTime = cd['common']['group']['what'].attrs['time']
	                    normalizedWhatTime = getNormalizedTimeString(whatTime, period)
	                    group.attrs['time'] = normalizedWhatTime
	                    logger.info('Normalized the /what/time of "{0}" to "{1}" ...'.format(whatTime, normalizedWhatTime))
	        # copy the common attributes
	        for attrName, attrVal in cd['common']['attributes'].items():
	            fo['/'].attrs[attrName] = attrVal
        # copy the datasets for the specific moment
        for dataset in cd['moments'][moment].keys():
            if(dataIndex == 1):
				fo.create_group(dataset)
            X1 = numpy.asarray(cd['moments'][moment][dataset]['data1']['data'])
            if(OffshoreFile is not None):
	            X2 = numpy.asarray(OffshoreFile[dataset + "/data1/data"])
	            # data outside the needed part are set to nodata of ODIM HDF5 format
	            if (offShore):
					X1[X2!=0]=255
	            else:
					X1[X2==0]=255
            dset = fo.create_dataset('/' + dataset + '/data%i/data'  %dataIndex, data=X1)
            fo.copy(cd['moments'][moment][dataset]['data1']['what'], fo['/' + dataset], name='/' + dataset + '/data%i/what'  %dataIndex)
            if(dataIndex == 1):fo.copy(cd['moments'][moment][dataset]['what'], fo['/' + dataset])
            if(dataIndex == 1):fo.copy(cd['moments'][moment][dataset]['where'], fo['/' + dataset])
			# Used to calculate the number of the radar bins in ranges interval from 5 km to 35 km and below 4 000 km height. 
            #elv = fo['/' + dataset + '/where'].attrs['elangle']
            #print "elevation is"
            #print elv 
            #sinelv = np.sin(np.deg2rad(elv))
            #satranges = (H/sinelv)
            #indxr = satranges[-1]*1./fo['/' + dataset + '/where'].attrs['rscale'] 
            #print "total number of points"
            #print len(np.rint(X1[:,10:min(np.rint(indxr),71)]).flatten())
            #print "number of a part"
            #print len(np.rint(X1[X1[:,10:min(np.rint(indxr),71)]==255]).flatten())
            #print "part is offShore"
            #print offShore
        dataIndex += 1
            
    fo.close()

def aggregate(fos, base, logger=None, excludeAngles=[], period=0, oneFile=False, OffshoreFile=None, offShore=False):
    """
    Aggregates a number of ODIM hdf5 file objects into files of single data type
    Each compound file has a filename composed of the argument "base" + "data type string" (eg. dBZ, V, W, dBuZ, ...)
    + ".hdf"

    Parameters
        fos : h5py file object array, containing the ODIM hdf5 file objects to be aggregated
        base :  string representing the basename for each compound filename
        logger : logging object
        excludeAngles : list of angles to be excluded from the aggregation
        period : Normalization period in minutes. If period is greater than zero the time in the
                 toplevel /what/time string will be normailzed according to 'period'. For example
                 if /what/time = '172156' and period is 5, /what/time will be changed to '172000'
    """
    global DATASETNAME
   
    try:
            # create a helper dictionary which contains all data necessary to
            # build up the compound files
        cd = initHelperDict()
        
        if DATASETNAME == 'scan':
                # Build up the common part of the helper dictionary ...
            constructHelperDictCommonForODIM_1_2(fos, cd)
            # Build up the Data part of the helper dictionary ...
            constructHelperDictDataForODIM_1_2(fos, cd, excludeAngles)
            # Now construct the compound files ...
            constructCompoundFilesForODIM_1_2(cd, base, period, logger)
        else:
            # build up the common part of the helper dictionary
            constructHelperDictCommon(fos, cd)
            # build up the Data part of the helper dictionary
            constructHelperDictData(fos, cd, excludeAngles)
            # Now construct the compound files
            if(oneFile):
			 	constructOneCompoundFile(cd, base, period, logger, OffshoreFile, offShore=offShore)
            else:
            	constructCompoundFiles(cd, base, period, logger)
    except Exception as error:
        if logger:
            logger.error(
                'Oops ... something went wrong while aggregating the odim hdf5 files ({0})...'.format(str(error)))
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


def constructHelperDictCommonForODIM_1_2(fos, cd):
    """
    Function which builts up a helper dictionary used for creating the compound files,
    when the input files are ODIM version 1.2.
    This function contributes to that part which is common for all datasets

    Parameters
        fos : h5py file object array, containing the ODIM hdf5 file objects (version 1.2)  to be aggregated
        cd : helper dictionary which contains all data necessary to build up the compound files

    """
    global ODIM_VERSION
    global ODIM_SOURCE
    # process the top level where if it exists
    if '/where' in fos[0]:
        where = {}
        for attrName in ['lon', 'lat', 'height']:
            if attrName in (fos[0])['/where'].attrs:
                where[attrName] = (fos[0])['/where'].attrs[attrName]
        cd['common']['group']['where'] = where
    # process the top level what if it exists
    if '/what' in fos[0]:
        what = {}
        for attrName in ['object', 'date', 'time']:
            if attrName in (fos[0])['/what'].attrs:
                what[attrName] = (fos[0])['/what'].attrs[attrName]
        what['version'] = numpy.string_(ODIM_VERSION)
        if ('WMO' in (fos[0])['/how'].attrs) and ((fos[0])['/how'].attrs['WMO'] in ODIM_SOURCE):
            what['source'] = numpy.string_(
                ODIM_SOURCE[(fos[0])['/how'].attrs['WMO']])
        cd['common']['group']['what'] = what
    # process the top level how if it exists
    if '/how' in fos[0]:
        how = {}
        for attrName in ['beamwidth', 'NEZH', 'radconstH', 'radconstV', 'NI', 'startepochs',  'system', 'software', 'wavelength', 'azmethod']:
            if attrName in (fos[0])['/how'].attrs:
                how[attrName] = (fos[0])['/how'].attrs[attrName]
        # get 'endepochs' from /how in the chronologic last file
        if 'endepochs' in (fos[-1])['/how'].attrs:
            how['endepochs'] = (fos[-1])['/how'].attrs['endepochs']
        cd['common']['group']['how'] = how
    # set the top level attributes
    cd['common']['attributes']['Conventions'] = numpy.string_(ODIM_CONVENTIONS)
    # pprint(cd)


def constructHelperDictDataForODIM_1_2(fos, cd, excludeAngles=[]):
    """
    Function which builts a helper dictionary used for creating the compound files,
    when the input files are ODIM version 1.2.
    This function contributes to the Dataset part of the dictionary

    Parameters
        fos : h5py file object array, containing the ODIM hdf5 file objects (version 1.2) to be aggregated
        cd : helper dictionary which contains all data necessary to build up the compound files
        excludeAngles : list of angles to be excluded from the aggregation

    """
    # dictionary to convert odim HDF5 quantity to well known radar moment abbreviations
    moments = {'DBZH': 'dBZ', 'VRAD': 'V', 'TH': 'dBuZ', 'KDP': 'KDP', 'ZDR': 'ZDR', 'PHIDP': 'PhiDP', 'RHOHV': 'RhoHV', 'WRAD':'W'}


    # itterate through all hdf5 file objects to build up the helper dictionary
    for fo in fos:
        # itterate through all scans
        scans = sorted([sc for sc in list(
            fo['/']) if sc.startswith('scan')], key=lambda s: int(s.lstrip('scan')))
        for sc in scans:
            # exclude angles if requested
            if isExcludedElevation(fo['/' + sc + '/where'].attrs['angle'], excludeAngles):
                # skip this dataset
                continue
            quantity = fo['/' + sc + '/what'].attrs['quantity']
            moment = moments[quantity]
            if moment in cd['moments'].keys():
                existingDatasets = sorted(
                    cd['moments'][moment].keys(), key=getDatasetNumber)
            else:
                cd['moments'][moment] = {}
                existingDatasets = []
            newDatasetName = getNewDatasetName(existingDatasets)
            cd['moments'][moment][newDatasetName] = {}
            cd['moments'][moment][newDatasetName]['data1'] = {}
            cd['moments'][moment][newDatasetName][
                'data1']['data'] = fo['/' + sc + '/data']
            cd['moments'][moment][newDatasetName]['data1']['what'] = {}
            #process /datasetN/data1/what
            for attrName in ['quantity', 'gain', 'offset', 'nodata', 'undetect']:
                if attrName in fo['/' + sc + '/what'].attrs:
                    if attrName == 'quantity' and fo['/' + sc + '/what'].attrs['quantity'] == 'DBZ':
                        cd['moments'][moment][newDatasetName]['data1'][
                            'what']['quantity'] = numpy.string_('DBZH')
                    else:
                        cd['moments'][moment][newDatasetName]['data1']['what'][
                            attrName] = fo['/' + sc + '/what'].attrs[attrName]
            #process /datasetN/what
            cd['moments'][moment][newDatasetName]['what'] = {}
            for attrName in ['product', 'startdate', 'starttime', 'enddate', 'endtime']:
                if attrName in fo['/' + sc + '/what'].attrs:
                    cd['moments'][moment][newDatasetName]['what'][
                        attrName] = fo['/' + sc + '/what'].attrs[attrName]
            #process /datasetN/where
            cd['moments'][moment][newDatasetName]['where'] = {}
            cd['moments'][moment][newDatasetName]['where'][
                'elangle'] = fo['/' + sc + '/where'].attrs['angle']
            cd['moments'][moment][newDatasetName][
                'where']['a1gate'] = numpy.int_(0)
            cd['moments'][moment][newDatasetName]['where'][
                'nbins'] = fo['/where'].attrs['xsize']
            cd['moments'][moment][newDatasetName][
                'where']['rstart'] = numpy.double(0.0)
            cd['moments'][moment][newDatasetName]['where'][
                'rscale'] = fo['/where'].attrs['xscale']
            cd['moments'][moment][newDatasetName]['where'][
                'nrays'] = fo['/where'].attrs['ysize']
            #process /datasetN/how
            if ('/' + sc + '/how') in fo:
                cd['moments'][moment][newDatasetName]['how'] = {}
                for attrName in ['startepochs', 'endepochs']:
                    if attrName in fo['/' + sc + '/how'].attrs:
                        cd['moments'][moment][newDatasetName]['how'][
                            attrName] = fo['/' + sc + '/how'].attrs[attrName]
                # make following attributes specific to the dataset as they can
                # differ from file to file
                for attrName in ['highprf', 'lowprf', 'rpm', 'pulsewidth']:
                    if attrName in fo['/how'].attrs:
                        cd['moments'][moment][newDatasetName]['how'][
                            attrName] = fo['/how'].attrs[attrName]
                            



def constructCompoundFilesForODIM_1_2(cd, base, period=0, logger=None):
    """
    Construct the compound files now from the helper dictionnary

    Parameters
        cd : helper dictionary which contains all data necessary to build up the compound files
        base :  string representing the basename for each compound filename
        period : Normalization period in minutes. If period is greater than zero the time in the
                 toplevel /what/time string will be normailzed according to 'period'. For example
                 if /what/time = '172156' and period is 5, /what/time will be changed to '172000'
        logger : logging object
    """
    # itterate through all the moments
    for moment in cd['moments'].keys():
        outputPath = base + '.' + moment + '.hdf'
        if logger:
            logger.info('Creating output file {0} for moment {1} ...'.format(
                outputPath, moment))
        fo = h5py.File(outputPath, 'w')
        # construct the common /where, /what and /how
        for commonGroup in ['where', 'what', 'how']:
            if commonGroup in cd['common']['group']:
                group = fo.create_group(commonGroup)
                for attrName, attrVal in cd['common']['group'][commonGroup].items():
                    group.attrs[attrName] = attrVal
                # Now treat the special case where you want to normalize /what/time
                # in the top level group
                if commonGroup == 'what' and ('time' in cd['common']['group']['what'].attrs.keys()) and (period >= 0):
                    whatTime = cd['common']['group']['what'].attrs['time']
                    normalizedWhatTime = getNormalizedTimeString(whatTime, period)
                    group.attrs['time'] = normalizedWhatTime
                    logger.info('Normalized the /what/time of "{0}" to "{1}" ...', whatTime, normalizedWhatTime)
        # construct the common attributes
        for attrName, attrVal in cd['common']['attributes'].items():
            fo['/'].attrs[attrName] = attrVal
        # construct the datasets for the specific moment
        for dataset in cd['moments'][moment].keys():
            # print "Creating group {0} ... ".format(dataset)
            # construct the dataset /datasetN/data1/data
            fo.create_group(dataset)
            data1Group = fo.create_group('/' + dataset + '/data1')
            datasetData = cd['moments'][moment][dataset]['data1']['data']
            fo['/' + dataset + '/data1'].create_dataset(
                "data", data=datasetData, chunks=datasetData.shape, compression=6)
            fo['/' + dataset +
                '/data1/data'].attrs['CLASS'] = numpy.string_('IMAGE')
            fo['/' + dataset +
                '/data1/data'].attrs['IMAGE_VERSION'] = numpy.string_('1.2')
            # construct the group /datasetN/data1/what
            data1whatGroup = data1Group.create_group('what')
            for attrName, attrVal in cd['moments'][moment][dataset]['data1']['what'].items():
                data1whatGroup.attrs[attrName] = attrVal
            # construct the groups /where, /what and /how in /dataset1
            for datasetGroup in ['where', 'what', 'how']:
                if datasetGroup in cd['moments'][moment][dataset] and (len(cd['moments'][moment][dataset][datasetGroup]) > 0):
                    dataGroup = fo.create_group(
                        '/' + dataset + '/' + datasetGroup)
                    for attrName, attrVal in cd['moments'][moment][dataset][datasetGroup].items():
                        dataGroup.attrs[attrName] = attrVal
        fo.close()
