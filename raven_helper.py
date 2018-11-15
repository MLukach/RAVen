import numpy as np
import logging.handlers
import socket
import h5py
import re
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.collections import PatchCollection
from mpl_toolkits.basemap import Basemap
	
"""
	Library of the Python modules programmed for the RAVen project.
	__________________________________________________________________
	Module requires some additional python libraries:
		wradlib, cartopy and shapely
	__________________________________________________________________
	Maryna Lukach, 2017 (RMI). 
"""

#-----------------------------------------------------------------------------------------
# General functions:

def create_logger(appbasename, logfile, options, maxLogFileSize = 1024*1024, maxLogFileRotation = 5):
    """
    Creates a logger with a rotating log file of approx. 'maxLogFileSize' 
    and a rotation of 'maxLogFileRotation'.
    returns the created logger object.
    """
    
    #set up a logger
    logger = logging.getLogger(appbasename)
    #set the basic loglevel
    if options.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)    
    
    #logging output is always directed to a rotating logfile
    
    #format of the logging output to logfiles
    logfileformat = '%%(asctime)s %s %s[%%(process)s] : %%(levelname)s : %%(message)s' % (socket.gethostname(), appbasename)
    #create a formatter for logfile output
    logfile_formatter = logging.Formatter(logfileformat)
    #create a handler for rotating logfiles
    lfh = logging.handlers.RotatingFileHandler(logfile, maxBytes = maxLogFileSize, backupCount = maxLogFileRotation)
    #attach corresponding formatter to the logfile handler
    lfh.setFormatter(logfile_formatter)
    #add logfile handler to logger
    logger.addHandler(lfh)
    
    #if the verbos option is set, logging output is directed to the console,
    #having a simpler output format
    if options.verbose:    
        #format of the logging output to the console
        consoleformat = '%s : %%(levelname)s : %%(message)s' % appbasename
        #create a formatter for the console output
        console_formatter = logging.Formatter(consoleformat)
        #create a handler for the console
        ch = logging.StreamHandler()
        #attach corresponding formatter to the console handler
        ch.setFormatter(console_formatter)
        #add logfile handler to logger
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

#-----------------------------------------------------------------------------------------
# Open and close HDF5 file and read dataset values functions:
	
def open_hdf5_file(InputPath, logger=None):
    """
    Opens the hdf5 file in the reading mode. Writes info to the log file. (from ohfa_helper.py)

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

def read_dataset(DfltInputPath, DatasetName = '/dataset1/quality'):
	"""
	Read dataset DatasetName from DfltInputPath file and return it as np.array 
	"""
	h5File = h5py.File(DfltInputPath, "r")
	print DfltInputPath
	print DatasetName
	quality = h5File[DatasetName]
	print "quality"
	print quality
	
	data = np.zeros(quality.shape,dtype = quality.dtype)
	
	quality.read_direct(data)
	
	h5File.close() 
	return(data)

def decode_value(data,gain,offset,nodata=None,noval=np.nan):
        """Decoded data into real value using a gain offset and nodata.
           By default NODATA is set to NAN values for convenience.        
        """
        inodata = (data == nodata)
        print "in red data inodata is"
        print len(inodata.flatten())
        data = offset + data * gain        
        if nodata is not None:
            data[inodata] = noval
        return(data)




def get_grid_parameters(f,nElevations):
	"""
	Gets parameters of the grid for each elevation
		input: 	f - file object of h5py library;
				nElevations - number of elevations in the file;
		output: (nElAngle,nRays,nBins,rScale,rStart)
				nElAngle - list of elevation angles
				nRays - list of number of rays per elevation
				nBins - list of number of bins per elevation
				rScale - list of range scales per elevation
				rStart - list of range starts per elevation
	"""
	nElAngle = []
	nRays = []
	nBins = []
	rScale = []
	rStart = []
	
	for i in range(nElevations):
		 where = f["dataset%d/where" % (i + 1)]
		 nRays.append(where.attrs["nrays"])
		 nBins.append(where.attrs["nbins"])
		 rScale.append(where.attrs["rscale"])
		 rStart.append(where.attrs["rstart"])
		 nElAngle.append(where.attrs["elangle"])
		 
	return(nElAngle,nRays,nBins,rScale,rStart)
	
def get_site_latlon(f):
	"""
	Gets the lon and lat coordinates of the radar site from the HDF5 file
		input: 	f - file object of h5py library;
		output: (lon, lat)
	"""
	where = f["where"]
	return(where.attrs["lon"],where.attrs["lat"])
	
def get_elevations_number(f):
	"""
	Gets the number of elevation the file given by the filename input variable
		input: 	filename -  string with the filename input variable;
		output: integer the number of elevations in the file
	"""
	DatasetIndex = 0
	NumExpr = re.compile('dataset\d+')
	for fname, fvalue in f.iteritems():
		IsDataset = NumExpr.match(fname)
		if IsDataset:
			DatasetIndex += 1	
	
	return(DatasetIndex)
	
def get_timestamp(FileName):
	"""get the timstamp from the name of the file
	   returns the date as a string and as a datetime object
	"""
	FileDate = re.search(r'\d{14}', FileName) 
	timestamp = dt.datetime.strptime(FileDate.group(),"%Y%m%d%H%M%S")
	return (FileDate.group(0),timestamp)
	
#-----------------------------------------------------------------------------------------
# Extra functions for on_off_shore module:

def inpolygon(polygon, xp, yp):
    """
    Gets the points xp and yp, that are inside of the polygon
    """
    return np.array([Point(x, y).intersects(polygon) for x, y in zip(xp, yp)],dtype=np.bool)
                    

#-----------------------------------------------------------------------------------------
# Plotting functions for vis_refl_vol2bird module

def plot_ppi(data, (radlon, radlat),
             site=(0, 0), proj=None, elev=0.,
             ax=None,basemap=None,vmin=None,vmax=None,cmap=None, norm=None ):
    """Plots a Plan Position Indicator (PPI). (modified from wradlib)

    .. versionchanged:: 0.6.0
       using osr objects instead of PROJ.4 strings as parameter

    The implementation of this plot routine is in cartesian axes and does all
    coordinate transforms beforehand. This allows zooming into the data as well
    as making it easier to plot additional data (like gauge locations) without
    having to convert them to the radar's polar coordinate system.


    Parameters
    ----------
    data : np.array
        The data to be plotted. It is assumed that the first dimension is over
        the azimuth angles, while the second dimension is over the range bins

    site : tuple
        Tuple of coordinates of the radar site.
        If `proj` is not used, this simply becomes the offset for the origin
        of the coordinate system.
        If `proj` is used, values must be given as (longitude, latitude)
        tuple of geographical coordinates.
    proj : osr spatial reference object
        GDAL OSR Spatial Reference Object describing projection
        If this parameter is not None, `site` must be set. Then the function
        will attempt to georeference the radar bins and display the PPI in the
        coordinate system defined by the projection string.
    elev : float or array of same shape as az
        Elevation angle of the scan or individual azimuths.
        May improve georeferencing coordinates for larger elevation angles.
    ax : matplotlib Axes object
        If given, the PPI will be plotted into this axes object. If None a
        new axes object will be created
    basemap : a basemap object for the reprojectio of the radar data

    Returns
    -------
    ax : matplotlib Axes object
        The axes object into which the PPI was plotted
    pm : matplotlib QuadMesh object
        The result of the pcolormesh operation. Necessary, if you want to
        add a colorbar to the plot.

    """
   
    (bmradlon, bmradlat) = basemap(radlon, radlat)

    # get the current axes.
    # this creates one, if there is none
    if ax is None:
        ax = pl.gca()

    # plot the colormesh
    pm = basemap.pcolormesh(bmradlon, bmradlat, data,vmin=vmin,vmax=vmax,cmap=cmap, norm=norm)
    

    return ax, pm

def make_map(projection='aeqd', resolution='i'):
	"""
	Plot the map for the Jabekke radar domain.
	"""
	fig, ax = plt.subplots()
	ax.set_axis_bgcolor('red')
	m = Basemap(projection='aeqd',lon_0=3.0642,lat_0=51.1917,width=212132.,height=212132.,resolution=None, ax=ax)
	m.shadedrelief()
	parallels = np.arange(50.5, 52.8, 1.)
	meridians = np.arange(0., 8., 2.)
	paral = m.drawparallels(parallels, labels=[1, 0, 0, 0],fontsize=10, color='black', ax=ax)
	for pi in paral.items():
            if len(pi[1][1])==1: pi[1][1][0].set_color('black')
	merid = m.drawmeridians(meridians, labels=[0, 0, 0, 1],fontsize=10, color='black', ax=ax)
	for mi in merid.items():
            
            if len(mi[1][1])==1: mi[1][1][0].set_color('black')
	return fig, ax, m

def drawShapesFromFile(filename,facecolor,edgecolor,m,ax):
	"""
	Draw the shapefiles.
	"""
    m.readshapefile(filename, 'temp', drawbounds = False)
    patches = []
    for info, shape in zip(m.temp_info, m.temp): patches.append(Polygon(np.array(shape), True) )
    ax.add_collection(PatchCollection(patches, facecolor=facecolor, edgecolor=edgecolor, linewidths=1,color='grey'))

def plot_refl(data,(radlon, radlat), coordinates, elv,ylabel = 'reflectivity (dBZ)', maxrange=35000,fdate='07-06-2016 15:20',ax=None, basemap=None, range50km=None, shore='offshore', fulldata=True):
	"""
	Plot refl for ground station - sattelite link approximation
	"""
	cMap = colors.ListedColormap(['lightgrey','cyan','deepskyblue','dodgerblue','blue','greenyellow','chartreuse','limegreen','green','darkgreen','yellow','orange','darkorange','red','darkred','magenta'])
	if (ylabel == 'reflectivity (dBZ)'):
		clevs = range(-16,32,3)
		norm = colors.BoundaryNorm(clevs, cMap.N)
		data[data==95.5] = -32.
		data[data==0.0] = -32.
		data[data>32.] = 32. 
		data = np.ma.masked_array(data, data == np.nan)
		data = np.ma.masked_array(data, data < -16.)
		vmin=-16
		vmax=32
	else:
		clevs = range(0,34,2)
		norm = colors.BoundaryNorm(clevs, cMap.N)
		data[data > 34.] = 34. 
		data = np.ma.masked_array(data, data == np.nan)
		data = np.ma.masked_array(data, data < 0.)
		vmin=0
		vmax=34
	
	ax, cf = plot_ppi(data, (radlon, radlat), vmin=vmin, vmax=vmax, cmap=cMap, norm=norm, basemap=basemap,ax=ax)
	if(fulldata):ofshroeText = ""
	else:
            if(shore is "offshore"):ofshroeText = ", off-shore part"
            else:ofshroeText = ", on-shore part"
	if(range50km):
            titleobj = plt.title("Jabbeke radar %.2f%s el. (50 km range%s)\n %s " %(elv,unichr(176),ofshroeText,fdate), color='black', size=12)
        else:
            titleobj = plt.title("Jabbeke radar %.2f%s el., %s %s " %(elv,unichr(176),ofshroeText,fdate), color='black',size=12)

	heatmap = plt.pcolor(data, cmap=cMap, norm=norm )
	cbar = basemap.colorbar(heatmap,"right", size="5%", pad='2%')
	cbar.ax.get_yaxis().set_ticks([])
	titleFontName = titleobj.get_fontname()
	for j, lab in enumerate(['-16','-13', '-10', '-7','-4','-1','2','5', '8', '11', '14','17','20','23','26','29']):
            cbar.ax.text(1.6, (1. * j + .5) / 16.0, lab, ha='center', va='center', color='black', family=titleFontName, size=10)
	cbar.ax.get_yaxis().labelpad = 30
	cbar.ax.set_ylabel(ylabel, rotation=90,color='black',size=10)
	(radarlon,radarlat) = basemap(coordinates[0],coordinates[1])
	basemap.plot(radarlon,radarlat,'wo',5)
	#plot bird radar ster with coordinates 51 31'58.28"N 2 57'17.57"E
	(brlon, brlat) = basemap(2.9548806,51.5328556)
	basemap.plot(brlon, brlat,'k*',3)
	plot_range_rings(([5000,35000]),(['5 km', '35 km']),(radarlon,radarlat), basemap=basemap, ax=ax, col="white")
	return(plt)


def plot_range_rings(range_rings, labels, (radlon, radlat), basemap=None,ax=None, col='k', ls='-', lw=2):
	"""
	Plot a series of range rings.
	Parameters
	----------
	range_rings : list
		List of locations in km to draw range rings.
	ax : Axis
		Axis to plot on.  None will use the current axis.
	col : str or value
		Color to use for range rings.
	ls : str
		Linestyle to use for range rings.
	"""
	for range_ring_location_km,label in zip(range_rings,labels):
		plot_range_ring(range_ring_location_km, label,(radlon, radlat),basemap=basemap, ax=ax, col=col, ls=ls, lw=lw)


def plot_range_ring(range_ring_location_km, label, (radlon, radlat), npts=100,basemap=None, ax=None, col='black', ls='-', lw=1):
        """
        Plot a single range ring.
        Parameters
        ----------
        range_ring_location_km : float
            Location of range ring in km.
        npts: int
            Number of points in the ring, higher for better resolution.
        ax : Axis
            Axis to plot on.  None will use the current axis.
        col : str or value
            Color to use for range rings.
        ls : str
            Linestyle to use for range rings.
        """
        theta = np.linspace(0, 2 * np.pi, npts)
        r = np.ones([npts], dtype=np.float32) * range_ring_location_km
        (xx,yy) = (radlon, radlat)
        x = xx + r * np.sin(theta)
        y = yy + r * np.cos(theta)
        basemap.plot(x, y, c=col, ls=ls,linewidth=1.6)
