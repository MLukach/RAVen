#! /usr/bin/env python
"""
 mtr_vbird.py
 Module generates MTR's for a given timestamp, based on the vertical profiles that are found in the input directory and containe the timestamp in the name of the files. 
 requiers wradlib.env to be sourced on gloin
 input: -h for help
		-o output directory for the MTR
		-i input directory with vertical bird prophile files
		-t timestamp for the generation of the MTR's 
		
 call python mrt_vbird.py -i ./h5.merged/offshore/ -o ./MTR/ -t YYYYmm(DD)

 ___________________________________
 
 Maryna Lukach (RMI), 2017 
 
"""
import sys,tarfile
import os
import glob
import datetime as dt
from optparse import OptionParser
import time
import pytz
import numpy as np
import re
import h5py
import pandas as pd
import math


def read_dataset(h5File, GroupName = '/dataset1/data1',undetval=0):
	"""
	Read dataset DatasetName from DfltInputPath file and return it as np.array 
	"""
	Dataset = h5File[GroupName + "/data"]
	Data = np.zeros(Dataset.shape,dtype = Dataset.dtype)
	Data = Dataset
	WhatGroup = GroupName + "/what"
	Data = decode_data(Data,h5File[WhatGroup].attrs["gain"],h5File[WhatGroup ].attrs["offset"],nodata=h5File[WhatGroup].attrs["nodata"], undetect=h5File[WhatGroup].attrs["undetect"],undetval=undetval)
         
	return(Data)

def theta_from_vectors(u,v):
	"""
	Calculate theta from u and v vector components
	"""
	Theta = np.arctan2(u, v) * 180 / np.pi
	if (np.isscalar(Theta)):
		if(Theta<0): Theta = 360 + Theta
	else:
		if(len(Theta < 0)):
			Theta[Theta<0] = 360 + Theta[Theta<0]
	return(Theta)

def decode_data(data,gain,offset,nodata=None,noval=np.nan,undetect=None,undetval=np.nan):
        """Decoded data into real value using a gain offset and nodata.
           By default NODATA is set to NAN values for convenience.        
        """
        data = offset + data * gain        
        if nodata is not None:
			data[data == (nodata + offset * gain)] = noval
        if undetect is not None:
			data[data == (undetect + offset * gain)] = undetval
		
        return(data)


def check_path(DirName, Create = False):
	"""
	Check if the given path exists and create it if required.
	"""
	if (os.path.exists(DirName) == False):
		if (Create) : os.system('mkdir -p %s' %(DirName))
		else:
			print('Directory %s does not exist and it is not allowed to creat it!' % DirName) 
			DirName = './'
	return DirName

def append_means_list(listOfMeans,a,strName):
	"""
	Add another mean value to the least of means
	"""
	nonNans = np.count_nonzero(~np.isnan(np.asarray(a[strName]))) 
	if(nonNans>0):listOfMeans.append(np.nansum(np.asarray(a[strName]))*1./nonNans)
	else:listOfMeans.append(np.nan)
	return listOfMeans

def nan_mean(array, axis=None):
	"""
	Calculate the mean of not NaN values
	"""
	marray = np.ma.masked_array(array,np.isnan(array))
	array = np.mean(marray,axis=axis)
	
	return array
	
def mean_theta(listToAppend,U,V,nonNans):
	"""
	Get mean theta
	"""
	if(nonNans > 0):
		MeanUmean = nan_mean(U)
		MeanVmean = nan_mean(V)
		MeanTheta = theta_from_vectors(MeanUmean,MeanVmean)
		listToAppend.append(MeanTheta)
	else:
		listToAppend.append(np.nan)
	return 	listToAppend			

def get_domain(word,directory):
	"""
	Find given word in the directory name
	"""
	Word = re.search(r'\b%s\b' %word,directory)
	if(Word):Domain = Word.group(0)
	else:Domain = None
	return Domain

def get_timestamp(FileName):
	"""get the timstamp from the name of the file
	   returns the date as a string and as a datetime object
	"""
	FileDate = re.search(r'\d{14}', FileName) 
	timestamp = dt.datetime.strptime(FileDate.group(),"%Y%m%d%H%M%S")
	return (FileDate,timestamp)

def get_nexttimestamp(VBDirFiles,FileIndex,days,FileDate,timestamp):
	
	# get next time stamp from the name of the next file in the ordered list of files
	if (FileIndex<len(VBDirFiles)-1):
		nextFileDate = re.search(r'\d{14}', VBDirFiles[FileIndex+1])
		nexTimestamp = dt.datetime.strptime(nextFileDate.group(),"%Y%m%d%H%M%S")
	elif (FileIndex==(len(VBDirFiles)-1)):
		nexTimestamp = dt.datetime.strptime(FileDate.group()[0:10] + "0000","%Y%m%d%H%M%S")
	# add one day if next file belongs to the next day
	if (timestamp.day != nexTimestamp.day):
		days += 1
	return (days, nexTimestamp)

def get_domain_from_directory(VBirdDirectory):
	Domain = get_domain("full_vp",VBirdDirectory)
	WindDataDirectory = VBirdDirectory
	fullDomain = True
	if(Domain is None):
		Domain = get_domain("onshore_vp",VBirdDirectory)
		if(Domain is None): 
			Domain = get_domain("offshore_vp",VBirdDirectory)
			if (Domain is None): Domain = ""
			else: 
				WindDataDirectory = VBirdDirectory.replace("offshore_vp","full_vp")
				fullDomain = False
		else: 
			WindDataDirectory = VBirdDirectory.replace("onshore_vp","full_vp")
			fullDomain = False
	return (Domain, fullDomain, WindDataDirectory)
					

def main():
	
	usage = """%prog -i path/to/input/directory/v -o path/to/output/directtory/Newfilename -t YYYYmm(DD)
	"""
    # determine basename of the application
	appbasename = os.path.basename(os.path.splitext(sys.argv[0])[0])
    
    # determine the basedir  of the application )
	basedir = os.path.dirname(os.path.abspath(sys.argv[0])) + os.sep
	
     # parsing options and arguments
	parser = OptionParser(usage=usage)
	parser.add_option("-d", "--debug", action="store_true", dest="debug", default=False, help="output debug messages during program execution")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="print messages (debug,info,warning,...) about the program execution to the console (stderr)")
	parser.add_option("-i", "--input directory with vertical bird prophile files", dest = "VBirdDirectory", default = "./tests/output/50kmvp.201609.single/", help = "path to the vbird files with the temperature profiles") 
	parser.add_option("-o", "--output directory for the MTR", dest = "MTRDirectory", default = "./tests/output/50kmvp.201609.single/mtr/", help = "path to the directory for the MTR files")
	parser.add_option("-t", "--date for MTR calculation", dest = "MTRDate", default = "20160915", help = "date for the selection of the MTR files")
	
	(options, files) = parser.parse_args()
	# get parameters values
	verbose = options.verbose
	debug = options.debug
	VBirdDirectory = options.VBirdDirectory
	MTRDirectory = options.MTRDirectory
	MTRDate = options.MTRDate
	
    # check if MTR's directory exists and create it if not
	check_path(MTRDirectory, Create = True)
	
	# get sorted list of files
	VBDirFiles = sorted(glob.glob(os.path.join(VBirdDirectory, 'bejab_vp_%s*.h5' %MTRDate)))
	
	utc = pytz.timezone("UTC")
	
	# predefine arrays (20 altitude levels between 200m and 4000m a.m.s.l.)
	SpecMTR = np.zeros([len(VBDirFiles), 20,1],dtype = float)
	SpecMTR_real = np.zeros([len(VBDirFiles), 20,1],dtype = float)
	NPoints = np.zeros([len(VBDirFiles), 20,1],dtype = float)*np.nan
	SpeedData = np.zeros([len(VBDirFiles), 20,1],dtype = float)*np.nan
	Correction = np.zeros([20,1],dtype = float)
	OpositeCorrection = np.zeros([20,1],dtype = float)
	
	# predefine arrays for MTR values per altitude interval between 200 and 800; between 800 and 1400; on an interval between 200 and 1400
	MTR200_800m = np.zeros([len(VBDirFiles), 1],dtype = float)
	MTR800_1400m = np.zeros([len(VBDirFiles), 1],dtype = float)
	MTR200_1400m = np.zeros([len(VBDirFiles), 1],dtype = float)
	
	# predefine arrays for speed values per altitude interval between 200 and 800; between 800 and 1400; on an interval between 200 and 1400
	speed200_800m = np.zeros([len(VBDirFiles), 1],dtype = float)
	speed800_1400m = np.zeros([len(VBDirFiles), 1],dtype = float)
	speed200_1400m = np.zeros([len(VBDirFiles), 1],dtype = float)
	
	FlyDirectionList = []
	
	# 12 timestamps per hour and 20 altitude levels between 200m and 4000m a.m.s.l.
	FlyVectorU = np.zeros([12, 20,1],dtype = float)
	FlyVectorV = np.zeros([12, 20,1],dtype = float)
	
	Theta = np.empty([20,1])*np.nan
	
	AllLevelsTheta = []
	
	DatesList = []
	DayList = []
	FlyDirectionList200_800m = []
	FlyDirectionList800_1400m = []
	FlyDirectionList200_1400m = []

	FlyDirectionProfileList = []
	FlyDateList = []
	ts = 0
	days = 0
	total_count_of_deleted = 0
	total_count_of_all = 0
	
	(Domain, fullDomain, WindDataDirectory) = get_domain_from_directory(VBirdDirectory)
	arc = "%s/%s_%s.tar.gz" %(WindDataDirectory,os.path.basename(os.path.normpath(WindDataDirectory)),MTRDate)
	wind_fiels_extracted = False
	if(os.path.isfile(arc)): 
		tar = tarfile.open(arc)
		tarlist = tar.getnames()
		if (not os.path.isfile("%s/%s" %(WindDataDirectory,tarlist[0]))):
			tar.extractall(path=WindDataDirectory)
			wind_fiels_extracted = True
		tar.close()
	
	# for each file in the sorted list of selected files
	if (len(VBDirFiles) != 0):
		
		for FileIndex in range(0,len(VBDirFiles)):
			FileName = VBDirFiles[FileIndex]
			# get the timstamp from the name of the file for this ans next file (if this one is not the last one)
			(FileDate,timestamp) = get_timestamp(FileName)
			(days, nexTimestamp) = get_nexttimestamp(VBDirFiles,FileIndex,days,FileDate,timestamp)
			
			# set the domain
			(Domain, fullDomain, WindDataDirectory) = get_domain_from_directory(VBirdDirectory)
					
			DatesList.append(timestamp.strftime("%Y%m%d%H%M%S"))
			DayList.append(timestamp.strftime("%Y%m%d%H0000"))
			
			# open the file
			f = h5py.File(FileName, 'r')
			# read Dencity as eta (cm^2/km^3)/11 (birds/cm^2) and NPoints outputs of the vol2bird algorithm
			DencityDataset = np.divide(read_dataset(f, GroupName = '/dataset1/data7'),11.0)
			
			
			ThisNPoints = read_dataset(f, GroupName = '/dataset1/data6')
			#print "ThisNPoints"
			#print ThisNPoints
			# for the on- and off- shore parts u and v vectors of full domain are used (read them from the correspondent file of the full domain)
			#print WindDataDirectory+"/bejab_vp_"+FileDate.group()+".h5"
			if(fullDomain is False):
				try:
					g = h5py.File(WindDataDirectory+"bejab_vp_"+FileDate.group()+".h5",'r')
					# in full file the dencity dataset is opend ('/dataset1/data11')
					FullDencity = read_dataset(g, GroupName = '/dataset1/data11')
					FullNPoints = read_dataset(g, GroupName = '/dataset1/data6')
				except Exception as error:
					print ('Oops, something went wrong (%s) opening of file %s, terminating ...' % (str(error), WindDataDirectory+"bejab_vp_"+FileDate.group()+".h5"))
				
				# Create right name of the opositefile (on for offshore and vise versas)
				if(FileName.find("off")>-1):oposite = FileName.replace("off","on")
				elif(FileName.find("on")>-1):oposite = FileName.replace("on","off")
				
				h = h5py.File(oposite, 'r')
				OpositeNPoints = read_dataset(h, GroupName = '/dataset1/data6')
				# For oposite file in place of dencity we also use eta
				OpositeDencity = np.divide(read_dataset(h, GroupName = '/dataset1/data7'),11.0)
				
				OtherNPoints = OpositeNPoints + ThisNPoints
				
				
			else:
				g = f
				OtherNPoints = ThisNPoints
				
				FullDencity = read_dataset(g, GroupName = '/dataset1/data11')
			
				FullNPoints = read_dataset(g, GroupName = '/dataset1/data6')
				
			
			# Clean up the data according to the FullDencity values
			# Set Dencity to 0 if FullDencity is 0 and the Dencity is nonzero and not a NaN
			DencityDataset[(FullDencity == 0.) & (np.isnan(DencityDataset) == False) & (DencityDataset > 0.)] = 0.
			# Set Dencity to NaN if FullDencity is NaN
			DencityDataset[np.isnan(FullDencity) == True] = np.nan
			
			if(fullDomain is False):
				# Clean up the oposite data according to the FullDencity values
				# Set Dencity to 0 if FullDencity is 0 and the Dencity is nonzero and not a NaN
				OpositeDencity[(FullDencity == 0.) & (np.isnan(OpositeDencity) == False) & (OpositeDencity > 0.)] = 0.
				# Set oposite Dencity to NaN if FullDencity is NaN
				OpositeDencity[np.isnan(FullDencity) == True] = np.nan
				
				# Correction based on the number of points (to have the sum of points of partial be not more than 3% away from the full_number of pointss)
				DencityDataset[np.divide(np.abs(FullNPoints - OtherNPoints),FullNPoints, where=np.logical_and(np.isfinite(FullNPoints),FullNPoints>0.))>0.03] = np.nan
				
			# read speed data
			Speed = read_dataset(g, GroupName = '/dataset1/data1')

			# read vector u
			FlyVectorU[int(timestamp.minute/5),:,:] = read_dataset(g, GroupName = '/dataset1/data12', undetval=np.nan)
			# read vector v
			FlyVectorV[int(timestamp.minute/5),:,:] = read_dataset(g, GroupName = '/dataset1/data13', undetval=np.nan) 
			
			if(fullDomain is False):
			
				Correction[FullNPoints>0] = np.divide(ThisNPoints[FullNPoints>0],FullNPoints[FullNPoints>0])
				OpositeCorrection[FullNPoints>0] = np.divide(OpositeNPoints[FullNPoints>0],FullNPoints[FullNPoints>0])
			
			# calculate SpecMTR fpr each altitude level
			SpecMTR[FileIndex,:] = np.divide(np.multiply(np.multiply(DencityDataset,Speed),3.6),5.)
			
			SpecMTR_real[FileIndex,:] = np.divide(np.multiply(np.multiply(FullDencity,Speed),3.6),5.)#FullDencity * Speed * 3.6/5
			if(fullDomain is False):
				OpositeSpecMTR = np.divide(np.multiply(np.multiply(OpositeDencity,Speed),3.6),5)#OpositeDencity * Speed * 3.6/5	
				OnOffMTR = np.multiply(OpositeCorrection,OpositeSpecMTR) + np.multiply(Correction,SpecMTR[FileIndex,:])
				SpecMTR[FileIndex,np.logical_and(np.isnan(OnOffMTR),
						np.divide(np.abs(SpecMTR_real[FileIndex,:]-np.multiply(Correction,SpecMTR[FileIndex,:])),SpecMTR_real[FileIndex,:]))] = np.nan
					
				SpecMTR[FileIndex,(np.divide(np.abs(SpecMTR_real[FileIndex,:]-OnOffMTR),SpecMTR_real[FileIndex,:], where=np.logical_and(np.isfinite(SpecMTR_real[FileIndex,:]),SpecMTR_real[FileIndex,:]>0.)) > 0.03)] = np.nan
				
				total_count_of_deleted += len(SpecMTR[FileIndex,(np.divide(np.abs(SpecMTR_real[FileIndex,:]-OnOffMTR),SpecMTR_real[FileIndex,:], where=np.logical_and(np.isfinite(SpecMTR_real[FileIndex,:]),SpecMTR_real[FileIndex,:]>0.)) > 0.03)] )
			
			total_count_of_all += len(SpecMTR[FileIndex,:])
			
			# set speed data
			SpeedData[FileIndex,:,:] = Speed
			
			# write down the nuber of poins			
			NPoints[FileIndex,:] = ThisNPoints
			
			# Get MRT as a sum of specific values
			MTR200_800m[FileIndex] = np.nansum(SpecMTR[FileIndex,1:4])
			MTR800_1400m[FileIndex] = np.nansum(SpecMTR[FileIndex,4:7])
			MTR200_1400m[FileIndex] = np.nansum(SpecMTR[FileIndex,1:7])
			
			# Get meanSpead for lower bin [200,400), upper bin [800, 1400) and full interval [200, 1400)
			speed200_800m[FileIndex] = nan_mean(Speed[1:4])
			speed800_1400m[FileIndex] = nan_mean(Speed[4:7])
			speed200_1400m[FileIndex] = nan_mean(Speed[1:7])
				
			# fly direction calculation
			ts += 1 
			if((int(timestamp.hour * timestamp.day) != int(nexTimestamp.hour * nexTimestamp.day)) | (FileIndex == (len(VBDirFiles)-1))):
				
				MeanU = nan_mean(FlyVectorU,axis=0)
				MeanU = MeanU.filled(np.nan)
				
				MeanV = nan_mean(FlyVectorV,axis=0)
				MeanV = MeanV.filled(np.nan)
				
				relevIndex = np.logical_and((~np.isnan(MeanU)),(~np.isnan(MeanV)))
				Theta[relevIndex] = theta_from_vectors(MeanU[relevIndex],MeanV[relevIndex])
				
				# write down flight direction
				FlyDirectionList.extend(Theta.flatten()[:].tolist())
				
				ProfList = [timestamp.strftime("%Y%m%d%H0000")]
				ProfList.extend(Theta.flatten()[:].tolist())
				
				FlyDirectionProfileList.append(ProfList)
				
				nonNans200_800m = len(relevIndex[relevIndex[1:4,:] == True])
				nonNans800_1400m = len(relevIndex[relevIndex[4:7,:] == True])
				nonNans200_1400m = len(relevIndex[relevIndex[1:7,:] == True])
				
				mean_theta(FlyDirectionList200_800m,MeanU[1:4,:],MeanV[1:4,:],nonNans200_800m)
				mean_theta(FlyDirectionList800_1400m,MeanU[4:7,:],MeanV[4:7,:],nonNans800_1400m)
				mean_theta(FlyDirectionList200_1400m,MeanU[1:7,:],MeanV[1:7,:],nonNans200_1400m)	
				FlyDateList.append(timestamp.strftime("%Y%m%d%H0000"))
				
				ts = 0
				Theta = np.empty([20,1])
			
			f.close()
			if(fullDomain is False):
				g.close()
				h.close()
		
		# write down to dataframe calculated lists		
		df1 = pd.DataFrame({"A_Date" : DatesList, "B_MTR_0.2_0.8km" : MTR200_800m[:,0].tolist(), "C_MTR_0.8_1.4km" : MTR800_1400m[:,0].tolist(), "D_MTR_0.2_1.4km" : MTR200_1400m[:,0].tolist(), "E_speed_0.2_0.8km" : speed200_800m[:,0].tolist(), "F_speed_0.8_1.4km" : speed800_1400m[:,0].tolist(), "G_speed_0.2_1.4km" : speed200_1400m[:,0].tolist()})
		df1.to_csv(MTRDirectory + "bejab_vbird_mtr_" + MTRDate + '_%s.csv' %Domain, index=False)
		
               
		tmpdf = pd.DataFrame({"A_dayHour" : DayList, "B_MTR_0.2_0.8km" : MTR200_800m[:,0].tolist(), "C_MTR_0.8_1.4km" : MTR800_1400m[:,0].tolist(), "D_MTR_0.2_1.4km" : MTR200_1400m[:,0].tolist(), "E_speed_0.2_0.8km" : speed200_800m[:,0].tolist(), "F_speed_0.8_1.4km" : speed800_1400m[:,0].tolist(), "G_speed_0.2_1.4km" : speed200_1400m[:,0].tolist()})
		Theta = FlyDirectionList
		
		tmpdf2 = pd.DataFrame({"A_dayHour":DayList,"MTR_0.0km":SpecMTR[:,0,0].tolist(),"Npoints_0.0km":NPoints[:,0,0].tolist(),"Speed_0.0km":SpeedData[:,0,0].tolist(),"MTR_0.2km":SpecMTR[:,1,0].tolist(),"Npoints_0.2km":NPoints[:,1,0].tolist(),"Speed_0.2km":SpeedData[:,1,0].tolist(),"MTR_0.4km":SpecMTR[:,2,0].tolist(),"Npoints_0.4km":NPoints[:,2,0].tolist(),"Speed_0.4km":SpeedData[:,2,0].tolist(),"MTR_0.6km":SpecMTR[:,3,0].tolist(),"Npoints_0.6km" : NPoints[:,3,0].tolist(),"Speed_0.6km":SpeedData[:,3,0].tolist(),"MTR_0.8km":SpecMTR[:,4,0].tolist(),"Npoints_0.8km":NPoints[:,4,0].tolist(),"Speed_0.8km":SpeedData[:,4,0].tolist(),"MTR_1.0km":SpecMTR[:,5,0].tolist(),"Npoints_1.0km":NPoints[:,5,0].tolist(),"Speed_1.0km":SpeedData[:,5,0].tolist(),"MTR_1.2km":SpecMTR[:,6,0].tolist(),"Npoints_1.2km":NPoints[:,6,0].tolist(),"Speed_1.2km":SpeedData[:,6,0].tolist(),"MTR_1.4km":SpecMTR[:,7,0].tolist(),"Npoints_1.4km":NPoints[:,7,0].tolist(),"Speed_1.4km":SpeedData[:,7,0].tolist(),"MTR_1.6km":SpecMTR[:,8,0].tolist(),"Npoints_1.6km":NPoints[:,8,0].tolist(),"Speed_1.6km":SpeedData[:,8,0].tolist(),"MTR_1.8km":SpecMTR[:,9,0].tolist(),"Npoints_1.8km":NPoints[:,9,0].tolist(),"Speed_1.8km":SpeedData[:,9,0].tolist(),"MTR_2.0km":SpecMTR[:,10,0].tolist(),"Npoints_2.0km":NPoints[:,10,0].tolist(),"Speed_2.0km":SpeedData[:,10,0].tolist(),"MTR_2.2km":SpecMTR[:,11,0].tolist(),"Npoints_2.2km":NPoints[:,11,0].tolist(),"Speed_2.2km":SpeedData[:,11,0].tolist(),"MTR_2.4km":SpecMTR[:,12,0].tolist(),"Npoints_2.4km":NPoints[:,12,0].tolist(),"Speed_2.4km":SpeedData[:,12,0].tolist(),"MTR_2.6km":SpecMTR[:,13,0].tolist(),"Npoints_2.6km":NPoints[:,13,0].tolist(),"Speed_2.6km":SpeedData[:,13,0].tolist(),"MTR_2.8km" : SpecMTR[:,14,0].tolist(),"Npoints_2.8km":NPoints[:,14,0].tolist(),"Speed_2.8km":SpeedData[:,14,0].tolist(),"MTR_3.0km":SpecMTR[:,15,0].tolist(),"Npoints_3.0km":NPoints[:,15,0].tolist(),"Speed_3.0km":SpeedData[:,15,0].tolist(),"MTR_3.2km":SpecMTR[:,16,0].tolist(),"Npoints_3.2km":NPoints[:,16,0].tolist(),"Speed_3.2km":SpeedData[:,16,0].tolist(),"MTR_3.4km":SpecMTR[:,17,0].tolist(),"Npoints_3.4km":NPoints[:,17,0].tolist(),"Speed_3.4km":SpeedData[:,17,0].tolist(),"MTR_3.6km":SpecMTR[:,18,0].tolist(),"Npoints_3.6km":NPoints[:,18,0].tolist(),"Speed_3.6km":SpeedData[:,18,0].tolist(),"MTR_3.8km":SpecMTR[:,19,0].tolist(),"Npoints_3.8km":NPoints[:,19,0].tolist(),"Speed_3.8km":SpeedData[:,19,0].tolist()})
		MeanDatesList = tmpdf.A_dayHour.unique()
		
		HourlyMeanMTR200_800m = []
		HourlyMeanMTR800_1400m = []
		HourlyMeanMTR200_1400m = []
		
		HourlyMeanSpeed200_800m = []
		HourlyMeanSpeed800_1400m = []
		HourlyMeanSpeed200_1400m = []
		HourlyMeanSpeed = []
		
		HourlytFilesNumber = []
		HourlyMeanSpecMTR = []
		HourlyMeanNpoints = []
		HourlyMeanDirection = []
		HourlyLevel = []
		MeanSpecMTRDate = []
		
		# create lists of mean data per day
		for dayh in MeanDatesList:
				
			a = tmpdf.query("A_dayHour=='%s'" %dayh)
			append_means_list(HourlyMeanMTR200_800m,a,'B_MTR_0.2_0.8km')
			append_means_list(HourlyMeanSpeed200_800m,a,'E_speed_0.2_0.8km')
			append_means_list(HourlyMeanMTR800_1400m,a,'C_MTR_0.8_1.4km')
			append_means_list(HourlyMeanSpeed800_1400m,a,'F_speed_0.8_1.4km')
			append_means_list(HourlyMeanMTR200_1400m,a,'D_MTR_0.2_1.4km')
			append_means_list(HourlyMeanSpeed200_1400m,a,'G_speed_0.2_1.4km')
			
			b = tmpdf2.query("A_dayHour=='%s'" %dayh)
			
			for level in (0.2*np.arange(0,20)):

				append_means_list(HourlyMeanSpecMTR,b,'MTR_%.1fkm' %level)
				append_means_list(HourlyMeanNpoints,b,'Npoints_%.1fkm' %level)
				append_means_list(HourlyMeanSpeed,b,'Speed_%.1fkm' %level)
				HourlyLevel.append('%d' %(level*1000))
				MeanSpecMTRDate.append(dayh)
				
			HourlytFilesNumber.append(len(a.index))
		
		# write down to dataframes lists of means per day and save them to csv-files.	
		df2 = pd.DataFrame({"A_Date" : MeanDatesList, "B_meanMTR_0.2_0.8km" : HourlyMeanMTR200_800m, "C_meanMTR_0.8_1.4km" : HourlyMeanMTR800_1400m, "D_meanMTR_0.2_1.4km" : HourlyMeanMTR200_1400m, "E_NumberOfFiles" : HourlytFilesNumber, "F_meanDirection_0.2_0.8km" : FlyDirectionList200_800m, "G_meanDirection_0.8_1.4km" : FlyDirectionList800_1400m, "H_meanDirection_0.2_1.4km" : FlyDirectionList200_1400m, "I_HourlyMeanSpeed_0.2_0.8km" : HourlyMeanSpeed200_800m, "J_HourlyMeanSpeed_0.8_1.4km" : HourlyMeanSpeed800_1400m, "K_HourlyMeanSpeed_0.2_1.4km" : HourlyMeanSpeed200_1400m})
		df2.to_csv(MTRDirectory + "bejab_vbird_avgmtr_" + MTRDate + '_%s.csv' %Domain, index=False)
		df0 = pd.DataFrame({"A_Date": DatesList,"SpecMTR0":SpecMTR[:,0,0].tolist(),"SpecMTR1":SpecMTR[:,1,0].tolist(),"SpecMTR2":SpecMTR[:,2,0].tolist(),"SpecMTR3":SpecMTR[:,3,0].tolist(),"SpecMTR4":SpecMTR[:,4,0].tolist(),"Npoints0":NPoints[:,0,0].tolist(),"Npoints1":NPoints[:,1,0].tolist(),"Npoints2":NPoints[:,2,0].tolist(),"Npoints3":NPoints[:,3,0].tolist(),"Npoints4":NPoints[:,4,0].tolist()})
		df0.to_csv(MTRDirectory + "bejab_vbird_profmtr_" + MTRDate + '_%s.csv' %Domain, index=False)
		df3 = pd.DataFrame( FlyDirectionProfileList)
		df3.to_csv(MTRDirectory + "bejab_vbird_profdirect_" + MTRDate + '_%s.csv' %Domain, index=False)
		df4 = pd.DataFrame({"A_Date" :MeanSpecMTRDate,"B_MTR_level" : HourlyLevel, "C_MTR_value" : HourlyMeanSpecMTR, "D_Npoints":HourlyMeanNpoints, "E_NumberOfFiles":np.repeat(np.array(HourlytFilesNumber),20,axis=0).flatten()[:].tolist(),"F_meanDirection":FlyDirectionList, "G_HourlyMeanSpeed":HourlyMeanSpeed})
		df4.to_csv(MTRDirectory + "bejab_vbird_avgmtralllevels_" + MTRDate + '_%s.csv' %Domain, index=False)
		
#-----------------------------------------------------------------------------------------
if __name__ == "__main__":
        main()
