#!/opt/baltrad/third_party/bin/python
import sys, os, glob, time, re, traceback
import odim_source
import _raveio
import _polarvolume
import _pyvol2bird
import h5py

VOL2BIRD_OUTPUTDIR="/tmp"

SOURCES=  ["bejab", "bewid", "bezav", 
    "dkbor", "dkrom", "dksin", "dkste", "dkvir", 
    "eehar", "eesur", 
    "hrbil", "hrosi",
    "fianj", "fiika", "fikes", "fikor", "fikuo",
    "filuo", "fiuta", "fivan", "fivim", 
    "frabb", "frale", "frbla", "frbol", "frbor", "frbou",
    "frcae", "frche", "frcol", "frgre", "frlep",
    "frmcl", "frmom", "frmtc", "frnan", "frnim",
    "fropo", "frpla", "frtou", "frtra", "frtre",
    "nldbl", "nldhl", 
    "searl", "sease", "sehud", "sekir", "sekkr", "selek",
    "selul", "seosu", "seovi", "sevar", "sevil", "seang",
    "silis", "sipas",
    "czska", "czbrd",
    "eslid", "esmad", "esval", "esmur", "esbar", "escor",
    "eszar", "essan", "essse", "esbad", "esmal", "essev",
    "esalm", "espma", "eslpa",
    "deemd", "deham", "deboo", "deros", "dehan", "dehnr", 
    "debln", "depro", "deess", "defld", "deumd", "deoft",
    "deneu", "dedrs", "denhb", "detur", "deeis", "deflg",
    "demuc", "desna", "demem",
    "plleg", "plram", "plpas", "plrze", "plpoz", "plswi",
    "plgda", "plbrz",
    "nober", "noosl", "nohgb", "norsa", "norst", "nobml",
    "noand", "nohas", "nosta", "nohur", 
    "ukcle", "ukham", "ukche", "ukcas", "ukpre", "uking",
    "ukcyg", "ukdud", "uklew", "ukcob", "ukhhd", "ukmun",
    "ukthu", "ukdea", "ukhmy", 
    "chalb", "chdol", "chlem",
    "skjav", "skkoj"]

def add_daylight(t):
    l = list(t)
    l[8] = -1
    return tuple(l)

# returns the archive time as a time tuple which can be fed to GetTimeTuple.
def GetTime():
  t = time.gmtime(time.time() - 15*60)
  return add_daylight(t)

def print_log(msg):
  tt = time.gmtime(time.time())
  tt = add_daylight(tt)
  logtime = time.strftime("%Y-%m-%d %H:%M", tt)
  print "%s: %s"%(logtime, msg)

def getNod(src):
    return odim_source.NODfromSource(src)

def get_filename(src):
  dd = src.date
  tt = src.time

  year = dd[:4]
  month = dd[4:6]
  day = dd[6:8]
  hour = tt[:2]
  minute = tt[2:4]

  return "%s_vp_%s%s%s%s%s00.h5"%(getNod(src), year, month, day, hour, minute)

def write_file_to_outdir(out_dir, outname, dst):
  fname = "%s/%s"%(out_dir, outname)
  #print fname
  rio = _raveio.new()
  rio.object=dst
  rio.save(fname)
  rio.close()
  print_log("Wrote %s"%fname)

def process_file(out_dir, in_name):
  print "in process_file"
  print in_name
  print h5py.File(in_name, 'r')

  vol = _raveio.open(in_name).object
  if _polarvolume.isPolarVolume(vol):
    if getNod(vol) in SOURCES:
      fname = get_filename(vol)     
      v2b = _pyvol2bird.new(vol)
      # calculate a vertical profile of birds
      vpr = v2b.vol2bird(vol)
      write_file_to_outdir(out_dir, fname, vpr)
  else:
    print_log("Not a polar volume. Ignoring")

def GetODCPath(tt):
  import odc_filesys
  return odc_filesys.GetPathFromTuple(tt)

def generate(out_dir, in_dir=None, strdate=None):
    tt = GetTime()

    start = time.time()

    # Which PVOLs do we already have before starting?
    path = in_dir
    if path == None:
      path = GetODCPath(tt)

    print_log("PATH=%s"%path)
    if(strdate is None):
		fstrs = sorted(glob.glob(os.path.join(path, '*vol*.h5')))
    else:
		fstrs = sorted(glob.glob(os.path.join(path, '%s*vol*.h5' %strdate)))
		
    for f in fstrs:
      FileDate = re.search(r'\d{12}', f)
      foutstr = sorted(glob.glob(os.path.join(out_dir, '*_vp_%s00.h5' %FileDate.group())))

      if (len(foutstr) == 0):
      	try:
      		print "process_file(%s, %s)" %(out_dir, f)
      		process_file(out_dir, f)
      	except Exception, e:
      		traceback.print_exc()

if __name__ == "__main__":
  from optparse import OptionParser
  usage = "usage: vol2bird_generate -o <output dir>"
  usage += ""
  parser = OptionParser(usage=usage)
  parser.add_option("-o", "--outdir", dest="out_dir",
                    default=VOL2BIRD_OUTPUTDIR,
                    help="Name of output directory to which to write.")
  parser.add_option("-i", "--indir", dest="in_dir",
                    default=None,
                    help="Name of input directory to scan for volumes.")
  parser.add_option("-t", "--datetime", dest="strdate",
                    default=None,
                    help="Date at the beginning of the name of input volumes.")
  (options, args) = parser.parse_args()
  generate(options.out_dir, options.in_dir,strdate=options.strdate)
