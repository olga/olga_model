import datetime
import numpy as np

# Read comma-separated (longname, shortname, lon, lat) file
class readLocations:    
    def __init__(self, txtfile):
        filein = np.genfromtxt(txtfile, dtype='str', delimiter=',', skip_header=0)
        # Read columns
        self.longName  = filein[:,0]
        self.shortName = filein[:,1]  
        self.lon       = np.array(filein[:,2], dtype=np.float32)  
        self.lat       = np.array(filein[:,3], dtype=np.float32)  
        self.type      = np.array(filein[:,4], dtype=np.int)  
        # Strip spaces
        self.longName  = [s.strip(' ') for s in self.longName]
        self.shortName = [s.strip(' ') for s in self.shortName]

# Main settings for OLGA
class olgaSettings:
    def __init__(self):
        # ----------------------------------
        # Local file system settings.
        # Full path to directory of this script. Append with '/' !!
        
        if(True): # Mint ---------------------
            self.olgaRoot     = '/home/bart/WRFnl/olga_model/' # Mint 
            self.domainRoot   = '/home/bart/WRFnl/olga_model/OLGA_NL/' # Mint 

            # The following directories are by default defined relative to the current directory. 
            # However, also absolute paths, at other disks/partitions/etc. are possible to store
            # e.g. the large input and output files somwhere else.
            self.wpsRoot      = self.domainRoot + 'WPS/'        # Path to root of WPS run directory
            self.wrfRoot      = self.domainRoot + 'WRF/'        # Path to root of WRF run directory
            self.olgaLogs     = self.domainRoot + 'logs/'       # Location to save logs
            self.figRoot      = '/home/scratch1/WRFnl/outputOLGA/' # Path to save OLGA output
            self.wrfDataRoot  = '/home/scratch1/WRFnl/outputWRF/'  # Path to store the WRF output
            self.gfsDataRoot  = '/home/scratch1/WRFnl/inputGFS/'   # Path to store the GFS data

        if(False): # Thunder -----------------
            self.olgaRoot     = '/scratch/mpi/mpiaes/m300241/WRFnl/olga_model/' # Thunder
            self.domainRoot   = '/scratch/mpi/mpiaes/m300241/WRFnl/olga_model/OLGA_NL/' # Thunder

            # The following directories are by default defined relative to the current directory. 
            # However, also absolute paths, at other disks/partitions/etc. are possible to store
            # e.g. the large input and output files somwhere else.
            self.wpsRoot      = self.domainRoot + 'WPS/'        # Path to root of WPS run directory
            self.wrfRoot      = self.domainRoot + 'WRF/'        # Path to root of WRF run directory
            self.olgaLogs     = self.domainRoot + 'logs/'       # Location to save logs
            self.figRoot      = self.domainRoot + 'outputOLGA/' # Path to save OLGA output
            self.wrfDataRoot  = self.domainRoot + 'outputWRF/'  # Path to store the WRF output
            self.gfsDataRoot  = self.olgaRoot   + 'inputGFS/'   # Path to store the GFS data

        # ----------------------------------
        # Read ASCII file with locations, here used for everything
        loc1 = readLocations(self.domainRoot + 'locations.txt')
        #loc2 = readLocations(self.domainRoot + 'airports.txt')

        # ----------------------------------
        # Computational settings. 
        self.mpiTasks    = 4 # Number of MPI tasks
        self.ompThreads  = 1 # Number of OpenMP threads

        # ----------------------------------
        # Number of domains. This can be less than the size of the arrays below
        # in which case only the first ndom are used (e.g. for quick testing of outer domain)
        self.ndom         = 1  

        # ----------------------------------
        # Time settings
        self.ttotal       = 24 # Total time to simulate [h]
        self.tslice       = 24 # Split 'ttotal' in 'tslice' chunks [h]
        self.dt_output    = ([60,60]) # 'history_interval' from namelist, per domain, in minutes

        # ----------------------------------
        # Manually specify times (UTC) over which to make PFD's and time series. Only if both times are within
        # one 'tslice', maps are made. BvS: add better description :)
        # Same 'tanalysis' is used for all domains!
        self.tanalysis    = ([2,21]) 

        # ----------------------------------
        # Main map settings per domain
        self.maps         = ([True,False]) # Make maps or not

        # By setting map_lat, map_lon, map_width and map_height to -1, OLGA automatically
        # tries to determine the best settings. Useful for setting up domains and first tests
        self.map_lat      = ([51.2,-1]) # Central latitude of map [deg]
        self.map_lon      = ([7.9,-1]) # Central longitude of map [deg]
        self.map_width    = ([690000,-1]) # Domain plot width [m]
        self.map_height   = ([680000,-1]) # Domain plot height [m]
        self.map_res      = (['i','i']) # Details of map (c=crude, l=low, i=interm, h=high)
        self.drawRivers   = ([True,False]) # Draw rivers
        self.drawCities   = ([True,False]) # Draw cities
        self.cityLoc      = ([loc1, loc1]) # Location of city/airport/.. markers on maps
        self.map_desc     = (['6x6km','2x2km'])

        # Figure settings. For maps, the height is determined from the map aspect ratio
        self.fig_width_px  = 650  # width of figure [px]
        self.fig_dpi       = 100  # resolution in [px/in]
        self.map_bottom_px = 20   # bottom margin [px]
        self.map_top_px    = 25   # top margin in [px]
        self.map_left_px   = 12   # left margin in [px]
        self.map_right_px  = 80   # right margin in [px]

        # ----------------------------------
        # Plot variables maps
        vars1 = (['pfd','wstar','zidry','cudepth','convection','clouds','rain','wind10m','wind1000m'])         
        vars2 = (['pfd','wstar','zidry','convection','clouds','cudepth'])         
        #vars1 = (['pfd'])         
        #vars2 = (['rr'])         
        self.map_vars     = ([vars1,vars2]) # variables to plot per domain

        # ----------------------------------
        # Settings soundings (Detailed settings are in src/skewtlogp.py)
        self.sounding     = ([True,False]) # Make soundings or not
        self.soundLoc     = ([loc1, loc1]) # Location of soundings

        # Settings time series (meteograms)
        self.meteogr      = ([True,False]) # Make meteograms or not
        self.meteogrLoc   = ([loc1, loc1]) # Location of meteograms

        # -----------------------------------------------------
        # Don't change below, some checking of input
        if(self.ttotal%self.tslice != 0):
            sys.exit('ttotal should be integerer multiple of tslice')

    # Use time structs to do calculations on time  
    def set_time(self,islice):
        self.islice      = islice
        self.abs_start   = datetime.datetime.strptime('%i %i %i %i'%(self.day,self.month,self.year,self.tstart),"%d %m %Y %H")
        self.t0          = self.tstart + (self.islice+0) * self.tslice
        self.t1          = self.tstart + (self.islice+1) * self.tslice
        self.startstruct = self.abs_start + datetime.timedelta(hours=self.t0)
        self.endstruct   = self.abs_start + datetime.timedelta(hours=self.t1)

