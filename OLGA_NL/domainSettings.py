import datetime

# Main settings for OLGA
class olgaSettings:
    def __init__(self):
        # ----------------------------------
        # Local file system settings.
        # Full path to directory of this script. Append with '/' !!
        #self.olgaRoot     = '/home/bart/WRFnl/olga_model/' # Mint 
        #self.domainRoot   = '/home/bart/WRFnl/olga_model/OLGA_NL/' # Mint 
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
        # Computational settings. 
        self.mpiTasks    = 4 # Number of MPI tasks
        self.ompThreads  = 1 # Number of OpenMP threads

        # ----------------------------------
        # Number of domains. This can be less than the size of the arrays below
        # in which case only the first ndom are used (e.g. for quick testing of outer domain)
        self.ndom         = 2  

        # ----------------------------------
        # Time settings
        self.ttotal       = 24 # Total time to simulate [h]
        self.tslice       = 24 # Split 'ttotal' in 'tslice' chunks [h]
        self.dt_output    = ([60,30]) # 'history_interval' from namelist, per domain, in minutes

        # ----------------------------------
        # Manually specify times (UTC) over which to make PFD's and time series. Only if both times are within
        # one 'tslice', maps are made. BvS: add better description :)
        # Same 'tanalysis' is used for all domains!
        self.tanalysis    = ([4,20])

        # ----------------------------------
        # Main map settings per domain
        self.maps         = ([True,True]) # Make soundings or not

        # By setting map_lat, map_lon, map_width and map_height to -1, OLGA automatically
        # tries to determine the best settings. Useful for setting up domains and first tests
        self.map_lat      = ([-1, -1]) # Central latitude of map [deg]
        self.map_lon      = ([-1, -1]) # Central longitude of map [deg]
        self.map_width    = ([-1, -1]) # Domain plot width [m]
        self.map_height   = ([-1, -1]) # Domain plot height [m]
        self.map_res      = (['l','i']) # Details of map (c=crude, l=low, i=interm, h=high)
        self.map_desc     = (['18x18km','6x6km']) 

        # ----------------------------------
        # Plot variables maps
        vars1 = (['pfd','swd','wstar','zidry','clouds','rr','wind10m','wind1000m'])         
        vars2 = (['pfd','wstar','zidry','cudepth'])         
        self.map_vars     = ([vars1,vars2]) # variables to plot per domain

        # ----------------------------------
        # Settings soundings (Detailed settings are in src/skewtlogp.py)
        self.sounding     = ([False,False]) # Make soundings or not
        sound_lat1        = ([self.map_lat[0]]) # Tmp arrays to populate locations, default=map center 
        sound_lon1        = ([self.map_lon[0]])  
        sound_lat2        = ([self.map_lat[1]]) 
        sound_lon2        = ([self.map_lon[1]]) 
        sound_name1       = (['center1'])
        sound_name2       = (['center2'])
        self.sound_name   = ([sound_name1,sound_name2]) # description of sounding location !! NO SPACES !! 
        self.sound_lat    = ([sound_lat1,sound_lat2]) # sounding latitudes 
        self.sound_lon    = ([sound_lon1,sound_lon2]) # sounding longitudes

        # Settings time series (meteograms)
        self.meteogr      = ([False,False]) # Make meteogram or not
        meteog_lat1       = ([self.map_lat[0]]) # Tmp arrays to populate locations, default=map center
        meteog_lon1       = ([self.map_lon[0]])  
        meteog_lat2       = ([self.map_lat[1]]) 
        meteog_lon2       = ([self.map_lon[1]])  
        meteog_name1      = (['center1'])
        meteog_name2      = (['center2'])
        self.meteogr_name = ([meteog_name1,meteog_name2]) # description of sounding location !! NO SPACES !! 
        self.meteogr_lat  = ([sound_lat1,sound_lat2]) # sounding latitudes 
        self.meteogr_lon  = ([sound_lon1,sound_lon2]) # sounding longitudes

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


