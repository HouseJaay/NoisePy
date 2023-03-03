import os
from obspy import UTCDateTime

# absolute path parameters
rootpath  = '/scratch/goxu/hao_shijie/huidong/phase1/'                               # root path for this data processing
#rootpath = '/media/shijie/NTU-HSJ-01/data/hawaii_new/'
CCFDIR    = os.path.join(rootpath,'CCF')                             # dir to store CC data
DATADIR   = os.path.join(rootpath,'phase1-h5-50Hz')                      # dir where noise data is located
STACKDIR = os.path.join(rootpath,'STACK')

# some control parameters
input_fmt   = 'h5'                                                          # string: 'h5', 'sac','mseed'
freq_norm   = 'rma'                                                         # 'no' for no whitening, or 'rma' for running-mean average, 'phase_only' for sign-bit normalization in freq domain.
time_norm   = 'rma'                                                          # 'no' for no normalization, or 'rma', 'one_bit' for normalization in time domain
cc_method   = 'xcorr'                                                       # 'xcorr' for pure cross correlation, 'deconv' for deconvolution; FOR "COHERENCY" PLEASE set freq_norm to "rma", time_norm to "no" and cc_method to "xcorr"
flag        = False                                                          # print intermediate variables and computing time for debugging purpose
acorr_only  = False                                                         # only perform auto-correlation
xcorr_only  = True                                                          # only perform cross-correlation or not
ncomp       = 1                                                             # 1 or 3 component data (needed to decide whether do rotation)

# station/instrument info for input_fmt=='sac' or 'mseed'
stationxml = False                                                          # station.XML file used to remove instrument response for SAC/miniseed data
rm_resp   = 'no'                                                            # select 'no' to not remove response and use 'inv','spectrum','RESP', or 'polozeros' to remove response
respdir   = os.path.join(rootpath,'resp')                                   # directory where resp files are located (required if rm_resp is neither 'no' nor 'inv')

# pre-processing parameters
cc_len    = 1*3600                                                            # basic unit of data length for fft (sec)
step      = 0.5*3600                                                             # overlapping between each cc_len (sec)
smooth_N  = 10                                                              # moving window length for time/freq domain normalization if selected (points)

# cross-correlation parameters
maxlag         = 6                                                        # lags of cross-correlation to save (sec)
substack       = True                                                       # True = smaller stacks within the time chunk. False: it will stack over inc_hours
                                                                            # for instance: substack=True, substack_len=cc_len means that you keep ALL of the correlations
                                                                            # if substack=True, substack_len=2*cc_len, then you pre-stack every 2 correlation windows.
substack_len   = cc_len                                                     # how long to stack over (for monitoring purpose): need to be multiples of cc_len
smoothspect_N  = 10                                                         # moving window length to smooth spectrum amplitude (points)
stack_method = 'rms-select' # rms-select or linear

# criteria for data selection
max_over_std = 10                                                           # threahold to remove window of bad signals: set it to 10*9 if prefer not to remove them

# maximum memory allowed per core in GB
MAX_MEM = 4.0

samp_freq = 50
dt = 1/samp_freq

freqmin, freqmax = 0.1, 5
loctype = "lonlat"

# Parameters for rms selection stacking
vmin, vmax = 1.2, 4.0 # used to calculate snr
maxperiod = 10
