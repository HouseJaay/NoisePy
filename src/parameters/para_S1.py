import os
import pandas as pd

# absolute path parameters
rootpath  = '/home/shijie/data/SL_01/'                               # root path for this data processing
CCFDIR    = os.path.join(rootpath,'CCF')                                    # dir to store CC data
DATADIR   = os.path.join(rootpath,'asdf_100Hz')                               # dir where noise data is located
local_data_path = os.path.join(rootpath,'2016_*')                           # absolute dir where SAC files are stored: this para is VERY IMPORTANT and has to be RIGHT if input_fmt is not h5 for asdf!!!
locations = os.path.join(DATADIR,'metadata', 'station.csv')                             # station info including network,station,channel,latitude,longitude,elevation: only needed when input_fmt is not h5 for asdf

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
# read station list
if input_fmt != 'h5':
    if not os.path.isfile(locations):
        raise ValueError('Abort! station info is needed for this script')
    locs = pd.read_csv(locations)

# pre-processing parameters
cc_len    = 200                                                            # basic unit of data length for fft (sec)
step      = 50                                                             # overlapping between each cc_len (sec)
smooth_N  = 10                                                              # moving window length for time/freq domain normalization if selected (points)

# cross-correlation parameters
maxlag         = 1.5                                                        # lags of cross-correlation to save (sec)
substack       = True                                                       # True = smaller stacks within the time chunk. False: it will stack over inc_hours
                                                                            # for instance: substack=True, substack_len=cc_len means that you keep ALL of the correlations
                                                                            # if substack=True, substack_len=2*cc_len, then you pre-stack every 2 correlation windows.
substack_len   = cc_len                                                     # how long to stack over (for monitoring purpose): need to be multiples of cc_len
smoothspect_N  = 10                                                         # moving window length to smooth spectrum amplitude (points)

# criteria for data selection
max_over_std = 10                                                           # threahold to remove window of bad signals: set it to 10*9 if prefer not to remove them

# maximum memory allowed per core in GB
MAX_MEM = 4.0

# load useful download info if start from ASDF
if input_fmt == 'h5':
    dfile = os.path.join(DATADIR,'download_info.txt')
    down_info = eval(open(dfile).read())
    samp_freq = down_info['samp_freq']
    freqmin   = down_info['freqmin']
    freqmax   = down_info['freqmax']
    start_date = down_info['start_date']
    end_date   = down_info['end_date']
    inc_hours  = down_info['inc_hours']
    ncomp      = down_info['ncomp']
else:   # sac or mseed format
    samp_freq = 20
    freqmin   = 0.05
    freqmax   = 2
    start_date = ["2016_07_01_0_0_0"]
    end_date   = ["2016_07_02_0_0_0"]
    inc_hours  = 24
dt = 1/samp_freq