import sys
from obspy.geodetics.base import gps2dist_azimuth
import pyasdf
import random
import logging
import os, glob
from os.path import join, exists
import numpy as np
import noise_module
from mpi4py import MPI
import importlib
from obspy import UTCDateTime


try:
    proj_path = sys.argv[1]
    proj_name = sys.argv[2]
except IndexError:
    print("Error, Need to provide project path")
    print("For example, noisepy_parameters.py is in /scratch/shijie001/data/hawaii/ccf1/")
    print("1st parameter: /scratch/shijie001/data/hawaii/")
    print("2cn parameter: ccf1")
    exit()

sys.path.append(proj_path)
par = importlib.import_module(proj_name + '.noisepy_parameters')

if par.stack_method == 'rms-select':
    par.substack = True
    par.substack_len = par.cc_len

def cut_trace(tr, starttime, endtime, cc_len, step, sps):
    """
    st: continuous noise data of one station, Trace format
    starttime: UTCDateTime
    enttime: UTCDateTime
    cc_len: seconds, length of one segment
    """
    nseg = int((endtime - starttime) / step)
    npts_seg = int(cc_len * sps)
    data = np.zeros(shape=(nseg, npts_seg))
    data_std = np.zeros(shape=nseg)
    data_t = np.zeros(shape=nseg)
    all_std = np.std(tr.data)
    for iseg in range(nseg):
        t0 = starttime + iseg*step
        data_t[iseg] = t0
        if t0 < tr.stats.starttime:
            continue
        elif t0 + cc_len > tr.stats.endtime:
            continue
        n0 = int((t0-tr.stats.starttime)*sps)
        data[iseg, :] = tr.data[n0:n0+npts_seg]
        data_std[iseg] = np.max(np.abs(data[iseg])) / all_std

    return data_std, data, data_t

#############################################
# Import parameters from parameter file
# Put parameter in the working directory

ccfdir = par.CCFDIR
stackdir = par.STACKDIR
input_fmt = par.input_fmt
flag = par.flag
DATADIR = par.DATADIR
MAX_MEM = par.MAX_MEM
ncomp = par.ncomp
cc_len = par.cc_len
step = par.step
samp_freq = par.samp_freq
acorr_only = par.acorr_only
max_over_std = par.max_over_std
dt = par.dt
vmin, vmax = par.vmin, par.vmax
maxperiod = par.maxperiod
loctype = par.loctype
maxlag = par.maxlag

fc_para = {
    "time_norm": par.time_norm,
    "freq_norm": par.freq_norm,
    "smooth_N": par.smooth_N,
    "dt": par.dt,
    "freqmin": par.freqmin,
    "freqmax": par.freqmax,
    "maxlag": par.maxlag,
    "cc_method": par.cc_method,
    "cc_len": par.cc_len,
    "substack": par.substack,
    "substack_len": par.substack_len,
    "smoothspect_N": par.smoothspect_N
}

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    if not os.path.isdir(ccfdir):
        os.mkdir(ccfdir)
    if not os.path.isdir(stackdir):
        os.mkdir(stackdir)
    files = sorted(glob.glob(os.path.join(DATADIR, '*.h5')))
    pairs = []
    for i in range(0, len(files)-1):
        for j in range(i+1, len(files)):
            pairs.append((files[i], files[j]))
    # optimize for IO
    random.shuffle(pairs)
else:
    pairs = None
pairs = comm.bcast(pairs, root=0)
npairs = len(pairs)
logging.basicConfig(filename=join(ccfdir, "node_%d.log"%rank), level=logging.WARNING)

for ick in range(rank, npairs, size):
    f1, f2 = pairs[ick][0], pairs[ick][1]
    logging.info("Start %s %s" % (f1, f2))
    ds1 = pyasdf.ASDFDataSet(f1, mode='r', mpi=False)
    ds2 = pyasdf.ASDFDataSet(f2, mode='r', mpi=False)
    sta1 = ds1.waveforms.list()[0]
    sta2 = ds2.waveforms.list()[0]
    tags1 = set(ds1.waveforms[sta1].get_waveform_tags())
    tags2 = set(ds2.waveforms[sta2].get_waveform_tags())
    tags = list(tags1.intersection(tags2))
    s_corr_all = []
    for tag in tags:
        t_start = UTCDateTime(*tuple(map(int, tag.split('_'))))
        t_end = t_start + par.data_len
        st1 = ds1.waveforms[sta1][tag]
        st2 = ds2.waveforms[sta2][tag]
        st1.merge(method=1, fill_value=0)
        st2.merge(method=1, fill_value=0)
        data_std1, data1, data_t = cut_trace(st1[0], t_start, t_end, cc_len, step, samp_freq)
        data_std2, data2, _ = cut_trace(st2[0], t_start, t_end, cc_len, step, samp_freq)
        

        white1 = noise_module.noise_processing(fc_para, data1)
        white1 = np.conjugate(white1)
        Nfft = white1.shape[1]
        Nfft2 = Nfft//2
        ind1 = np.where((data_std1<max_over_std)&(data_std1>0))[0]

        white2 = noise_module.noise_processing(fc_para, data2)
        ind2 = np.where((data_std2 < max_over_std) & (data_std2 > 0))[0]
        ind = np.intersect1d(ind1, ind2)
        if len(ind) < 1:
            print("Don't have enough data %s, %s" % (sta1, sta2))
            continue

        s_corr, t_corr, n_corr = noise_module.correlate(
            white1[ind, :Nfft2], white2[ind, :Nfft2], fc_para, Nfft, data_t)
        s_corr_all.append(s_corr)
    if not s_corr_all:
        continue
    
    t = np.arange(-Nfft2+1, Nfft2)*dt
    ind = np.where(np.abs(t) <= maxlag)[0]
    npts_corr = len(ind)
    s_corr_all = np.concatenate(s_corr_all, axis=0)
    s_corr_all = s_corr_all.reshape((-1, npts_corr))
    
    if loctype == "lonlat":
        cord1 = ds1.get_all_coordinates()[sta1]
        cord2 = ds2.get_all_coordinates()[sta2]
        lat1, lon1 = cord1['latitude'], cord1['longitude']
        lat2, lon2 = cord2['latitude'], cord2['longitude']
        dist, _, __ = gps2dist_azimuth(lat1, lon1, lat2, lon2)
        dist /= 1000

    if par.stack_method == 'rms-select':
        idx1 = int(dist / vmax * (1 / dt))
        idx2 = int(dist / vmin * (1 / dt))
        idx3 = int(idx2 + maxperiod * 4 / dt)
        stacked, nstack = noise_module.rms_selection_stack(s_corr_all, idx1, idx2, idx3, 1)
        data_type = "Allstack_rmss"
    elif par.stack_method == 'linear':
        stacked = np.sum(s_corr_all, axis=0)
        nstack = s_corr_all.shape[0]
        data_type = "Allstack_linear"

    # output
    cdir = join(stackdir, sta1)
    fo = join(cdir, "%s_%s.h5" % (sta1, sta2))
    if not exists(cdir):
        os.mkdir(cdir)
    dso = pyasdf.ASDFDataSet(fo, mode='w', mpi=False)
    stack_parameters = {'ngood': nstack, 'dist': dist, 'dt': dt, 'maxlag': maxlag,
                        'lonS':lon1, 'latS':lat1, 'lonR':lon2, 'latR':lat2}
    dso.add_auxiliary_data(data=stacked, data_type=data_type, path='ZZ', parameters=stack_parameters)

comm.barrier()
if rank == 0:
    sys.exit()