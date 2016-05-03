#! /usr/bin/env python

import numpy as np
import argparse
import sys
import infodata
from time import strftime
from pypulsar.formats import psrfits
from pypulsar.formats import filterbank
import find_radar_mod
import rfifind
import subprocess
import glob

def channels_to_mask(data, nchannels, frequencytomask, bandwidth):
    """
    Finds out which channels are affected by radar based on the given
    radar radio frequency and the bandwidth of the radar. 
    Input: raw data object, number of channels (960 for PALFA),
           list [radar radio freqiencies], corresponding 
           list [bandwidths]
    Output: Int Array of channels to mask.        
    """
    lenchannel = (data.freqs.max()-data.freqs.min())/float(nchannels)
    channelstomask = np.array(())
    channel = int((frequencytomask-data.freqs.min())/lenchannel)
    if (bandwidth != 0.0):
        channelsperbandwidth = int(bandwidth/lenchannel)+1
        channelstomask = np.append(channelstomask, np.linspace(channel-channelsperbandwidth/2, 
                                   channel+channelsperbandwidth/2, channelsperbandwidth+1))
    else:
        channelstomask = np.append(channelstomask, np.linspace(channel-25, channel+25, 41))
    return np.unique(channelstomask.astype('int'))        

def chans_for_timeseries(channelstomask, nchannels):
    """
    Outputs an array of channels not containing the radar. 
    These channels will be masked out.
    """
    nchans = range(0, nchannels)
    for chan in channelstomask:
        nchans.remove(chan)   
    return np.array(nchans).astype('int32')

def write_mask_for_timeseries(old_mask, new_zap_chans, outbasenm):
    """
    Writes a mask file to get the timeseries of only the channels containing the radar.
    This way we can extract narrow band radar signals which are not very bright.
    """
    new_zap_ints = np.array(())
    new_zap_chans_per_int = np.array(())
    num_chans_per_int = np.array(())
    write_mask(old_mask, new_zap_chans, new_zap_ints, new_zap_chans_per_int, num_chans_per_int, outbasenm)
     
def make_timeseries(data, old_mask, frequenciestomask, bandwidth, nchannels, outbasenm):
    """
    Generates a timeseries of specific channels using Presto's prepdata.
    """
    for ii in range(len(frequenciestomask)):
        channelstomask = channels_to_mask(data, nchannels, 
                         float(frequenciestomask[ii]), float(bandwidth[ii]))
        new_zap_chans = chans_for_timeseries(channelstomask, nchannels) 
        write_mask_for_timeseries(old_mask, new_zap_chans, outbasenm+frequenciestomask[ii])
        subprocess.call(['cp', '%s_rfifind.bytemask'%outbasenm, '%s%s_new_rfifind.bytemask'%(outbasenm,frequenciestomask[ii])])
        subprocess.call(['cp', '%s_rfifind.ps'%outbasenm, '%s%s_new_rfifind.ps'%(outbasenm,frequenciestomask[ii])])
        subprocess.call(['cp', '%s_rfifind.inf'%outbasenm, '%s%s_new_rfifind.inf'%(outbasenm,frequenciestomask[ii])])
        subprocess.call(['cp', '%s_rfifind.rfi'%outbasenm, '%s%s_new_rfifind.rfi'%(outbasenm,frequenciestomask[ii])])
        subprocess.call(['cp', '%s_rfifind.stats'%outbasenm, '%s%s_new_rfifind.stats'%(outbasenm,frequenciestomask[ii])])
        subprocess.call(['prepdata', '-nobary', '-mask', '%s%s_new_rfifind.mask'%(outbasenm,frequenciestomask[ii]),
                          '-o', '%s%s'%(outbasenm,frequenciestomask[ii]), '%s.fits'%outbasenm])
        subprocess.call(['rm', '%s%s_new_rfifind.bytemask'%(outbasenm,frequenciestomask[ii])])
        subprocess.call(['rm', '%s%s_new_rfifind.ps'%(outbasenm,frequenciestomask[ii])])
        subprocess.call(['rm', '%s%s_new_rfifind.inf'%(outbasenm,frequenciestomask[ii])])
        subprocess.call(['rm', '%s%s_new_rfifind.rfi'%(outbasenm,frequenciestomask[ii])])
        subprocess.call(['rm', '%s%s_new_rfifind.stats'%(outbasenm,frequenciestomask[ii])])
    
def chans_per_int_with_radar(rawdatafile, inf, frequenciestomask, bandwidth, threshold, 
                    winlen, nchannels, start, old_mask, outbasenm):
    """
    Identifies the intervals contaminated by radar based on the original radar removal algorithm. 
    """
    masked_intervals = [] 
    for ii in range(len(frequenciestomask)):
        print 'radar frequency: %s MHz'%frequenciestomask[ii]
        rad_data = np.fromfile(outbasenm+'%s.dat'%frequenciestomask[ii], dtype = np.float32, count=-1)
        maxpows, medpows = find_radar_mod.compute_maxpows(rad_data, inf, float(winlen[ii]))
        rad_data, mask, blocks_to_mask = find_radar_mod.apply_mask_to_maxpows(rad_data, 
                                                maxpows, medpows, float(threshold[ii]))
        masked_intervals.append(find_radar_mod.write_radar_intervals(rad_data, outbasenm
                                                            +'%s'%frequenciestomask[ii]))
    return masked_intervals

def merge_intervals(masked_intervals, outbasenm):
    txtfiles = []
    clipbinsfiles = glob.glob("*radar_samples.txt")
    if clipbinsfiles:
        for ii in range(len(clipbinsfiles)):
            #txtfiles.append(np.loadtxt(outbasenm+'%s_radar_samples.txt'%frequenciestomask[ii], delimiter = ':'))
            txtfiles.append(np.loadtxt(clipbinsfiles[ii], delimiter = ':'))
    clipbinsfile = np.concatenate((txtfiles), axis=0)
    count = 0
    if len(clipbinsfile):
        order = np.lexsort(clipbinsfile.T)
        clipbinsfile = clipbinsfile[order]
        indices=[]
        for i in range(len(clipbinsfile)-1):
            if (clipbinsfile[i+1][0]==clipbinsfile[i][0]) and (clipbinsfile[i+1][1]==clipbinsfile[i][1]):
                indices.append(i)
        clipbinsfile = np.delete(clipbinsfile, np.asarray(indices), axis=0)
        for i in range(len(clipbinsfile)):
            count += clipbinsfile[i][1]-clipbinsfile[i][0]
        print 'number of masked intervals: %s'%len(clipbinsfile)
    with open(outbasenm+"_merged_radar_samples.txt",'w') as ff:
        ff.write("# Samples containing the radar\n")
        ff.write("# Intervals are samples to remove 'start:stop' (inclusive!)\n")
        ff.write("# First sample number is 0\n")
        ff.write("# Number of data samples to remove: %d of %d (%.2g %%)\n" %
                 (count, len(masked_intervals[1]), 100.0*count/len(masked_intervals[1])))
        np.savetxt(ff, clipbinsfile, "%d", delimiter=':')

def zap_chans(full_channels_to_mask, old_zap_chans):
    """
    Zaps entire channels.
    Input: array of channels that need to be completely zapped.
           array of channels already detected by rfifind.
    Output: array of input channels combined.
    """
    new_zap_chans = np.unique(np.append(old_zap_chans,full_channels_to_mask)).astype('int32') 
    return new_zap_chans

def zap_chans_per_int(old_zap_chans_per_int, ints_to_mask, mask_chans_per_int):
    """
    intervals to mask correspond to respective radar_chans
    """
    new_zap_chans_per_int = old_zap_chans_per_int
    num_chans_per_int = np.array((), dtype='int32')
    j=0
    for ii in ints_to_mask:
        new_zap_chans_per_int[ii] = np.unique(np.append(new_zap_chans_per_int[ii], 
                                              mask_chans_per_int[j])).astype('int32')
        j+=1
    for jj in new_zap_chans_per_int:
        num_chans_per_int = np.append(num_chans_per_int, len(jj)).astype('int32')
    return new_zap_chans_per_int, num_chans_per_int

def zap_ints(full_ints_to_mask, old_zap_ints):
    """
    Zaps entire intervals.
    Input: array of intervals that need to be completely zapped.
           array of intervals already detected by rfifind.
    Output: array of above intervals combined.
    """
    new_zap_ints = np.unique(np.append(old_zap_ints, full_ints_to_mask)).astype('int32')
    return new_zap_ints    
        
def read_mask(old_mask):
    """
    Reads in relevant information from an rfifind.mask file.
    """
    mask_zap_chans = np.array(list(old_mask.mask_zap_chans), dtype=np.int32)
    mask_zap_ints = old_mask.mask_zap_ints
    mask_zap_chans_per_int = old_mask.mask_zap_chans_per_int
    return mask_zap_chans, mask_zap_ints, mask_zap_chans_per_int

def write_mask(old_mask, new_zap_chans, new_zap_ints, new_zap_chans_per_int, 
               num_chans_per_int, outbasenm):
    """
    Writes a new .mask file in the same format as rfifind.mask.
    """
    info_array = np.array((old_mask.time_sig, old_mask.freq_sig, old_mask.MJD, old_mask.dtint,
                           old_mask.lofreq, old_mask.df), dtype='float64')
    with open(outbasenm+"_new_rfifind.mask",'wb') as f:
        info_array.tofile(f)
        np.array((old_mask.nchan, old_mask.nint, old_mask.ptsperint)).tofile(f)
        np.array((len(new_zap_chans)), dtype='int32').tofile(f)
        new_zap_chans.tofile(f)
        np.array((len(new_zap_ints)), dtype='int32').tofile(f)
        new_zap_ints.tofile(f)
        num_chans_per_int.tofile(f)
        for ii in new_zap_chans_per_int:
            ii.tofile(f)

def main():
    fn = args.infn
    print strftime("%Y-%m-%d %H:%M:%S")
    # Read in the raw data.(Probably not necessary anymore. Look into this)
    if fn.endswith(".fil"):
        filetype = "filterbank"
        rawdatafile = filterbank.filterbank(fn)
        if not args.outbasenm:
            outbasenm = fn[:-4]
        else:
            outbasenm = args.outbasenm
    elif fn.endswith(".fits"):
        filetype = "psrfits"
        rawdatafile = psrfits.PsrfitsFile(fn)
        if not args.outbasenm:
            outbasenm = fn[:-5]
        else:
            outbasenm = args.outbasenm
    else:
        raise ValueError("Cannot recognize data file type from\
                          extension. (Only '.fits' and '.fil'\
                          are supported.)")
    #Read data
    if filetype == 'psrfits':
        inffn = fn[:-5]+'_rfifind.inf'
    else:
        inffn = fn[:-4]+'.inf'
    inf = infodata.infodata(inffn)
    # Read in the original rfifind.mask information.
    old_mask = rfifind.rfifind(fn[:-5]+"_rfifind.mask")
    old_zap_chans, old_zap_ints, old_zap_chans_per_int = read_mask(old_mask)
    
    # Now make the timeseries only containing each of the radar signals individually.   
    make_timeseries(rawdatafile, old_mask, args.frequenciestomask, args.bandwidth, 
                    args.nchannels, outbasenm)
    start = 0
    ints_with_rfi_to_mask = []
    mask_zap_chans_per_int = []
    # Identify intervals contaminated by radar using the original radar removal algorithm.
    masked_intervals = chans_per_int_with_radar(rawdatafile,inf, 
                                      args.frequenciestomask, args.bandwidth, args.threshold, 
                                      args.winlen, args.nchannels, start, old_mask, outbasenm)
    merge_intervals(masked_intervals, outbasenm)
    print strftime("%Y-%m-%d %H:%M:%S")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Find radar in PALFA data.")
    parser.add_argument("-o", "--outbasenm", dest='outbasenm', default=None,
                        help="Base name of the output file containing " \
                             "list of intervals to clip. " \
                             "(Default: <inputfn>+'_radar_samples.txt'")
    parser.add_argument("-f", "--radar-frequencies", nargs = '+', dest='frequenciestomask', default=[],
                        help="List of radar frequencies that need to be masked. " \
                             "(Default: [])")
    parser.add_argument("-b", "--bandwidth", nargs = '+', dest='bandwidth', default=[],
                        help="Bandwidth of the radars in the same order as radar frequencies."\
                             "If you don't know the bandwidth of one of the signals put 0 for"\
                             "that entry and the default bandwidth will be selected."\
                             "Make sure that the order of radar frequencies match the"\
                             " associated bandwidths. Default: 40 channels wide.")
    parser.add_argument("-w", "--window-length", nargs = '+', dest='winlen', default=0.25,
                        help="Duration of window length in seconds for each radar frequency. " \
                             "(Default: 0.25 s)")
    parser.add_argument("-t", "--threshold", nargs = '+', dest='threshold', default=3.0,
                        help="Threshold of max power / median power " \
                             "above which an interval is considered " \
                             "to be contaminated by the radar. " \
                             "(Default: 3)")
    parser.add_argument("-c", "--nchannels", dest='nchannels', default=960, type=int,
                        help="Number of channels to divide the frequencies into. PALFA nchannels(960)" \
                             "(Default: 960)")
    parser.add_argument("infn",
                        help="PsrFits raw data file to use " \
                             "to identify intervals containing the radar.")
    args = parser.parse_args()
    if (len(args.frequenciestomask)!=len(args.bandwidth)) or (len(args.frequenciestomask)!=
                      len(args.winlen))or (len(args.frequenciestomask)!=len(args.threshold)):
        raise ValueError("The number of radar frequencies inputted must equal the "\
                         "number of bandwidths.")                                                      
    main()

