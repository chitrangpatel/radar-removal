#! /usr/bin/env python
import sys
import os.path
import argparse
import warnings

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import infodata

DTYPE = 'float32'

def compute_maxpows(data, inf, winlen):
    nblocks = int(np.ceil(inf.N*inf.dt/winlen))
    nbin = inf.N/float(nblocks)

    # Compute maximum Fourier power for each block of time series data
    maxpows = np.zeros(nblocks)
    medpows = np.zeros(nblocks)
    for iblock in xrange(nblocks):
        istart = int(np.round(iblock*nbin))
        iend = int(np.round((iblock+1)*nbin))
        powers = np.abs(np.fft.rfft(data[istart:iend]))**2
        maxpows[iblock] = np.max(powers[1:])
        medpows[iblock] = np.median(powers[1:])
        imax = np.argmax(powers[1:])
    maxpows[-1] = np.median(maxpows[:-1])
    # Scale the maximum powers
    maxpows /= np.median(maxpows)
    #maxpows /=np.median(powers)
    return maxpows, medpows

def compute_minpows(data, inf, winlen):
    nblocks = int(np.ceil(inf.N*inf.dt/winlen))
    nbin = inf.N/float(nblocks)

    # Compute maximum Fourier power for each block of time series data
    minpows = np.zeros(nblocks)
    medpows = np.zeros(nblocks)
    variance = np.zeros(nblocks)
    for iblock in xrange(nblocks):
        istart = int(np.round(iblock*nbin))
        iend = int(np.round((iblock+1)*nbin))
        #powers = np.abs(np.fft.rfft(data[istart:iend]))**2
        powers = np.abs(data[istart:iend])
        minpows[iblock] = np.min(powers[1:])
        medpows[iblock] = np.median(powers[1:])
        variance[iblock] = np.var(powers[1:])
        imin = np.argmin(powers[1:])
    minpows[-1] = np.median(minpows[:-1]) 
    minpows -= medpows
    # Scale the minimum powers
    minpows /= np.median(minpows)
    variance /= np.median(variance)
    print minpows[:25]
    print variance[:25]
    #return minpows, medpows
    return variance, medpows

def apply_mask_to_maxpows(data, maxpows, medpows, thresh):
    nblocks = len(maxpows)
    nbin = len(data)/float(nblocks)

    # Generate a mask based on block-wise Fourier power
    # being larger than a threshold
    mask = np.zeros(len(data), dtype=bool)
    block_to_mask = []
    for iblock in xrange(len(maxpows)):
        istart = int(np.round(iblock*nbin))
        iend = int(np.round((iblock+1)*nbin))
        if maxpows[iblock] > thresh:
            mask[istart:iend] = 1
            block_to_mask.append(iblock)
    masked = np.ma.masked_array(data, mask=mask, fill_value=np.median(data))
    return masked, mask, block_to_mask

def apply_mask_to_minpows(data, minpows, medpows, thresh):
    nblocks = len(minpows)
    nbin = len(data)/float(nblocks)

    # Generate a mask based on block-wise Fourier power
    # being larger than a threshold
    mask = np.zeros(len(data), dtype=bool)
    block_to_mask = []
    for iblock in xrange(len(minpows)):
        istart = int(np.round(iblock*nbin))
        iend = int(np.round((iblock+1)*nbin))
        if minpows[iblock] > thresh:
            mask[istart:iend] = 1
            block_to_mask.append(iblock)
    masked = np.ma.masked_array(data, mask=mask, fill_value=np.median(data))
    return masked, mask, block_to_mask


def write_masked_dm0_timeseries(masked, outbasenm):
    maskedfn = outbasenm+"_noradar.dat"
    (masked.filled(fill_value=np.ma.median(masked))).tofile(maskedfn)
    # Create soft link for inf file
    inflink = maskedfn[:-4]+".inf"
    if not os.path.exists(inflink):
        os.symlink(outbasenm+".inf", inflink)


def write_radar_intervals(masked, outbasenm):
    # Invert the mask
    inverted = invert_mask(masked) 
    slices = np.ma.flatnotmasked_contiguous(inverted)
    
    if slices:
        towrite = np.empty((len(slices), 2))
        for ii, badslice in enumerate(slices):
            towrite[ii] = (badslice.start-1000, badslice.stop+1000)
    #radar_samples = np.flatnonzero(np.ma.getmaskarray(data))
    with open(outbasenm+"_radar_samples.txt", 'w') as ff:
        ff.write("# Samples containing the radar\n")
        ff.write("# Intervals are samples to remove 'start:stop' (inclusive!)\n")
        ff.write("# First sample number is 0\n")
        ff.write("# Number of data samples to remove: %d of %d (%.2g %%)\n" %
                 (inverted.count(), len(inverted), 100.0*inverted.count()/len(inverted)))
        if slices:
            np.savetxt(ff, towrite, "%d", delimiter=':')
    return inverted

def invert_mask(masked):
    return np.ma.masked_array(masked.data, mask=~masked.mask)


def validate_timeseries(inf):
    if inf.bary:
        raise ValueError("Radar removal should only be performed " \
                         "using topocentric time series as input!")
    if inf.DM != 0:
        raise ValueError("Radar removal should only be performed " \
                         "using a DM=0 pc/cc time series as input!")


def plot_timeseries_comparison(masked, inf):
    fig = plt.figure(figsize=(16,4))

    times = np.arange(len(masked))*inf.dt
    warnings.warn("Only plotting every 10th point of time series.")
    plt.plot(times[::10], masked.data[::10], 'k-', drawstyle='steps-post', 
             label='Time series', zorder=1)

    inverted = invert_mask(masked) 
    slices = np.ma.flatnotmasked_contiguous(inverted)
    if slices:
        for ii, badslice in enumerate(slices):
            if ii == 0:
                label='Radar indentified'
            else:
                label="_nolabel"
            tstart = inf.dt*(badslice.start)
            tstop = inf.dt*(badslice.stop-1)
            plt.axvspan(tstart, tstop, alpha=0.5, 
                        fc='r', ec='none', zorder=0, label=label)
    
    plt.figtext(0.02, 0.02,
                "Frac. of data masked: %.2f %%" % ((len(masked)-masked.count())/float(len(masked))*100), 
                size='x-small')
    plt.figtext(0.02, 0.05, inf.basenm, size='x-small')

    plt.xlabel("Time (s)")
    plt.ylabel("Intensity")
    plt.xlim(0, times.max()+inf.dt)

    plt.subplots_adjust(bottom=0.15, left=0.075, right=0.98)


def plot_powerspec_comparison(masked, inf):
    fig = plt.figure(figsize=(16,4))

    # Make diagnostic plots
    freqs = np.fft.fftfreq(int(inf.N), inf.dt)
    nfft = int(inf.N/2)
    freqs = freqs[:nfft]
    plt.plot(freqs, np.abs(np.fft.rfft(masked.data)[:nfft])**2, 'r-')
    plt.plot(freqs, np.abs(np.fft.rfft(masked.filled(fill_value=np.ma.median(masked)))[:nfft])**2, 'k-')
    plt.xlabel("Freq (Hz)")
    plt.ylabel("Raw Power")
    plt.xscale('log')
    plt.yscale('log')
    plt.figtext(0.02, 0.02,
                "Frac. of data masked: %.2f %%" % ((len(masked)-masked.count())/float(len(masked))*100), 
                size='x-small')
    plt.figtext(0.02, 0.05, inf.basenm, size='x-small')

    plt.xlim(0.066, 2000)

    plt.subplots_adjust(bottom=0.15, left=0.075, right=0.98)


def read_datfile(fn):
    # Get input dat/inf file names
    if not fn.endswith(".dat"):
        raise ValueError("Input file must be a PRESTO '.dat' file!")
    inffn = fn[:-4]+'.inf'
    
    # Read inf/dat files
    inf = infodata.infodata(inffn) 
    rawdata = np.fromfile(fn, dtype=DTYPE, count=inf.N) 
    #validate_timeseries(inf)

    return rawdata, inf


def main():
    fn = args.infn
    if args.outbasenm is None:
        outbasenm = fn[:-4] # Trim off '.dat'
    else:
        outbasenm = args.outbasenm

    rawdata, inf = read_datfile(fn)

    maxpows, medpows = compute_maxpows(rawdata, inf, args.winlen)
    data, mask, blocks_to_mask = apply_mask_to_maxpows(rawdata, maxpows, medpows, args.maxthresh)
    
    minpows, medpows = compute_minpows(data, inf, args.winlen)
    data, mask, blocks_to_mask = apply_mask_to_minpows(data, minpows, medpows, args.minthresh)
    print minpows
    print "Frac. of data masked: %.2f %%" % ((len(data)-data.count())/float(len(data))*100), 

    # Write out masked file
    #write_masked_dm0_timeseries(data, outbasenm)

    # Write out list of masked samples
    write_radar_intervals(data, outbasenm)

    if args.plot_timeseries:
        plot_timeseries_comparison(data, inf)
        try:
            plt.savefig(outbasenm+"_timeseries_compare.png")
        except:
            pass

    if args.plot_powerspec:
        plot_powerspec_comparison(data, inf)
        try:
            plt.savefig(outbasenm+"_powerspec_compare.png")
        except:
            pass    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Find radar in PALFA data.")
    parser.add_argument("-o", "--outbasenm", dest='outbasenm', default=None,
                        help="Base name of the output file containing " \
                             "list of intervals to clip. " \
                             "(Default: <inputfn>+'_radar_samples.txt'")
    parser.add_argument("-w", "--window-length", dest='winlen', default=0.25,
                        type=float,
                        help="Duration of window length in seconds. " \
                             "(Default: 0.25 s)")
    parser.add_argument("-t", "--minthreshold", dest='minthresh', default=1.5,
                        type=float,
                        help="Threshold of (min power -median power) / median power " \
                             "above which an interval is considered " \
                             "to be contaminated by the dips in the time series. " \
                             "(Default: 1.5)")
    parser.add_argument("-T", "--maxthreshold", dest='maxthresh', default=3,
                        type=float,
                        help="Threshold of max power / median power " \
                             "above which an interval is considered " \
                             "to be contaminated by the radar. " \
                             "(Default: 3)")
    parser.add_argument("--no-timeseries-plot", dest="plot_timeseries", action="store_false",
                        help="Do not plot time series comparison. " \
                             "(Default: create time series plot)")
    parser.add_argument("--no-powerspec-plot", dest="plot_powerspec", action="store_false",
                        help="Do not plot power spectrum comparison. " \
                             "(Default: create power spectrum plot)")
    parser.add_argument("infn",
                        help="Topocentric time series at DM=0 pc/cc to use " \
                             "to identify intervals containing the radar.")
    args = parser.parse_args()
    main()
