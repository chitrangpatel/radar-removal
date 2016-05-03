# radar-removal

A finely tuned algorithm to remove the Arecibo Punta Borinquen carsr radar operating at 
1275 and 1332 MHz.
This can also be tuned to remove any number of specified radar radio frequencies.

Basic idea:

1) This algorithm generates a time series of a specified bandwidth around specified radio 
   frequencies.
2) It feeds this time series to the bias monitor signal removal algorithm where it considers
   a time window of a specified length, computes the fft of that time window and identifies 
   bright periodic sub-pulse signals. It does this over the entire time series. Once the 
   brightest signals are identified for each window over the time series it picks out the 
   windows that are contaminated by the radar( > the specified threshold) after normalising 
   with respect to the median of these brightest signals.     
3) It outputs a text file containing all the range of time samples that are contaminated 
   by the radar. This file is used by our pipeline to mask out all the time samples identified 
   bu this algorithm over the entire observation bandwidth (not just the specified radar bandwidth). 

General usage:

python new_radar.py -f 1275.0 1332.0 -b 10.0 10.0 -w 0.25 0.25 -t 2.5 2.5 \\ 
                    -c 960 puppi_57505_J1905+0414_0167_subs_0001.fits 

-f : radar frequencies in MHz
-b : bandwidth in MHz to generate timeseries around each radar frequency
-w : The window length in seconds to use for calculating the fft and identifying 
     the samples corroupted by the radar. 
-t : Signal to noise threshold above which the time window is considered contaminated 
     by the radar.
-c : number of channels (960 for the PALFA survey) 
