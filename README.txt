# PhaseAnalysis Toolbox 

General description - phase analysis matlab toolbox was composed to investigate 
spike-phase preference in different neuronal populations of the hippocampus 
and brainstem during waking behavior and subsequent sleep. This program reads 
neuralynx CSC files and Plexon Spike Sorted files into matlab. 

Dependent functions - to run this program, there are required subfunctions within
the subfunctions folder that allow the proper read and structuring of the data.
Loading files: stateLetter2NumberConverter, Nlx2MatCSC, ts, tsd, replay_LoadSpikes
Hilbert for theta: InstSig_theta, InstSig
Phase statistics: circstat toolbox (circ_mean;circ_rtest;circ_r)

 
 