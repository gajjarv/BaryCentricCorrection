Post-detection barycentric correction for narrowband SETI. 

Code uses TEMPO to calculate expected Doppler velocity towards the source position from the input file header. 
It takes SIGPROC formatted filterbank file as an input and outputs barycentrically corrected SIGPROC file as well. 

Requirements:
-----------------
1. SIGPROC
2. TEMPO
3. gsl


Install:
------------------
1. Edit Makefile to connect SIGPROC and GSL
2. make clean
3. make

Use:
-----------------
> barycentre_seti "Input file" -verbose > "Output file"

Limitations:
-----------------
1. At this point, only works with 32-bit SIGPROC filterbank file ~~with positive foff~~
2. ~~Negative and positive relative velocity both fixes spectra-to-spectra~~
3. ~~Only negative relative velocity corrects the individual spectra~~
 
Testing:
-----------------
Testing for this code is inside directory "pipeline_tests"

Python notebook to create four filterbank files using setigen. 

nVelDec.fil : Negative Velocity which is decreasing spectra-to-spectra

nVelInc.fil : Negative Velocity which is increasing spectra-to-spectra

pVelDec.fil : Positive Velocity which is decreasing spectra-to-spectra

pVelInc.fil : Positive Velocity which is increasing spectra-to-spectra 

Ones these files are created using the python code, use following code to create four plots.

python3 barycentric_correction_code_check.py 

------------
## Post-detection Barycentric Correction

The Earth's movement around its axis and the Sun introduces a Doppler effect, which causes a shift in the frequency and arrival time of radio signals. To compensate for this relative motion and measure the precise emitted frequency and arrival time of radio signals, the barycentric correction is a crucial technique. For instance, pulsar timing requires barycentric correction to accurately compare observations taken at different epochs and telescopes. In the past, very few narrowband SETI surveys attempted to compare their observations across different epochs and telescopes, but for simultaneous observations, barycentric correction is imperative. The figure demonstrates the significant drift introduced by changes in the relative motion in the line of sight towards the beam centre for a length of 30 minutes of observations. Therefore, correction is necessary to compare detections from the two LOFAR stations.

![Narrowband drifting from Barycentric frequency](/path/to/figure)

Typically, barycentric corrections are introduced by adjusting the local oscillator at the time of observations. However, in this study, we record beamformed baseband voltages during observations and produce three different temporal and spectral resolution data products during post-processing. This is because we are interested in searching for a wide variety of signals, including narrowband signals, broadband transient signals, and wide-band pulsating signals. As introducing local oscillator shifts during observations can impact our other signal searches, we correct for barycentric drift after the channelization and detection of the baseband voltages for narrowband signal searches. We have developed a novel barycentric correction code specifically designed for high-spectral resolution SIGPROC filterbank products for technosignature searches. The code utilizes the TEMPO routine to calculate relative velocity towards the observing targets at both locations, thus allowing for precise correction of the barycentric drift. The code is publicly available for use [here](https://github.com/gajjarv/BaryCentricCorrection).

### Algorithm outline

The Doppler effect on a radio-emitting source can be stated as,

f_{em} = f_{obs} \left( 1 \pm \frac{v_{rel}}{c} \right).


Here, `f_{obs}` is the observed frequency, `f_{em}` is the emitted frequency, `v_{rel}` is the velocity of the source relative to the observer, and `c` is the speed of light. In this case, the Doppler effect only depends on the velocity of the source relative to the observer. For a source moving towards the observer we can consider `- v_{rel}`, while for a source moving away from the source we consider `+ v_{rel}`. For example, for the TIC 27677846 which is one of the target observed on July 15th, 2021 from both the LOFAR stations simultaneously. The expected relative velocity (`v_{rel}/c`) towards the source will be `-7.669 × 10^{-5}` and `-7.496 × 10^{-5}` causing relative shift (`f_{em} - f_{obs}`) of 11.504 kHz and 11.246 kHz for a transmitted signal at 150 MHz observed at the Sweden and Ireland stations, respectively. These are significant shifts and distinct at the two stations which needs to be corrected in order to compare this same signal observed at the two stations.

