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

# Post-detection Barycentric Correction

The movement of the Earth around its axis and the Sun introduces a Doppler effect that causes radio signals' frequency and arrival time to shift. In the case of a radio-emitting source, this effect can be described by the following equation:

$$f_{obs} = f_{em} \left( 1 - v_{rel} \right),$$

where $f_{obs}$ is the observed frequency, $f_{em}$ is the emitted frequency (or barycentric frequency), and $v_{rel}$ is the velocity of the source relative to the observer, normalized by the speed of light. The Doppler effect only depends on the velocity of the source relative to the observer, and for a source moving towards the observer, we can consider $-v_{rel}$, while for a source moving away from the observer, we consider $+v_{rel}$. The relative velocity between the observer and the source will differ for different observing epochs, the location of the source in the sky, and the geographical location of the telescope.

For instance, the target TIC\,27677846 was observed on July 15th, 2021, from both the LOFAR stations simultaneously. The expected relative velocity ($v_{rel}$) towards the source was $-7.669\times10^{-5}$ and $-7.496\times10^{-5}$, causing a relative shift ($f_{em}-f_{obs}$) of +11.504\,kHz and +11.246\,kHz for a hypothetical ETI signal transmitted at a constant frequency of 150\,MHz observed at the Sweden and Ireland stations, respectively. These are significant shifts that are distinct at the two stations, and they need to be corrected to compare the same signal observed at the two stations.

![Barycentric Frequency Example](https://github.com/gajjarv/BaryCentricCorrection/blob/master/pipeline_tests/Barycentric_frequency_drift.png)

The plot depicts the Doppler drift of a narrowband signal in the topocentric observing frame at four different observing epochs. Simulated waterfalls with narrowband signals observed from the Irish LOFAR station towards the direction of TIC\,27677846 are shown for different times of the year in blue. The expected sign and direction of change of the relative velocity are labeled at the top of each plot. It is assumed that a hypothetical narrowband ETI signal is transmitted at a constant frequency of 1420\,MHz (with zero drift rate). As shown, the same signal is observed at different frequencies and drift rates depending on the sign and direction of change of the relative velocity at different epochs of observations. For instance, in the first panel, the relative velocity is positive and increases with observing time. Therefore, the observing frequency has been shifted to a lower frequency (as described by Equation \ref{eq:doppler}), and it continues to shift to even lower frequencies with time.

In order to accurately measure the precise emitted frequency and arrival time of radio signals and compensate for the relative motion caused by the Earth's movement, barycentric correction is a crucial technique. Pulsar timing, for example, requires barycentric correction to accurately compare observations taken at different epochs and telescopes. Although few narrowband SETI surveys in the past attempted to compare their observations across different epochs and telescopes, barycentric correction is imperative for simultaneous observations.

Typically, barycentric corrections are introduced by adjusting the local oscillator during observations. However, in our study, we record beamformed baseband voltages during observations and produce three different data products with varying temporal and spectral resolutions during post-processing. We are interested in searching for a wide variety of signals, including narrowband signals, broadband transient signals, and wide-band pulsating signals. Introducing local oscillator shifts during observations can impact our other signal searches. Therefore, we correct for barycentric drift after the channelization and detection of the baseband voltages for narrowband signal searches.

The movement of the Earth around its axis and the Sun introduces a Doppler effect that causes radio signals' frequency and arrival time to shift. Figure \ref{fig:bary_simulated_drift} depicts the Doppler drift of a narrowband signal in the topocentric observing frame at four different observing epochs. Simulated waterfalls with narrowband signals observed from the Irish LOFAR station towards the direction of TIC\,27677846 are shown for different times of the year in blue. The expected sign and direction of change of the relative velocity are labeled at the top of each plot. It is assumed that a hypothetical narrowband ETI signal is transmitted at a constant frequency of 1420\,MHz (with zero drift rate). As shown, the same signal is observed at different frequencies and drift rates depending on the sign and direction of change of the relative velocity at different epochs of observations. For instance, in the first panel, the relative velocity is positive and increases with observing time. Therefore, the observing frequency has been shifted to a lower frequency (as described by Equation \ref{eq:doppler}), and it continues to shift to even lower frequencies with time.

We have developed a novel barycentric correction code specifically designed for high-spectral resolution SIGPROC filterbank products for technosignature searches. The code uses the TEMPO routine to calculate relative velocity towards the observing targets at both locations, thus allowing for precise correction of the barycentric drift.The Doppler shift caused by the relative motion between the transmitter and receiver will change over time as it is observed. Over a longer time frame, these shifts will exhibit a sinusoidal curve with a sidereal year period. Over a shorter time frame, the same pattern (superimposed on the yearly pattern) will be visible, but with a sidereal day period. Consequently, the relative velocity will change during the observation period, thereby altering the observed frequency of the received narrowband signal. This leads to a drift in the narrowband signal observed by the observer. Figure \ref{fig:bary_simulated_drift} displays examples of observed drifts at four different epochs for the same narrowband signal source observed from the same location. It is evident that if the relative velocity is positive and increasing with time (leftmost plot in Figure \ref{fig:bary_simulated_drift}), the signal, which is stationary in the barycentric frame, will drift towards lower frequencies as time progresses in the topocentric frame, as per Equation \ref{eq:doppler}. For our case, we assume that the frequency channels in a given filterbank file are ordered in descending order, with the highest frequency channel ($f_{1}$) for each time sample being the first channel.

The goal of our tool is to shift every frequency channel from the observing frame to the actual emitted frequency frame after correcting for the barycentric relative velocity. This is done to remove any additional narrowband signal drift introduced by the relative velocity. We aim to keep the first channel frequency of all time samples the same in the barycentric frame, thus relative shifts between spectra are needed to apply for each time sample corresponding to the inferred relative velocities.

Our tool measures the relative velocity at each time sample towards a given direction in the sky from a given telescope at the time of observations. Let's assume an observing scenario where $v_{rel}>0$. Following Equation (1), we can state that the emitted frequency (or barycentric frequency) will be higher than the observed frequency. That means that observations stored in the first topocentric frequency channel correspond to a higher barycentric frequency. Thus, this spectrum needs to be shifted to a higher frequency. If the relative velocity increases with time, the consecutive time sample's first channel barycentric frequency will be slightly higher than the previous sample's first channel emitted frequency. Our tool thus shifts the spectra of each time sample towards higher frequency such that the first channel's emitted frequency matches across all time samples. Due to these shifts, we either replace empty channels at the edge of the spectra with zeros in case of squeeze or drop extra channels in case of expansion.

As given in Equation \ref{eq:doppler}, the Doppler frequency shifts are frequency-dependent, impacting higher frequencies more relative to lower frequencies. In other words, for spectra where frequencies are ordered from higher to lower frequency, the first channel will be shifted relatively more compared to the last frequency channel. To compensate for this, we either expand or squeeze spectra as shown in Figure \ref{fig:spectra_expands_squeeze}. Figure \ref{fig:code_outline} outlines the logical flow of the code for a case of an input filterbank file with a descending order of frequency. By comparing Figure \ref{fig:bary_simulated_drift} and the code outline in Figure \ref{fig:code_outline}, we can consider one of the cases where the relative velocity is negative and increasing in absolute value with time. In this case, we need to shift consecutive spectra to lower and lower frequencies to match up their first frequency channel. Furthermore, for any negative relative velocity (either increasing or decreasing with time), we need to squeeze the individual spectra as shown in Figure \ref{fig:spectra_expands_squeeze}. Similarly, the same can be consider for the case of positive relative velocity. 

# Pipeline outline

![Pipeline output](https://github.com/gajjarv/BaryCentricCorrection/blob/master/pipeline_tests/code_outline.png)

An outline of the post-detection barycentric correction algorithm for an input filterbank file in SIGPROC format. For this case, the input filterbank file has a descending order in frequency, and $v_{rel}$ represents the relative velocity between the transmitter and observer. The algorithm considers two cases depending on whether the relative velocity is positive or negative, which indicates whether the source is moving away from or towards the observer, respectively. Each of these cases is further divided into two where the absolute value of the relative velocity can either increase or decrease, requiring the spectra to be shifted to either the higher or lower frequency end. For all cases with $+v_{rel}$, each spectra is expanded, and for $-v_{rel}$, each spectra is squeezed. The code then writes each of these spectra into another SIGPROC filterbank file, which will have each channel frequency closely corrected to the barycentric frame of reference.
