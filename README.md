# IIA_SF1_Data_Analysis


## Description
This project introduces signal modelling/processing techniques and applies them to audio and musical signals. 
Firstly, non-parametric methods based on transforms such as the Discrete Fourier Transform (DFT) and its fast variant the FFT are studied and experiments are carried out with windowing, frequency resolution and noise reduction. 
Furthermore, I have implemented various probabilistic inference methods and explored a wide variety of parametric models with estimation using least squares, maximum likelihood and Bayesian techniques. I analysed the constant model, the linear trend model, the autoregressive (AR) model and the sinusoidal model in detail and studied the model choice within likelihood and Bayesian probabilistic settings as well. I have then successfully applied these techniques to audio signals, using models to perform packet loss concealment and interpolation of missing data in audio, as well as constrained interpolation for clipped and/or heavily quantised signals. During the four-week project, I have solidified our understanding of signal analysis using filters and windowing functions, explored noise reduction methods such as Wiener filter and power subtraction filters, and mastered the probabilistic inference methods to perform lost audio recovery and time-series audio prediction.


## Results
Examples of tasks accomplished are demonstrated below:

### 1. Amplitude Modulation
The windowed spectrum on audio files with amplitude modulation is shown in the figure below:
(a) An + B:
![alt text](https://github.com/YoujingYu99/IIA_SF1_Data_Analysis/blob/master/plots/an1.jpg?raw=true)

(b) 1 + β sin(φ n):
![alt text](https://github.com/YoujingYu99/IIA_SF1_Data_Analysis/blob/master/plots/an2.jpg?raw=true)

(c) a_{n-1} + v_n:
![alt text](https://github.com/YoujingYu99/IIA_SF1_Data_Analysis/blob/master/plots/an3.jpg?raw=true)

### 2. Denoising Audios
In an attempt to denoise noisy audios, two algorithms were proposed. 

#### Automatic Calculation of Noise Power Spectrum
The procedure is described as below:
    1. Plot waveform and determine a threshold for noise. 
    2. Check the wave amplitude at the time instances before, at and after the current time. If all of these have an amplitude below the threshold, this period is deemed as when no audio activity is present.
    3. Collect all periods with no audio (which is pure noise) and calculate the power spectrum.
    
The waveform for a sample audio before and after filtering is shown below:
![alt text](https://github.com/YoujingYu99/IIA_SF1_Data_Analysis/blob/master/plots/female_soft_time.jpg?raw=true)
    
#### Automatic Block-wise Calculation of Noise Power Spectrum
The procedure is described as below:
    1. Partition the audio file into blocks, each of block length N.
    2. Check the time instance before, at and after the current time. If all of these have an amplitude below the threshold, this period is deemed as when no audio activity is present.
    3. Collect all periods with no audio (which is pure noise) and calculate the power spectrum.
    4. Repeat for all blocks.
    
The waveform for a sample audio before and after filtering is shown below:
![alt text](https://github.com/YoujingYu99/IIA_SF1_Data_Analysis/blob/master/plots/male_vary_time.jpg?raw=true)


### 3.Finding the Signal

The task is to find the onset of a signal embedded in noisy signals. The posterior probabilities are calculated and plotted below:
![alt text](https://github.com/YoujingYu99/IIA_SF1_Data_Analysis/blob/master/plots/needle_posterior_probs.jpg?raw=true)


### 4. Recover Lost Signals
In this task, we aim to recover the lost signals. Two algorithms are proposed, namely the forward-backward algorithm and the Bayesian interpolation method.

#### Forward-Backward Algorithm
The forward-backward algorithm is described as below:

    1. Determine the sections in the lossy audio where there are no audio signals, their starting positions and ending positions.
    2. For each lost section, extract the audio section before and after the period.
    3. Apply AR model choice on the audio sections (by calculating the marginal likelihoods) to determine the most probable model order and generate ML and MAP estimates of the parameter θ.
    4. From the ML and MAP estimates, use the forward prediction formula, backward prediction formula and a weighting coefficient α to generate the lost signal.
    5. Map the generated signal onto the original lossy signal. 
    6. Iterate for all lost sections.
    7. Load original signal and calculate the mean-squared-error (MSE).


The original, lossy and regenerated waveforms are plotted below:
![alt text](https://github.com/YoujingYu99/IIA_SF1_Data_Analysis/blob/master/plots/AR_armst_recovery.jpg?raw=true)

#### Bayesian Interpolation
The Bayesian interpolation algorithm is described as below:
    1. Determine the sections in the lossy audio where there are no audio signals, their starting positions and ending positions.
    2. For each lost section, extract the audio section before and after the period.
    3. Apply AR model choice on the audio sections to determine the most probable model order and generate ML and MAP estimates of the parameter θ.
    4. Generate matrix A, where the ML or MAP estimated filter coefficients of the all-pole filter are mapped into the matrix.
    5. Partition matrix A by columns according to the lengths of the audio section before the lost period, the period with no signal present and the audio section after the period.
    6. Reregenerate lost data. Assign the regenerated data back into the lossy audio.
    7. Iterate for all lost sections.
    8. Load original signal and calculate the mean-squared-error (MSE).

The original, lossy and regenerated waveforms are plotted below:
![alt text](https://github.com/YoujingYu99/IIA_SF1_Data_Analysis/blob/master/plots/AR_armst_recovery_2.jpg?raw=true)


#### Sinusoidal Model
The sinusoidal model can be used to derive a Bayesian version of the discrete fourier transform (DFT) for noisy signals as well. The algorithm is described as below:

    1. A single-sided frequency spectrum of the audio is plotted.
    2. Peaks and their locations in the frequency spectrum are identified as the most significant components. The minimum peak prominence and minimum peak separation is set to 3.2 and 0.1 to filter out the noise.
    3. There are a total of 27 peaks identified, which are the most significant frequencies. 
    4. Perform Bayesian estimation on each single frequency component, where the signal amplitude (in the log domain) is obtained from the peak height, the prior mean is specified as the inverse of the signal amplitude, and the prior variance is chosen as 1. From this step we obtain the regularised relative coefficients of Bayesian amplitude of each frequency identified.
    5. Now we consider a multiple frequency model in which the sinusoids have frequencies set equal to the normalised most significant frequencies just identified, and the amplitude set as the Bayesian amplitude of each frequency identified in the previous step. This sets up the design matrix G.
    6. Reconstruct the signal from the multiple frequency model and listen to it to confirm that the MAP reconstructed signal sound almost identical to the original signal and are hence correctly regenerated.
    
The original, lossy and regenerated waveforms are plotted below:
![alt text](https://github.com/YoujingYu99/IIA_SF1_Data_Analysis/blob/master/plots/sinu_organ_prior.jpg?raw=true)

