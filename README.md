# IR_Analysis_Synthesis

Apology
------
This repository is under construction - check back soon if code is missing or broken or contact James Traer (jtraer@mit.edu) if you have any urgent requests or if you have found major problems - Last updated 2017/June/19

## Overview

The matlab code in this repository is designed to measure Acoustic and Vibrational Impulse Responses (IRs) with low-intensity broadcasts, analyze the IRs and synthesize IRs that either replicate, or violate, the measured statistics of the measured IRs. To make a measurement a specific audio signal must be broadcast through the environment or object to be measured and re-recorded. The transfer-function of the apparatus can be measured with the same toolbox by placing the speaker and microphone directly adjacent and making a measurement.  This transfer function can then be removed from the other recordings.  The major work is done by the following scripts:
* MakeGolay.m - Synthesizes the required audio signal to be broadcast during measurement
* IR_Extract.m - Extracts the IR from the recorded bnroadcast (it must have access to the original signal)
* IR_Analysis.m - Analyzes the subsequent IRs to measure:
    * Spectra
    * Time for reverberation to decay 60dB (RT60) as a function of frequency
    * Direct-to-reverberant signal power ratio (DRR) as a function of frequency
    * Kurtosis (i.e. sparsity of echoes) as a function of time after the first arrival
    * Modes (frequency, onset power and RT60 of each identified). 
* IR_Statistics.m - Collates the results from IR_Analysis.m, groups them according to user specified Meta-tags, and plots the mean and variation across the groups

### Wrapper functions

The major files (IR_Extract.m,IR_Analysis.m and IR_Statistics.m) are all wrapper functions that navigate through a directory tree passing audio files and matlab data structures to a set of subroutines which do all the real processing. The subroutines create a directory with the same name as the recorded audio from a given measurement and they write their results therein as a set of plots (.jpg), audio files (.wav) and matlab structures (.mat).  The pipeline is as follows:
* IR_Extract.m feeds audio files (i.e the measurement) into hExtrct.m which saves the results into a structure (H.mat) 
    * IR_Extract.m is designed to be fast and simple, but it does require manual input to select the start and end times of the recorded IR and it prompts the user to manually label the data.  The labels are required for later sorting and plotting and the labelling process is done first to minimize any chance of the user forgetting when and how a measurement was recorded.
* IR_Analysis loads the aforementioned structures, feeds them into two analysis scripts. hPrp.m extracts frequency depedent RT60 and DRR and hExtrctMds.m searches for peaks in the spectrum to identify modes. The new properties are added to the pre-existing matlab structure and saved in a file named H_<NUMBER_OF_FREQUENCY_SUBBANDS>.mat.
    * IR_Analysis.m is slow and, with many IRs and many frequency channels it can take some time.  It requires no user input and can be run as a background process.
* IR_Statistics loads the aforementioned structures of analyzed IRs and extracts a list of Meta-tags.  Both the data and tags are passed to hPltStts.m which plots the results grouped according to the tags, to WrtDT2HTML.m which writes an html-file containing the IRs and information about them. 
    * IR_Statistics is relatively fast and may prompt the user to label data if, and only if, a label needed for pltting is missing. 

In many cases it might be easier to initially forego the wrappers and work directly with hExtrct.m, hPrp.m, hExtrctMds.m and hPltStts. The wrappers become helpful when the data set becomes large but initially they may create more confusion than they are worth.

### Meta-data

A meta-data file for each IR is written to a file Meta.txt within the directory of outputs for that file.  It is simply a text file of the format

VARIBALE_NAME_1   VALUE
VARIBALE_NAME_2   VALUE

When IR_Extract is run it looks for Meta.txt and reads what is in it.  If the user has queried information not in Meta.txt, the user is prompted to enter it manually and then the information is written to Meta.txt.  If teh information is already there, matlab reads the data, adds the information to the data-structure 'H.Meta' and moves on. Thus Meta data need only be entered once.  If IR_Extract.m is run again (with different parameters) and Meta.txt exists the user is not queried.

Similarly, the start and end time indices of the recorded IR are stored in a file labeled IR_Start_End.txt, and do not need to be re-entered if IR_Extract is run multiple times.

### Inputs

The Data required to run any of the wrappers is entered in a subroutine.  An example is given here as Input_Example.m. This specifies a structure (R) of paths to recorded audio files and the raw Golay sequences used in the broadcasts. A structure of Calibration recordings (C) - this can be empty for preliminary analyses. A structure of Meta-tags (Mt) - anything listed here will be queried by IR_Extract.m and will be used by IR_Statistics.m to group and plot data, and a structure of
an optional structure of amendments (Amnd) which further specifies how IRs are to be grouped in IR_Statistics.m.  

All the wrapper functions IR_Extract.m, IR_Analysis.m and IR_Statistics.m read the same input file.

### Known issues

* IR_Statistics.m has very recently gone through a big overhaul and - though functional - it is over-engineered, poorly commented and is easily broken by real-world data. 
* hExtrctMds.m (which is designed to find spectral peaks amd thereby select modes) is new and is both clunky, slow and not fully tested.  Modal analysis is still experimental. Similarly all the Plotting subroutines that plot modes are mostly ugly and untested.
* IR_Analysis.m produces a bloated structure with many fields saving overlapping data, much of which is probably useless.
* Currently some of the calibration steps are being done at the plotting stage.  Ideally this should all be done in IR_Analysis.m, such that the oututs are all calibrated and can be taken at face value later.  However, the calibration steps are easier to check and debug in IR_Statistics.m (which runs faster).  At some point, when I am convinced they are bug free and doing what they should, I plan to move them into IR_Analysis.


# Detailed Instructions

## To measure an IR
To measure IRs broadcast a Golay sequence (created by running MakeGolaySequence.m - audio parameters can be specified therein) out of a loudspeaker and re-recording with a microphone. The broadcast and recorded signal should be synced to the same clock and hence a multitrack recorder is ideal (we used a Tascam DR40) but a laptop can also be used using standard audio software, as long as the same program controls both the broadcast and recording. The measurement technique is robust to noise and does not require silence nor loud broadcast volumes, however, the technique is sensitive to motion of thae apparatus (sturdy tripods are recommended) and distortion from clipping in the recoding or overdriving the speaker. Though background noise is not fatal it does contribute to the noise floor and loud broadcasts are desirable if possible.  If not possible, the noise floor can be diminished by lengthening the broadcast time.

The IR is extracted from the recorded track by ExtractIR.m (input path and analysis parameters specified in the input file). 

## To analyze measured IRs

## To Calibrate the Apparatus

To account for distortions due to the equipment Calibration measurements can be made and included in the analysis - (if no calibration measurements are included the IR produced has included the electrical and mechanical filtering of the speaker and microphone).  To make calibration measurements use the same speaker, microphone and recorder as used in IR measurements and broadcast the same Golay sequence in as close to an anechoic environment as possible (e.g. a room damped with absorbent foam, curtains and carpeting, or a grassy field will do).  As most speakers have different transfer functions in different directions of broadcast, calibration measurements should be made with the speaker and microphone at a range of positions. We reccomend at minimum:
* microphone facing speaker front
* microphone facing speaker left
* microphone facing speaker back
* microphone facing speaker right
* microphone facing speaker top
* microphone facing speaker bottom
In all cases the microphone should face the speaker and the speaker-microphone separation distance should be constant. ExtractIR.m will accept an arbitrary number of calibration files (if calibration recordings are repeated the results will be averaged) at arbitrary positions. Metadata is required which ExtractIR.m reads from text files (Meta.txt -in the same directory as the audio).  If no text file is present ExtractIR.m will query the user to enter the information at the matlab command prompt and write the necessary file. Thereafter this data is read automatically.

## To plot the statistics of a set of analyzed IRs



