# IR_Analysis_Synthesis
Apology
------
This repository is under construction - check back soon if code is missing or broken or contact James Traer (jtraer@mit.edu) if you have any urgent requests or if you would like to be informed when I finally have the code uploaded - 2016/Nov/08
Measuring Impulse Responses (IRs)
--------------------------------
To measure IRs broadcast a Golay sequence (created by running MakeGolaySequence.m - audio parameters can be specified therein) out of a loudspeaker and re-recording with a microphone. The broadcast and recorded signal should be synced to the same clock and hence a multitrack recorder is ideal (we used a Tascam DR40) but a laptop can also be used using standard audio software, as long as the same program controls both the broadcast and recording. The measurement technique is robust to noise and does not require silence nor loud broadcast volumes, however, the technique is sensitive to motion of thae apparatus (sturdy tripods are recommended) and distortion from clipping in the recoding or overdriving the speaker. 

The IR is extracted from the recorded track by ExtractIR.m (input path and analysis parameters specified therein). 

To account for distortions due to the equipment Calibration measurements can be made and included in the analysis - (if no calibration measurements are included the IR produced has included the electrical and mechanical filtering of the speaker and microphone).  To make calibration measurements use the same speaker, microphone and recorder as used in IR measurements and broadcast the same Golay sequence in as close to an anechoic environment as possible (e.g. a room damped with absorbent foam, curtains and carpeting, or a grassy field will do).  As most speakers have different transfer functions in different directions of broadcast, calibration measurements should be made with the speaker and microphone at a range of positions. We reccomend at minimum:
-- microphone facing speaker front
-- microphone facing speaker left
-- microphone facing speaker back
-- microphone facing speaker right
-- microphone facing speaker top
-- microphone facing speaker bottom
In all cases the microphone should face the speaker and the speaker-microphone separation distance should be constant. ExtractIR.m will accept an arbitrary number of calibration files (if calibration recordings are repeated the results will be averaged) at arbitrary positions. Metadata is required which ExtractIR.m reads from text files (Meta.txt -in the same directory as the audio).  If no text file is present ExtractIR.m will query the user to enter the information at the matlab command prompt and write the necessary file. Thereafter this data is read automatically.
