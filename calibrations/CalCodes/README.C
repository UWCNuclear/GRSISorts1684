// 60Co  ../AnalysisTrees/analysis44483_sum.root
// 152Eu  ../AnalysisTrees/analysis44484_sum.root
// 56Co  ../AnalysisTrees/analysis44485_sum.root
// alpha  ../FragmentTrees/fragment44503_sum.root

First you need to add subrun together with gadd!

You first need to run ./SetupConfFile to get a good EmptyGains.cal with the right addresses at the righ channel numbers and set all offsets to zero and all gains to 1 (uncomment the commented lines).

There is a README file for the Ge calibrations, but to run it you just need to run the script

./energy.sh

Which calls the two calibration scripts linear_energy.C and quad_energy.C. When you run the script just follow the instructions in the terminal and it will do the calibration. Linear_energy.C doesn’t really do a linear calibration, you use a run with 60Co data (44483) and it does a quick calibration which quad_energy.C needs to work and do the actual calibration. I should comment that you may need to change the fit ranges, the binning is very different in the old TIGRESS DAQ versus the New TIGRESS DAQ so it may not quite work as intended. Just change the bin range though and it should be fine. Once the scripts have run you will get a file called

quad_energy_coeff.txt.

Copy this into the file CalibrationParameters.h and then make the script Conf_Setup.C (and CalFileConstructor.C?), the command to compile it is commented at the top of the script. Then run the compiled scripts;

./SetupConfFile
./ConstructCalibrationFile

This will make the CalibrationFile you need for grsisort.

I have also included a script for the silicon calibration

alphacal.cxx  

./ACal fragment.root EmptyGains.cal ACal.root

You need to compile it (the command is in the script). It will output the calibration parameters which you need to copy into CalibrationParameters.h, again remake Conf_Setup.C and run it and ConstructCalibrationFile and you will have your calibration file. 

