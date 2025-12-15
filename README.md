# TrimRunner.py
Code for running trim via a python script. Uses trimUtils.py to actually generate the TRIM input. Each trim instance runs in its own temporary folder. This project is designed to simulate primaries 5-10% above the highest energy you will need, and then in the parser script we save tracks from lower energy particles (primaries after collisions) to avoid needing to do thousands of individual TRIM sims. Output is the standard TRIM ASCII output, gzip'd.

# CollisionParser.cpp
C++ code for parsing TRIM COLLISON.txt.gz (sic) files. Saves output as a root file. Included is a unix Makefile as well as CMakelists and make.bat script for running on windows. While the TRIM sims are not "binned" in energy, we impose an artificial binning in the output of this code so limit the number of individual tracks stored at each primary energy. These files can get very large, so it is recommended to set the min/max energy you care about, as well as "bin" size and max number of tracks per "bin" to keep the output manageable. NOTE: this is just a library of TRIM tracks, it does NOT incorporate efficiencies--some recoils near threshold may produce 0 vacancies, so a separate efficiency curve must be calculated/used.

# trimEfficiencyCalculator.py
Code for calculating the vacancy production efficiency as a function of energy. The output is a CSV file. The min/max energy ranges can be supplied, but not that we stop calculating efficiencies when it is greater than 0.99999. To calcualte efficiencies we run individual TRIM sims at specific energies, parse the output files to find the number of collisions producing a vacancy, and divide that by the number of primaries to calculate the efficiency. The TRIM sims are deleted after this calculation is completed.

