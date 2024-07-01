# mm_phase_stim
Repository for vStr Phase Stim project

1. Run "generateStimPhases.m" to generate "stim_phases.mat" inside each session folder.
"stim_phases.mat" contains the phase (correction of phase reversal is already applied) in each of the 4 frequency bands at the time of stimulation.

2. Run "generateCellStimResponse.m" to generate "<session_name>\_stim_response.mat" inside each session folder. This file contains pre and post stim window firing rates and spiking latencies for neurons, calculated for each stimulation epoch.

3. Run "generateSpikingPhaseRelationshipsWithSameBinCount.m" to generate the main dependent variable (Firing-rate modulation by phase). This requires 1. and 2. to have run prior. Generates "<session_name>\_phase_response_5_bins.mat" and ".png". The number of bins can be changed in the script as a parameter.

4. 
