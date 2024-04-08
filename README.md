# mm_phase_stim
Repository for vStr Phase Stim project

1. Run "getStimPhases.m" to generate "stim_phases.mat" inside each session folder.
"stim_phases.mat" contains the phase (correction of input inversion is applied here) in each of the 4 frequency bands at the time of stimulation.

2. Run "getCellStimResponse.m" to generate "<session_name>\_stim_response.mat" inside each session folder. This file contains pre and post stim window firing rates and spiking latencies for neurons, calculated for each stimulation epoch.

3. Run "getSpikingPhaseRelationshipsWithSameBinCount.m" to generate the main dependent variable (Firing-rate modulation by phase). This requires 1. and 2. to have run prior. Generates "<session_name>\_phase_response_5_bins.mat" and ".png". The number of bins can be changed in the script as a parameter.

4.  Run "getNonStimSpikePhases.m" to calculate the phases at the time of spikes that occurred outside of stimulation period, to generate "<session_name>\_nonstim_spk_phases.mat". The phases are binned into 25 bins and are later used to approximate the change in firing-rate expected from phase-change alone (in the absence of any stimulation). The spike-phases are also used to calculate spike-phase locking in neurons outside of stimulation periods. ECHT is used to generate phases (correction of input inversion is applied here).

5. Run "getNonStimPLV_causal.m" to quantify spike phaselocking using PLV using spikes that occurred outside of stimulation period, to generate "<session_name>\_spike_phaselock_causal_plv.mat". Requires 4. to be run prior.

6. Run "searchForLFPArtifacts.m" to generate image dump of how stim affect phase estimation by vanilla Hilbert transform.

7. Run "getPSDandBroadBandfit.m" to generate session-wise PSD and 1/f fit as calculated by IRASA as well as FOOOF


# Matlab path
```
addpath(genpath('mm_phase_stim\code-matlab\shared'));

addpath('mm_phase_stim\fieldtrip');

addpath(genpath('mm_phase_stim\fieldtrip\fileio\'));

addpath(genpath('mm_phase_stim\fieldtrip\utilities'));

addpath(genpath('mm_phase_stim\fieldtrip\contrib\spike'));

addpath(genpath('mm_phase_stim\fieldtrip\specest'));

addpath(genpath('mm_phase_stim\fieldtrip\preproc'));

addpath(genpath('mm_phase_stim\fieldtrip\external\brainstorm')); 

addpath(genpath('mm_phase_stim\code-matlab\graveyard'));

addpath(genpath('mm_phase_stim\code-matlab\ECHT'));

```
