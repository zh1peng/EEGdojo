addpath(genpath('/media/NAS/misc/matlab_toolbox/eeglab2023.1'))
params = struct;
params.ChanLocation = '/media/NAS/misc/matlab_toolbox/eeglab2023.1/plugins/dipfit/standard_BEM/elec/standard_1005.elc';
params.WriteSidecars = false;

neuracle2bids('/media/NAS/EEGdata/BEAT/sourcedata', params);