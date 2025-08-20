function prep(params)
%PREP Run the EEG preprocessing pipeline.
%   Usage: prep(params);
%   params is a struct with fields corresponding to the properties of
%   EEGdojo.Params

    [~, baseName, ~] = fileparts(params.filename);
    outputfile = [baseName '_preprocessed.set'];

    EEG = pop_loadset('filename', params.filename, 'filepath', params.filepath);
    p = EEGdojo.Params(params);

    pipe = EEGdojo.Pipeline(EEG, p, fullfile(p.outputpath, 'preproc.log'));

    pipe = pipe.addStep(@EEGdojo.preprocess.remove_channels, 'labels', p.chan2remove);
    pipe = pipe.addStep(@EEGdojo.preprocess.downsample, 'freq', p.DownsamplingRate);
    pipe = pipe.addStep(@EEGdojo.preprocess.crop_by_markers, p.StartMarker, p.EndMarker, 'PadSec', p.PadTime);
    pipe = pipe.addStep(@EEGdojo.preprocess.filter, 'HPfreq', p.HPfreq, 'LPfreq', p.LPfreq);
    pipe = pipe.addStep(@EEGdojo.preprocess.cleanline, 'freq', p.PowerLineFrequency);
    pipe = pipe.addStep(@EEGdojo.preprocess.remove_bad_channels, 'EOGLabels', p.eogchan, ...
        'FlatlineCriterion', p.FlatLineCriterion, 'ChannelCriterion', p.ChannelCriterion, 'LineNoiseCriterion', p.LineNoiseCriterion);
    pipe = pipe.addStep(@EEGdojo.preprocess.reref, 'excludeLabels', p.eogchan);
    pipe = pipe.addStep(@EEGdojo.preprocess.run_ica, 'EOGLabels', p.eogchan, 'ICLabelRejectMask', p.ICLabel);
    pipe = pipe.addStep(@EEGdojo.preprocess.interpolate);

    pipe = pipe.run();

    pop_saveset(pipe.EEG, 'filename', outputfile, 'filepath', p.outputpath);

end
