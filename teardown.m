function teardown()
    % TEARDOWN - Remove the EEGdojo toolbox from the MATLAB path.
    rmpath(genpath(fileparts(mfilename('fullpath'))));
end