
function X = preprocess(X,eogchannels,badchannels,fs)
% All the usual EEG preprocessing, except epoching and epoch rejection as
% there are not events to epoch with natural stimuli. Instead, bad
% data is set = 0 in the continuous stream, which makes sense when
% computing covariance matrices but maybe not for other purposes. 
% Code adapted from eegisc.m version 0.05 dated 08/05/16.

debug = 0;   % turn this on to show data before/after preprocessing. 
stdThresh=4; % samples larger this many STDs will be marked as outliers
HPcutoff =0.5; % HP filter cut-off frequequency in Hz

% pick your preferred high-pass filter
[b,a,k]=butter(5,HPcutoff/fs*2,'high'); sos = zp2sos(b,a,k);
[bn,an]=butter(4,[58 61]./(fs/2),'stop'); 
          
[T,D,N]=size(X); 

% Preprocess data for all N subjects
for i=1:N
    
    data = X(:,:,i);

    % remove starting offset to avoid filter trancient
    data = data-repmat(data(1,:),T,1);
    
    % show the original data
    if debug, subplot(2,1,1); imagesc((1:T)/fs,1:D,data'); title(['Subject ' num2str(i)]); end

    % high-pass filter
    data = sosfilt(sos,data);          
    
    % regress out eye-movements;
    data = data - data(:,eogchannels) * (data(:,eogchannels)\data); 

    % detect outliers above stdThresh per channel; 
    data(abs(data)>stdThresh*repmat(std(data),[T 1])) = NaN;
    
    % remove 40ms before and after;
    h=[1; zeros(round(0.04*fs)-1,1)];    
    data = filter(h,1,flipud(filter(h,1,flipud(data))));
    
    % Mark outliers as 0, to avoid NaN coding and to discount noisy channels
    data(isnan(data))=0;

    % also zero out bad channels
    data(:,badchannels{i})=0; 
    
    % show the result of all this
    if debug, subplot(2,1,2); imagesc((1:T)/fs,1:D,data'); caxis([-100 100]); xlabel('Time (s)'); drawnow; end

    X(:,:,i) = data;
    
end

% discard the eog channels
X = X(:,setdiff(1:D,eogchannels),:);

