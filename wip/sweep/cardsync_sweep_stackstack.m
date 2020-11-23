
% rename?: cardsync_sweep_interstack.m
% and cardsync_sweep ---> cardsync_sweep_intrastack.m ?


%% Run recon_binned_sweep.m until apodization is complete

% Then ...


%% Stack-stack cardsync
% Maximise cross-correlation


%% Load Masks and Resample

% for iStk = 1:nStack
for iStk = 1
   
    N = load_untouch_nii( strcat( S(iStk).maskHeartFile(1:end-7), '_swp3d_apod.nii.gz' ) );
    M{iStk} = N.img;

    clear N
    
end


%% Load Cardiac Phases

% for iStk = 1:nStack
for iStk = 1
    
    thetaFrameSwpBins = [];
    thetaFrameSwpBins = cell2mat( S(iStk).thetaFrameSwpBins );
    thetaFrameSwpBins = thetaFrameSwpBins(:);
    
end