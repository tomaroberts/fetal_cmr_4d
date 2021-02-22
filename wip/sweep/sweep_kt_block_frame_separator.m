function [thetaFrameMulti, idxFrameMulti] = sweep_kt_block_frame_separator( S, PARAM, iStk )

% Separates k-t SWEEP reconstucted data with overlapping blocks into
% datasets the same size as the original FOV
%
% - S     = cardsync structure
% - PARAM = output from ktrecon code
% - iStk  = index of stack
%
%

if nargin == 2
    iStk = 1;
end

% Load
swpWindows = PARAM.Sweep.swpWindows(:);
thetaFrame = cell2mat( S(iStk).thetaFrame );

% Count number of frame reconstructions
uv       = unique( swpWindows );
n_uv     = histc( swpWindows, uv );
n_uv_max = max(n_uv); % max number of frame duplicates

% initalise
frameNaNs  = nan( 1, max(swpWindows) );
for ii = 1:n_uv_max
    idxFrameMulti{1,ii}   = frameNaNs;
    thetaFrameMulti{1,ii} = frameNaNs;
end

% Find how many times each frame location was reconstructed and assign
% frames with multiple reconstructions to thetaFrameMulti
for target_val = 1:max(swpWindows)

    idx_duplicate_frames = find( swpWindows == target_val );
    num_duplicate_frames = numel( idx_duplicate_frames );   
    
    for n = 1:num_duplicate_frames
        idxFrameMulti{1,n}(target_val)   = idx_duplicate_frames(n);
        thetaFrameMulti{1,n}(target_val) = thetaFrame( idx_duplicate_frames(n) );
    end
    
end

% end fn
end

