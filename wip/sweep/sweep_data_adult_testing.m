
% Testing number of frames in sweep acquisition required for accurate heart
% rate estimate from x-f data

% - Adult data
% - Code adapted from estimate_heartrate_xf.m

cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Adult\2019_03_28_Adult_Heart\swp_3d_python');

xtRcn = load_untouch_nii('IMG_3D_swp_s11_rlt_ab.nii.gz');
xtRcn = double(xtRcn.img);

% implay_RR(xtRcn);

lastFrame = [1020 1040 1060 1080 1100 1120 1140 1200]; % change number of frames used to determine HR

for ff = lastFrame
    
    
    subsetFrames = 1001:ff; % subset of frames to perform temporal FT
    imSeq = xtRcn(:,:,subsetFrames);
    
    
    %% Draw heart roi
    if ~exist('roi','var')
        imtar(imSeq(:,:, round(median(1:numel(subsetFrames))) ) );
        roi = roipoly;
        close all;
    end
    
    
    %% fn inputs:
    hrRange         = [40 110]; % typical adult
    useUpsampleFreq = true;
    frameDuration   = 0.0660; %taken from s11_rlt_parameters.mat
    visibleStr      = 'on';
    
    % TODO: convert following to optional input
    peakProminenceRatioMax = 2;
    peakWidthRatioMin      = 2;
    f0maxDiff              = 0.2; % Hz
    
    
    %% Below taken from estimate_heartrate_xf.m
    
    
    % Dimensions
    
    [nX,nY,nDyn,nChan] = size(imSeq);
    
    
    % Realtime Image Sequence
    
    if nChan > 1  % ensure single-channel images
        imSeq = sqrt( sum( imSeq.^2, 4 ) );
    end
    
    
    % Timing
    
    minHR = min( hrRange );  % bpm
    maxHR = max( hrRange );  % bpm
    
    minRR = 60/maxHR;  % s
    maxRR = 60/minHR;  % s
    
    minFreq = 1/maxRR; % Hz
    maxFreq = 1/minRR; % Hz
    
    if useUpsampleFreq
        dt = 0.0001;                   % target temporal resolution for estimation of heartrate
        df = 1/(maxRR) - 1/(maxRR+dt); % frequency resolution at min expected fundamental frequency
        % TODO: update to read optional input argument df instead of using two lines above
        nF = range(calc_freq( nDyn, frameDuration )) / df;  % number of x-f frames to reconstruct
    else
        nF = nDyn;
    end
    
    
    % Check against ROI
    
    if isempty( roi )
        
        warning('No ROI specified, using full FOV'),
        
        roi = true(nX,nY);
        
        nF = nDyn;  % don't upscale frequency resolution for full-FOV ROI since
        % reconstructing many more frames in x-f space for full-FOV
        % ROI is computationally demanding
        
    else
        
        nVox = sum(roi(:));
        
        while nF*nVox > 1e8
            nF = nF/10;
        end
        
    end
    
    
    % x-t to x-f
    
    xt2xf = @( xt ) fftshift( fft( xt, [], 3 ), 3 );
    
    
    %% Get Image Signal in ROI
    
    
    xtRoi = nan( sum(roi(:)), 1, nDyn );
    
    for iF = 1:nDyn
        imFrame = imSeq(:,:,iF);
        xtRoi(:,:,iF) = imFrame(roi);
    end
    
    
    %% Get x-f Magnitude in ROI
    
    
    xfRoi = abs( xt2xf( padarray( xtRoi, [0,0,ceil((nF-nDyn)/2)] ) ) );
    
    xfMeanSig = nan(1,size(xfRoi,3));
    
    for iF = 1:length(xfMeanSig)
        xfMeanSig(iF) = mean( xfRoi(:,:,iF) );
    end
    
    xfMeanSig = smooth( xfMeanSig, nF/nDyn, 'moving' );  % remove ringing introdcued by zero-padding x-t space in time to upscale frequency resolution
    
    
    %% Identify Peaks in Frequency Spectrum
    
    
    % frequencies
    
    f = calc_freq( size(xfRoi,3), frameDuration );
    
    
    % find frequencies in range of fundamental (f0)
    
    indf0pos = find(f<minFreq,1,'last'):find(f>maxFreq,1,'first');
    indf0neg = find(f<-maxFreq,1,'last'):find(f>-minFreq,1,'first');
    [ ~, iF0pos, ~, peakProminenceF0pos ] = findpeaks( xfMeanSig(indf0pos), 'NPeaks', 1, 'SortStr', 'descend' );
    [ ~, ~, peakWidthF0pos ]              = findpeaks( xfMeanSig(indf0pos), f(indf0pos), 'NPeaks', 1, 'SortStr', 'descend' );
    [ ~, iF0neg, ~, peakProminenceF0neg ] = findpeaks( xfMeanSig(indf0neg), 'NPeaks', 1, 'SortStr', 'descend' );
    [ ~, ~, peakWidthF0neg ]              = findpeaks( xfMeanSig(indf0neg), f(indf0neg), 'NPeaks', 1, 'SortStr', 'descend' );
    f0pos = abs( f(indf0pos(iF0pos)) );
    f0neg = abs( f(indf0neg(iF0neg)) );
    fPeak = [ indf0neg(iF0neg) indf0pos(iF0pos) ];
    
    
    %% Identify Fundamental Frequencies
    
    peakLocDiff = abs( f0neg - f0pos );
    
    peakProminence = mean([peakProminenceF0pos,peakProminenceF0neg]);
    if isempty( peakProminenceF0pos )
        peakProminenceF0pos = 0;
    end
    if isempty( peakProminenceF0neg )
        peakProminenceF0neg = 0;
    end
    
    peakWidth = mean([peakWidthF0pos,peakWidthF0neg]);
    if isempty( peakWidthF0pos )
        peakWidthF0pos = 0;
    end
    if isempty( peakWidthF0neg )
        peakWidthF0neg = 0;
    end
    
    % Resolve discrepency between peak locations
    if ( peakProminenceF0pos > peakProminenceRatioMax * peakProminenceF0neg ) %&& ( peakWidthRatioMin * peakWidthF0pos < peakWidthF0neg )
        fundamental = f0pos;
        peakProminence = peakProminenceF0pos;
        peakWidth      = peakWidthF0pos;
    elseif ( peakProminenceF0neg > peakProminenceRatioMax * peakProminenceF0pos ) %&& ( peakWidthRatioMin * peakWidthF0neg < peakWidthF0pos )
        fundamental = f0neg;
        peakProminence = peakProminenceF0neg;
        peakWidth      = peakWidthF0neg;
    else
        if abs( f0neg - f0pos ) <= f0maxDiff
            fundamental  = [ f0neg f0pos ];
        elseif ( peakProminenceF0pos > peakProminenceF0neg ) %&& ( peakWidthF0pos < peakWidthF0neg )
            fundamental = f0pos;
            peakProminence = peakProminenceF0pos;
            peakWidth      = peakWidthF0pos;
        elseif ( peakProminenceF0neg > peakProminenceF0pos ) %&& ( peakWidthF0neg < peakWidthF0pos )
            fundamental = f0neg;
            peakProminence = peakProminenceF0neg;
            peakWidth      = peakWidthF0neg;
        else
            fundamental  = [ f0neg f0pos ];
        end
    end
    
    if isempty( fundamental )
        error( '' )
    elseif length( fundamental ) == 1
        fundamental = [ fundamental, NaN ];
    end
    
    relativePeakProminence = peakProminence / max(xfMeanSig);
    
    
    %% Calculate Heartrate
    
    
    refFreq    = fundamental;
    
    rrInterval = mean( 1 ./ refFreq( ~isnan( refFreq ) ) );
    
    
    %% Calculate Cardiac Trigger Times
    
    
    nTrigger = ceil( nDyn * frameDuration / rrInterval );
    
    triggerTime = rrInterval * (0:nTrigger);
    
    
    %% Verbose
    
    
    % if ( isVerbose )
    
    figure( 'Name', 'heart_rate_estimate', 'Visible', visibleStr );
    
    plot( f, xfMeanSig ),  % plot mean x-f magnitude in ROI v. frequency
    hold on
    plot( f( fPeak), xfMeanSig( fPeak ), 'o')
    title('Heart Rate Estimation');
    xlabel('frequency (Hz)');
    ylabel('mean-signal_{ROI} (a.u.)');
    legend('signal','peaks'),
    hAx = gca;
    hAx.YTick = 0;
    a = axis;
    text(a(1)+0.05*(a(2)-a(1)),a(4)*0.9, sprintf('mean R-R interval = %.1f ms', rrInterval * 1000 ) ),
    if ( ~isempty( fundamental ) )
        text(a(1)+0.05*(a(2)-a(1)),a(4)*0.8, sprintf('f_{0} = %.3f, %.3f Hz => %.1f, %.1f ms', fundamental(1), fundamental(2), 1000/fundamental(1), 1000/fundamental(2) ) ),
    end
    
    axes('Position',[0.2 0.2 0.2 0.2])
    findpeaks( xfMeanSig(indf0neg), f(indf0neg), 'NPeaks', 1, 'SortStr', 'descend', 'Annotate', 'extents' );
    legend off
    
    axes('Position',[0.675 0.2 0.2 0.2])
    findpeaks( xfMeanSig(indf0pos), f(indf0pos), 'NPeaks', 1, 'SortStr', 'descend', 'Annotate', 'extents' );
    legend off
    
    if any(~roi(:))  % show ROI
        
        axes('Position',[.6 .425 .35 .35])
        im = abs( mean( imSeq, 3 ) );
        B = bwboundaries( roi, 'noholes' );  % NOTE: 'noholes' option added to avoid crashing in MATLAB2017a
        imshow( im, [ prctile(im(:),1), prctile(im(:),99)] )
        hold on
        for iB = 1:numel(B)
            line(B{iB}(:,2),B{iB}(:,1),'LineWidth',1,'Color','c')
        end
        
    end
    
    % end
    
    
    % end last frame loop
end