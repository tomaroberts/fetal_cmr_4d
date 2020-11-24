function S = cardsync_sweep_intrastack( S, varargin )
%CARDSYNC_SWEEP_INTRASTACK  Estimate heart rates for Sweep volume and sort data.
%
%   S = CARDSYNC_SWEEP_INTRASTACK( S ) 
%
%   CARDSYNC_SWEEP_INTRASTACK( ..., 'name', value ) specifies optional input argument. 
%   See code for name-value pairs.
%
%   See also PREPROC.

%   tar (t.roberts@kcl.ac.uk)


%% Optional Input Argument Default Values

default.resultsDir      = pwd;       % path of directory to save results
default.isVerbose       = true;


%% Parse Input

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end


addRequired(   p, 'S', ... 
    @(x) validateattributes( x, {'struct'}, {'vector'}, mfilename) );

add_param_fn(   p, 'resultsdir', default.resultsDir, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, S, varargin{:} );

resultsDir      = p.Results.resultsdir;
isVerbose       = p.Results.verbose;


%% Setup

nStack      = numel( S );
pngfigDir   = fullfile( resultsDir, 'figs', 'png' );
matfigDir   = fullfile( resultsDir, 'figs', 'fig' );


%% Load Parameters

for iStk = 1:nStack
    
    M       = matfile( S(iStk).rltParamFile );
    P(iStk) = M.PARAM;
    
    nSW(iStk) = P(iStk).Sweep.numSwpWindows;
    SWw(iStk) = P(iStk).Sweep.swpWinFullWidth;
    SWs(iStk) = P(iStk).Sweep.swpWinStride;
%     SWo = zeros(1,nSW);
    nSL(iStk) = max( P(iStk).Sweep.swpWindows(:) );
    
end


%% tRR of All Sweep Windows

xRange = 1:sum( nSW );
yRange = nan( nStack, sum(nSW) );
yCtr   = 1;

figure; hold on;
for iStk = 1:nStack

    yRange( iStk, yCtr:yCtr+S(iStk).nLoc-1 ) = cell2mat([S(iStk).tRR]);
    yCtr = yCtr + S(iStk).nLoc;

    plot( xRange, yRange, 'o' );
end
legend(S.desc);
title('tRR');
xlabel('Sweep Window Index');
ylabel('tRR (seconds)');
    
% Save Fig
hFig = gcf;
hFig.Name = sprintf( 'sweep_windows_rr_intervals' );
save_figs( pngfigDir, hFig, matfigDir )

if ~isVerbose
    close( hFig )
end


%% Mean Heart Rate based on multi-tRR Estimates from Overlapping Sweep Windows

% Collate Duplicate tRR Estimates at each Sweep Location

for iStk = 1:nStack
    
    % Assign tRR to Every Frame
    S(iStk).tRRSwpVec = repelem( cell2mat([S(iStk).tRR]), SWw(iStk) );
    
    % Collate Duplicate tRR Values at Every Location in Sweep Volume
    for iSwpLoc = 1:nSL(iStk)
    
        iDuplicateFrames           = find( P(iStk).Sweep.swpWindows(:) == iSwpLoc );
        S(iStk).tRRSwpLoc{iSwpLoc} = S(iStk).tRRSwpVec( iDuplicateFrames );
    
    end
   
    clear iDuplicateFrames
    
end


%% Plot Mean Heart Rate based on Multi-Sweep Window tRR Calculations

for iStk = 1:nStack

    % Calculate Max Number of Duplicate Frames (ugly method)
    for ii = 1:nSL(iStk)
        d(ii) = numel( cell2mat(S(iStk).tRRSwpLoc(ii)) );
    end
    maxDupl = max(d);

    duplicateHR = nan( maxDupl, length(S(iStk).tRRSwpLoc) );

    % This is very ugly...
    for iDupl = 1:maxDupl

        idx2get = find( cellfun('size',S(iStk).tRRSwpLoc,2) >= iDupl);
        for iSwpLoc = 1:numel(idx2get)
            duplicateHR(iDupl,idx2get(iSwpLoc)) =  60 ./ S(iStk).tRRSwpLoc{idx2get(iSwpLoc)}(iDupl);
        end

    end

    % Mean/Max/Min Heart Rate for Every Frame
    S(iStk).hrMeanSwpLoc = mean( duplicateHR, 1, 'omitnan' );
    S(iStk).hrMaxSwpLoc  =  max( duplicateHR );
    S(iStk).hrMinSwpLoc  =  min( duplicateHR );

    figure; hold on;
    plot(S(iStk).hrMeanSwpLoc, '-r'); % bodge for legend
    area(S(iStk).hrMaxSwpLoc,  'FaceColor',[0.9 0.9 0.9]);
    area(S(iStk).hrMinSwpLoc,  'FaceColor','white');
    plot(S(iStk).hrMeanSwpLoc, '-r');
    
    AX = axis;
    axis([1 length(S(iStk).tRRSwpLoc) round(AX(4)*0.8) AX(4)]);
    title([S(iStk).desc ' - Sweep Windows Heart Rate Variation']);
    xlabel('Sweep z-Location/Frame Number');
    ylabel('Heart Rate [bpm]');
    legend('Mean','Min/Max', 'Location', 'SouthWest');
    
    % Save Fig
    hFig = gcf;
    hFig.Name = sprintf( 'sweep_windows_hr_variation_%s', S(iStk).desc );
    save_figs( pngfigDir, hFig, matfigDir )
    
    if ~isVerbose
        close( hFig )
    end

end


%% Calculate Phase of Every Frame in Sweep Volume
% Sweep locations (SwpLoc) binned according to number of Sweep Windows and
% Sweep Window Stride

% Reminder:
% No. Sweep Bins = No. Sweep Locations / No. Sweep Windows + (2 * Sweep Win Stride)
% e.g.:
% nSL (no. loca) = 1088
% nSW (no. win)  = 32
% SWs (stride)   = 32
%                = 34 bins (due to stride being one-third of width)

for iStk = 1:nStack

    nBins(iStk) = nSW(iStk) + 2; 
    % TODO: 
    % - not sure this is will work for all strides
    % - Might be problematic if Sweep Window Width not divisible by Stride    
    
    % Bin Frame Times (set frame 1 @ t = 0)
    tF = [];
    tF = S(iStk).tFrameSwpLoc - S(iStk).tFrameSwpLoc(1);
    tF = reshape( tF, [], nBins(iStk) ); % bin width = stride length
    
    frameDuration = P(iStk).Timing.frameDuration;
    
    % Bin Cardiac Phases
    thetaF = [];

    % Mean RR for each Bin
    tRRMean = [];
    tRRMean = 60 ./ S(iStk).hrMeanSwpLoc( 1:SWs(iStk):end ); % simply get first value for each bin

    
    % Calculate Cardiac Phases for Each Window
    for iRR = 1:nBins(iStk)
        
        S(iStk).nTrigger{iRR} = ceil( SWs(iStk) * frameDuration / tRRMean(iRR) );
        S(iStk).tRTriggerMean{iRR} = tF(1,iRR) + tRRMean(iRR) * ( 0:S(iStk).nTrigger{iRR} );
        
        [ ~, S(iStk).cardPhaseFraction(:,iRR) ] = calc_cardiac_timing( tF(:,iRR), S(iStk).tRTriggerMean{iRR} );
    
    end
    
    
    % Offset Cardiac Phases by Phase of Preceding Window and deltaPhase of Current Window
    deltaPhase  = 0;
    phaseOffset = 0;

    for iRR = 2:nBins(iStk) % First iRR does not need offset
        
        % Change in Cardiac Phase for Current Window
        deltaPhase(iRR)  = S(iStk).cardPhaseFraction(2,iRR) - S(iStk).cardPhaseFraction(1,iRR); 
        
        % Phase of Previous Frame
        phaseOffset(iRR) = S(iStk).cardPhaseFraction(end,iRR-1) + deltaPhase(iRR); % phase of previous frame
        
        % Apply Offset to Phases of Current Window
        S(iStk).cardPhaseFraction(:,iRR) = S(iStk).cardPhaseFraction(:,iRR) + phaseOffset(iRR); % offset phases of new window

    end
    
    
    % Convert to 2PI Range
    thetaFrameSwpBins = wrapTo2Pi( 2 * pi * S(iStk).cardPhaseFraction );
    
    % Assign to S
    for iRR = 1:nBins(iStk)
        S(iStk).thetaFrameSwpBins{iRR} = thetaFrameSwpBins(:,iRR);
    end

    % Plot thetaFrame Graph
    thetaFrameSwpBins = thetaFrameSwpBins(:);

    figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    plot( thetaFrameSwpBins, 'ok', 'MarkerFaceColor','k' );
    plot( 1:SWs(iStk):nSL(iStk), thetaFrameSwpBins(1:SWs(iStk):nSL(iStk)), 'ok', 'Markersize', 15 );
    xticks( 1:SWs(iStk):nSL(iStk)+SWs(iStk) );
    grid; grid minor;
    axis([1 nSL(iStk)+SWs(iStk) -0.1 2*pi+0.1]);
    ylabel('Cardiac Phase (theta)');
    xlabel('Sweep Volume Location (z-position/frame index)');
    legend('Frame','First Frame of Bin','Location','NorthWest');
    
    % Save Fig
    hFig = gcf;
    hFig.Name = sprintf( 'sweep_intrastack_cardiac_phases_%s', S(iStk).desc );
    save_figs( pngfigDir, hFig, matfigDir )
    
    if ~isVerbose
        close( hFig )
    end

end


%% Save Results to Text Files

% TODO: Do I need something like this here?

% fid = fopen( fullfile( resultsDir, 'mean_rrinterval.txt' ), 'w' );
% fprintf( fid, '%.6f ', mean( cell2mat( [ S.tRR ] ) ) );
% fclose( fid );
% fid = fopen( fullfile( resultsDir, 'rrintervals.txt' ), 'w' );
% fprintf( fid, '%.6f ', cell2mat( [ S.tRR ] ) );
% fclose( fid );
% fid = fopen( fullfile( resultsDir, 'cardphases_intraslice_cardsync.txt' ), 'w' );
% fprintf( fid, '%.6f ', cell2mat( [ S.thetaFrame ] ) );
% fclose( fid );


%% Save Results to .mat File

save( fullfile( resultsDir, 'results_cardsync_sweep_intrastack.mat' ), 'S', '-v7.3' );


% cardsync_sweep_intrastack(...)
end