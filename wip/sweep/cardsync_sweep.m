function S = cardsync_sweep( S, varargin )
%CARDSYNC_SWEEP  Estimate heart rates for Sweep volume and sort data.
%
%   S = CARDSYNC_SWEEP( S ) 
%
%   CARDSYNC_SWEEP( ..., 'name', value ) specifies optional input argument. 
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
    
    % Get Required Parameters
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
    for iSwpLoca = 1:nSL(iStk)
    
        iDuplicateFrames             = find( P(iStk).Sweep.swpWindows(:) == iSwpLoca );
        S(iStk).tRRSwpLoca{iSwpLoca} = S(iStk).tRRSwpVec( iDuplicateFrames );
    
    end
   
    clear iDuplicateFrames
    
end


%% Plot Mean Heart Rate based on Multi-Sweep Window tRR Calculations

for iStk = 1:nStack

    % Calculate Max Number of Duplicate Frames (ugly method)
    for ii = 1:nSL(iStk)
        d(ii) = numel( cell2mat(S(iStk).tRRSwpLoca(ii)) );
    end
    maxDupl = max(d);

    duplicateHR = nan( maxDupl, length(S(iStk).tRRSwpLoca) );

    % This is very ugly...
    for iDupl = 1:maxDupl

        idx2get = find( cellfun('size',S(iStk).tRRSwpLoca,2) >= iDupl);
        for iSwpLoca = 1:numel(idx2get)
            duplicateHR(iDupl,idx2get(iSwpLoca)) =  60 ./ S(iStk).tRRSwpLoca{idx2get(iSwpLoca)}(iDupl);
        end

    end

    % Mean/Max/Min Heart Rate for Every Frame
    hrMeanSwpLoca = mean( duplicateHR, 1, 'omitnan' );
    hrMaxSwpLoca  =  max( duplicateHR );
    hrMinSwpLoca  =  min( duplicateHR );

    % Plot Mean, Min, Max HRs
    figure; hold on;
    plot(hrMeanSwpLoca, '-r'); % bodge for legend
    area(hrMaxSwpLoca,  'FaceColor',[0.9 0.9 0.9]);
    area(hrMinSwpLoca,  'FaceColor','white');
    plot(hrMeanSwpLoca, '-r');
    
    AX = axis;
    axis([1 length(S(iStk).tRRSwpLoca) round(AX(4)*0.8) AX(4)]);
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
% Sweep locations (SwpLoca) binned according to number of Sweep Windows and
% Sweep Window Stride

% Reminder:
% No. Sweep Bins = No. Sweep Locations / No. Sweep Windows + (2 * Sweep Win Stride)
% e.g.:
% nSL (no. loca) = 1088
% nSW (no. win)  = 32
% SWs (stride)   = 32
%                = 34 bins (due to stride being one-third of width)

for iStk = 1:nStack

    % Bin Frame Times (set frame 1 @ t = 0)
    tF = [];
    tF = S(iStk).tFrameSwpLoca - S(iStk).tFrameSwpLoca(1);
    tF = reshape( tF, [], nSW(iStk) );
    
    frameDuration = P(iStk).Timing.frameDuration;
    
    % Bin Cardiac Phases
    thetaF = [];

    % Mean RR for Each Bin
    tRRmean = [];
    tRRmean = 60 ./ hrMeanSwpLoca( 1:nSW(iStk):end );

    % Calculate Phase for Each Window
    for iRR = 1:nSW(iStk)

        nTrigger(iRR) = ceil( nSW(iStk) * frameDuration / tRRmean(iRR) );
        tRTriggermean(:,iRR) = tF(1,iRR) + tRRmean(iRR) * (0:nTrigger(iRR));

    %     % NB: this is slightly different to deltaPhase below
    %     % Here, the 33rd frame begins is offset according to RR of previous
    %     % window
    %     if iRR > 1
    %         deltaPhase = cPF(2,iRR-1) - cPF(1,iRR-1);
    %         cPF(1,iRR) = cPF(end,iRR-1) + deltaPhase;
    %     end

        [ ~, cardPhaseFraction(:,iRR) ] = calc_cardiac_timing( tF(:,iRR), tRTriggermean(:,iRR) );

    end

    % Offset cardPhase by Phase of Preceding Window and deltaPhase of Current Window
    deltaPhase  = 0;
    phaseOffset = 0;

    for iRR = 2:nSW(iStk) % First iRR does not need offset
        
        % Change in Cardiac Phase for Current Window
        deltaPhase(iRR)  = cardPhaseFraction(2,iRR) - cardPhaseFraction(1,iRR); 
        
        % Phase of Previous Frame
        phaseOffset(iRR) = cardPhaseFraction(end,iRR-1) + deltaPhase(iRR); % phase of previous frame
        
        % Apply Offset to Phases of Current Window
        cardPhaseFraction(:,iRR) = cardPhaseFraction(:,iRR) + phaseOffset(iRR); % offset phases of new window

    end
    
    % Convert to 2PI Range
    thetaFrameSwpBins = wrapTo2Pi( 2 * pi * cardPhaseFraction );
    
    % Assign to S
    for iRR = 1:nSW(iStk)
        S(iStk).thetaFrameSwpBins{iRR} = 2 * pi * thetaFrameSwpBins(:,iRR);
    end

    % Plot thetaFrame Graph
    thetaFrameSwpBins = thetaFrameSwpBins(:);

    figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    plot( thetaFrameSwpBins, 'ok', 'MarkerFaceColor','k' );
    plot( 1:nSW(iStk):nSL(iStk), thetaFrameSwpBins(1:nSW(iStk):nSL(iStk)), 'ok', 'Markersize', 15 );
    xticks( 1:nSW(iStk):nSL(iStk) );
    grid; grid minor;
    axis([1 nSL(iStk) -0.1 2*pi+0.1]);
    ylabel('Cardiac Phase (theta)');
    xlabel('Sweep Volume Location (z-position/frame index)');
    legend('Frame','First Frame of Bin','Location','NorthWest');
    
    % Save Fig
    hFig = gcf;
    hFig.Name = sprintf( 'sweep_cardiac_phase_binning_%s', S(iStk).desc );
    save_figs( pngfigDir, hFig, matfigDir )
    
    if ~isVerbose
        close( hFig )
    end

end





% cardsync_sweep(...)
end