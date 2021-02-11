%% FCMR cardsync simulation


%% Scan Simulation Setup
tStudy           = 3.211601000000000e+04;
tSeries          = 3.560101000000000e+04;

% seriesDuration   = 60.993000030517580;
numLoc           = 7;
locOrder         = 1:numLoc;
numDyn           = 96;
numDynDummy      = 0;
numFrame         = numLoc * numDyn; 
frameDuration    = 0.073;
seriesDuration   = ( numLoc * numDynDummy * frameDuration ) + ( numLoc * numDyn * frameDuration );
% nb: seriesDuration ^^ does not match exactly with scanner value for M2D

iStk             = 1;
S(iStk).desc     = 'sim-stack';

mean_tRR         = 0.45; % seconds
simulationMethod = 'HRvariation';


%% Sweep Window Setup

SWP(iStk).windowWidth   = 96;
SWP(iStk).windowStride  = 32;
SWP(iStk).windowCentre  = SWP(iStk).windowWidth/2 : SWP(iStk).windowStride : SWP(iStk).windowWidth*numLoc;

swpWinCentreMax    = SWP(iStk).windowCentre + SWP(iStk).windowWidth/2;   
swpWinWithinBounds = swpWinCentreMax <= (numLoc * numDyn);
numWithinBounds    = numel(find(swpWinWithinBounds));
SWP(iStk).windowCentre(numWithinBounds+1:end) = [];

% SWP(iStk).numSwpWin     = (numLoc * numDyn) / SWP(iStk).windowStride;
SWP(iStk).numSwpWin     = numel( SWP(iStk).windowCentre );

fn_windowIdx = @(winCentre, winWidth) winCentre-(winWidth/2)+1:winCentre+(winWidth/2);

for iWin = 1:SWP(iStk).numSwpWin
    SWP(iStk).frameIdxArray{1,iWin} = fn_windowIdx( SWP(iStk).windowCentre(iWin), SWP(iStk).windowWidth );
end

SWP(iStk).maxM2DFrameIdx = max( cell2mat(SWP(iStk).frameIdxArray) );
SWP(iStk).maxSwpFrameIdx = SWP(iStk).numSwpWin * SWP(iStk).windowStride;


%% Simulated Heart Beat

switch simulationMethod
    
    case 'HRconstant'
        S(iStk).tRR   = mean_tRR .* ones(1,numLoc);
        SWP(iStk).tRR = mean_tRR .* ones(1,SWP.numSwpWin);
    
        B = mean_tRR .* ones(1,numFrame);
        
    case 'HRvariation'

        nSamples             = numFrame;
        x                    = 1:nSamples;
        periodicityFactor    = 0.002;
        var_tRR              = 0.1; % 10% variation based on eyeballing existing HRs 
        smoothingKernelWidth = 50;
        
        A = mean_tRR + ( var_tRR*mean_tRR .* cos( 2*pi*periodicityFactor*x + 2*pi*rand ) + ( var_tRR*mean_tRR * randn(1,nSamples) ) );
        B = smoothdata( A, 'gaussian', smoothingKernelWidth );
        % figure; plot(x,A,x,B,'-bx');
        
        % Calculate Mean tRR values (Discretize based on numLoc)
        % M2D
        nSamplesPerSlices = floor( numel(x)/numLoc );
        true_tRR_re       = reshape( B(1:nSamplesPerSlices*numLoc), nSamplesPerSlices, [] );
        S(iStk).tRR       = mean( true_tRR_re, 1 );
        
        C = nan(size(x));
        C(1:nSamplesPerSlices*numLoc) = repelem( S(iStk).tRR, nSamplesPerSlices );
        
        % SWEEP
        nSamplesPerSlices = floor( numel(x)/SWP(iStk).numSwpWin );
        true_tRR_re       = reshape( B(1:nSamplesPerSlices*SWP(iStk).numSwpWin), nSamplesPerSlices, [] );
        SWP(iStk).tRR     = mean( true_tRR_re, 1 );
        
        D = nan(size(x));
        D(1:nSamplesPerSlices*SWP(iStk).numSwpWin) = repelem( SWP(iStk).tRR, nSamplesPerSlices );
        
        figure; plot( x,B,'-k.', x,C,'-b', x,D,'-r', 'LineWidth', 1.5 );
%         figure; plot( x,C,'-b', x,D,'-r', 'LineWidth', 2 );
        AX = axis; axis([AX(1) AX(2) 0 AX(4)*1.25]);
        xlabel('Frame Number'); ylabel('Heart Rate [seconds]');
        
        legend('Simulated HR', 'M2D HR', 'SWEEP HR');
        
end

% Ground Truth Heart Beat
SIM.true_tRR = B;

SIM.thetaFrame(1) = 2 * pi * rand;
for iFrame = 1:numFrame-1
    deltaThetaFrame          = frameDuration / SIM.true_tRR(iFrame);
    SIM.thetaFrame(iFrame+1) = SIM.thetaFrame(iFrame) + (2 * pi * deltaThetaFrame);
end

SIM.thetaFrame = wrapTo2Pi( SIM.thetaFrame );
SIM.desc       = 'True-HR-simulation';
SIM.nLoc       = numLoc;
SIM.numFrame   = numFrame;


%% Calculate Slice Timing
sliceDuration    = seriesDuration / numLoc;
sliceStartOffset = numDynDummy * frameDuration;

sliceTimeOffset = nan(1,numLoc);
for iLoc = locOrder(:).'  % NOTE: for loop index values should be row vector
    sliceTimeOffset(iLoc) = ( tSeries + double(iLoc-1) * sliceDuration + sliceStartOffset - tStudy );  % seconds
end

sliceTime        = sliceTimeOffset(:)';
sliceTimeOffset  = sliceTimeOffset(:)' - tSeries;

% Struct
P.Timing.tStudy           = tStudy;
P.Timing.tSeries          = tSeries;
P.Timing.seriesDuration   = seriesDuration;
P.Timing.numLoc           = numLoc;
P.Timing.numFrame           = numFrame;
P.Timing.locOrder         = locOrder;
P.Timing.numDyn           = numDyn;
P.Timing.numDynDummy      = numDynDummy;
P.Timing.frameDuration    = frameDuration;
P.Timing.sliceDuration    = sliceDuration;
P.Timing.sliceStartOffset = sliceStartOffset;
P.Timing.sliceTime        = sliceTime;
P.Timing.sliceTimeOffset  = sliceTimeOffset;


%% Calculate Frame Timing

S(iStk).nLoc = P.Timing.numLoc;

for iLoc = 1:S(iStk).nLoc

    % Dynamic Image Series
    S(iStk).frameDuration   = P.Timing.frameDuration;
    S(iStk).tFrame{iLoc}    = P.Timing.sliceTime(iLoc) + S(iStk).frameDuration * (0:(P.Timing.numDyn -1));

end

% Map tFrame to Sweep Windows
tFrameM2DVec = cell2mat(S(iStk).tFrame);
for i = 1:SWP(iStk).maxM2DFrameIdx

    current_tFrame = tFrameM2DVec(i);
    
    for iWin = 1:SWP(iStk).numSwpWin
        for iFrame = 1:SWP(iStk).windowWidth            
            if ( SWP.frameIdxArray{1,iWin}(iFrame) == i )
                SWP.tFrame{1,iWin}(iFrame) = current_tFrame;
            end   
        end
    end
    
end


%% Simulate Cardiac Timing

for iLoc = 1:S(iStk).nLoc    
    rrInterval  = S(iStk).tRR(iLoc);
    nTrigger    = ceil( P.Timing.numDyn * P.Timing.frameDuration / rrInterval );
    triggerTime = rrInterval * (0:nTrigger);
    S(iStk).tRTrigger{iLoc} = triggerTime + S(iStk).tFrame{iLoc}(1); 
end

for iWin = 1:SWP(iStk).numSwpWin    
    rrInterval  = SWP(iStk).tRR(iWin);
    nTrigger    = ceil( SWP(iStk).windowWidth * P.Timing.frameDuration / rrInterval );
    triggerTime = rrInterval * (0:nTrigger);
    SWP(iStk).tRTrigger{iWin} = triggerTime + SWP(iStk).tFrame{iWin}(1); 
end


%% Estimate Cardiac Phases

% Calculate cardiac phase 
for iLoc = 1:S(iStk).nLoc    
    if ~any( isnan( S(iStk).tRTrigger{iLoc} ) )
        [ ~, cardPhaseFraction ] = calc_cardiac_timing( S(iStk).tFrame{iLoc}, S(iStk).tRTrigger{iLoc} );
        S(iStk).thetaFrame{iLoc} = 2 * pi * cardPhaseFraction;
    else
        S(iStk).thetaFrame{iLoc} = nan( size( S(iStk).tFrame{iLoc} ) );
    end
end

for iWin = 1:SWP(iStk).numSwpWin     
    if ~any( isnan( SWP(iStk).tRTrigger{iWin} ) )
        [ ~, cardPhaseFraction ] = calc_cardiac_timing( SWP(iStk).tFrame{iWin}, SWP(iStk).tRTrigger{iWin} );
        SWP(iStk).thetaFrame{iWin} = 2 * pi * cardPhaseFraction;
    else
        SWP(iStk).thetaFrame{iWin} = nan( size( SWP(iStk).tFrame{iWin} ) );
    end
end


%% Multiple Sweep Estimates of Cardiac Phase
% Map Sweep Cardiac Phase Estimates to Single Frames
idxFrameSWPVec   = cell2mat( SWP(iStk).frameIdxArray );
thetaFrameSWPVec = cell2mat( SWP(iStk).thetaFrame );

maxNumEstimates = max( histc( idxFrameSWPVec, unique(idxFrameSWPVec) ) );
SWP(iStk).thetaFrameMultiEstimates = nan( maxNumEstimates, SWP.maxM2DFrameIdx );

for i = 1:SWP(iStk).maxM2DFrameIdx

    current_tFrames             = find( idxFrameSWPVec == i );
    current_thetaFrameEstimates = thetaFrameSWPVec( current_tFrames );
    
    SWP(iStk).thetaFrameMultiEstimates( 1:numel(current_thetaFrameEstimates), i ) = current_thetaFrameEstimates';
    
end

SWP(iStk).thetaFrameMean = nanmean( SWP(iStk).thetaFrameMultiEstimates );

% Plot
figure; hold on;
for iEst = 1:maxNumEstimates
    h(iEst) = plot( SWP(iStk).thetaFrameMultiEstimates(iEst,:), 'o' );
    set(h(iEst), 'MarkerFaceColor', get(h(iEst),'Color'));
end
% plot( SWP(iStk).thetaFrameMean, 'ok', 'MarkerFaceColor','k' );

xticks( 1:SWP(iStk).windowStride:SWP(iStk).maxM2DFrameIdx );
grid; grid minor;
axis([1 SWP(iStk).maxM2DFrameIdx+10 -0.1 2*pi+0.1]);

title( S(iStk).desc );
ylabel('Cardiac Phase (theta)');
xlabel('Frame Index');


%% Save Sm2d struct
Sm2d          = S;
Sm2d.desc     = 'M2D-simulation';
Sm2d.numFrame = numFrame;


%% Save Sswp struct
Sswp      = SWP;
Sswp.desc = 'SWEEP-simulation';
Sswp.nLoc = SWP(iStk).numSwpWin;


%% Sweep Method to Calculate thetaFrame

% Offset Cardiac Phases by Phase of Preceding Window and deltaPhase of Current Window
deltaPhase  = 0;
phaseOffset = 0;

for iRR = 2:S(iStk).nLoc % First iRR does not need offset

    % Change in Cardiac Phase for Current Window
    deltaPhase(iRR)  = S(iStk).thetaFrame{1,iRR}(2) - S(iStk).thetaFrame{1,iRR}(1); 

    % Phase of Previous Frame
    phaseOffset(iRR) = S(iStk).thetaFrame{1,iRR-1}(end) + deltaPhase(iRR); % phase of previous frame

    % Apply Offset to Phases of Current Window
    S(iStk).thetaFrame{1,iRR} = S(iStk).thetaFrame{1,iRR} + phaseOffset(iRR); % offset phases of new window

end

% Reassign to S
% S(iStk).thetaFrame = S(iStk).thetaFrameSwp;

% Wrap to 2PI range
for iRR = 1:S(iStk).nLoc
    S(iStk).thetaFrame{1,iRR} = wrapTo2Pi( S(iStk).thetaFrame{1,iRR} );
end


%% Save Sm2d_adjusted struct
Sm2d_adj      = S;
Sm2d_adj.desc = 'M2D-simulation-adjusted';


%% Plot Cardiac Phases

% plot_thetaframe( SIM, 'sim' );
% plot_thetaframe( Sm2d );
% plot_thetaframe( Sm2d_adj );
% plot_thetaframe( Sswp );


%% Comparison Plots

figure('units','normalized','outerposition',[0 0 1 1]);
xAX = Sswp(iStk).maxSwpFrameIdx;
ctr = 1;

for ii = [1,3]

    subplot(2,2,ctr);
    
    if ii == 1
        Scomp         = Sm2d;
        numLoc2Plot   = Sm2d(iStk).nLoc;
        numFrame2Plot = Sm2d(iStk).numFrame;
    elseif ii == 2
        Scomp         = Sswp;
        numLoc2Plot   = Sswp(iStk).numSwpWin;
        numFrame2Plot = Sswp(iStk).maxSwpFrameIdx;
    elseif ii == 3
        Scomp         = SWPopt;
        numLoc2Plot   = SWPopt(iStk).numSwpWin;
        numFrame2Plot = SWPopt(iStk).maxSwpFrameIdx;
    end
    
    Scomp.thetaFrameVec = cell2mat( Scomp.thetaFrame );

%     figure; 
    hold on;

    plot( SIM.thetaFrame,     'ok', 'MarkerFaceColor','k' );
    plot( Scomp.thetaFrameVec, 'or'  );

%     numFrame2Plot = numel( SIM.thetaFrame );

    % plot errors
    for iFrame = 1:numFrame2Plot
        thetaFrameDifference(iFrame) = Scomp.thetaFrameVec(iFrame)-SIM.thetaFrame(iFrame);

        if abs( thetaFrameDifference(iFrame) ) < pi
            plot([iFrame iFrame],[Scomp.thetaFrameVec(iFrame) SIM.thetaFrame(iFrame)],'r');
        elseif abs( thetaFrameDifference(iFrame) ) > pi
            if SIM.thetaFrame(iFrame) > Scomp.thetaFrameVec(iFrame)
                plot([iFrame iFrame],[SIM.thetaFrame(iFrame) 2*pi],'r');
                plot([iFrame iFrame],[0 Scomp.thetaFrameVec(iFrame)],'r');
            elseif Scomp.thetaFrameVec(iFrame) > SIM.thetaFrame(iFrame)
                plot([iFrame iFrame],[Scomp.thetaFrameVec(iFrame) 2*pi],'r');
                plot([iFrame iFrame],[0 SIM.thetaFrame(iFrame)],'r');
            end       
        end

    end

    xticks( 1:numFrame2Plot/numLoc2Plot:numFrame2Plot );
    grid; grid minor;
    axis([1 xAX+10 -0.1 2*pi+0.1]);
%     legend( SIM.desc, Scomp.desc );

    ctr = ctr+1;
    
    
    % Error Plot
%     figure;

    subplot(2,2,ctr);
 
    plot( abs( wrapToPi(thetaFrameDifference) ), 'o' );

    xticks( 1:numFrame2Plot/numLoc2Plot:numFrame2Plot );
    grid; grid minor;
    axis([1 xAX+10 -0.1 pi+0.1]);
    xlabel('Frame Number'); ylabel('Cardiac Phase Error: (Estimate - Simulation) [radians]');
    
    ctr = ctr+1;
    
end