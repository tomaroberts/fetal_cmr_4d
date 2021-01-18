function plot_thetaframe( S, acqType )

%% PLOT_THETAFRAME
%
% Plot cardiac phase assigned to every frame
%
% Requires:
% - S structure normally in: /reconDir/cardsync/results_cardsync_*_.mat
%
%
% TAR (t.roberts@kcl.ac.uk)


%% Init

if nargin < 2
    acqType = false;
end
nStack   = numel( S );
 

%% Plot thetaFrame for all Loc

for iStk = 1:nStack

    thetaFrame = []; nFrame = []; nLoc = [];
    
    if strcmp(acqType, 'swp')
        thetaFrame = S(iStk).thetaFrameSwpLoc;
        nFrame     = numel( S(iStk).thetaFrameSwpLoc );
        nLoc       = S(iStk).nLoc;
    elseif strcmp(acqType, 'sim')
        thetaFrame = S(iStk).thetaFrame;
        nFrame     = numel( S(iStk).thetaFrame );
        nLoc       = S(iStk).nLoc;
    else
        thetaFrame = cell2mat( S(iStk).thetaFrame );
        nFrame     = numel( cell2mat( S(iStk).thetaFrame ) );
        nLoc       = S(iStk).nLoc;
    end
    
    figure('units','normalized','outerposition',[0 0 1 1]); 
    hold on;
    
    plot( thetaFrame, 'ok', 'MarkerFaceColor','k' );
    plot( 1:nFrame/nLoc:nFrame, thetaFrame( 1:nFrame/nLoc:nFrame ), 'ok', 'Markersize', 15 );
    
    xticks( 1:nFrame/nLoc:nFrame );
    grid; grid minor;
    axis([1 nFrame+10 -0.1 2*pi+0.1]);
    
    title( S(iStk).desc );
    ylabel('Cardiac Phase (theta)');
    xlabel('Frame Index');
    
    if ~acqType   
        legend('Frame','First Frame of Slice','Location','NorthWest');
    else
        legend('Frame','First Frame of Bin','Location','NorthWest');
    end
        
end


% plot_thetaframe(...)
end