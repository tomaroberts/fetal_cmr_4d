%% Plot cardsync variables
%
% - Created for examining timing of Sweep data
%
%
% tar (t.roberts@kcl.ac.uk)


%% Load
seriesNos = [20 22 23 24];
nStack    = numel(seriesNos);

reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Sweep\c_fcmr333';
cd(reconDir);

% cardsync_intra
load( fullfile( reconDir, 'cardsync', 'results_cardsync_intraslice.mat' ) );

totalNumSweepWin = sum([S.nLoc]);

for iStk = 1:nStack
    
    M       = matfile( S(iStk).rltParamFile );
    P(iStk) = M.PARAM;

end


%% tFrame - time of every frame

numDyn = P(1).Encoding.NrDyn(1);

% figure; hold on;
% 
% % Linear Sweep tFrames
% plot( S(iStk).tFrameSwpLoca, 'k-' );
% 
% % Sweep Windows
% plot( 1:numDyn, S(iStk).tFrame{1,1} );
% for iLoc = 2:S(iStk).nLoc
%     xRange = ( (iLoc-1) * numDyn ) + 1:...
%         ( (iLoc) * numDyn );
%     plot( xRange, S(iStk).tFrame{1,iLoc});
% end
% title(['Stack ID: ' S(iStk).desc]);
% xlabel('Frame Index');
% ylabel('Time');


%% tRR - estimated RR interval of every reconstruction (each Sweep Window)

xRange = 1:totalNumSweepWin;
yRange = nan(nStack,totalNumSweepWin);
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
ylabel('Mean HR (seconds)');


%% tRTrigger - R-wave/trigger times within a Sweep Window
% These are given by linear increment of tRR

% tRRvec = cell2mat([S(iStk).tRR]);
% 
% iWinA   = 10; 
% iWinB   = 11;
% 
% % RR
% tRRWinA = tRRvec(iWinA);
% tRRWinB = tRRvec(iWinB);
% tRTriggerA = cell2mat([S(1).tRTrigger(1,iWinA)]);
% tRTriggerB = cell2mat([S(1).tRTrigger(1,iWinB)]);
% 
% % Frame times
% tFrameWinA = cell2mat([S(1).tFrame(iWinA)]);
% tFrameWinB = cell2mat([S(1).tFrame(iWinB)]);
% 
% % Plot
% figure; hold on;
% plot([tFrameWinA,tFrameWinB],'o');
% 
% %%% BELOW --- trying to plot tRR alongside tFrame for visualisation
% % % AX     = axis;
% % % yMin   = AX(3); 
% % % xAxis  = yMin * ones(size([tFrameWinA,tFrameWinB]));
% % % yAxisA = linspace( min(tRTriggerA), max(tRTriggerA), length(tFrameWinA) );
% % % yAxisB = linspace( min(tRTriggerB), max(tRTriggerB), length(tFrameWinB) );
% % % yAxis  = [yAxisA yAxisB];
% % % 
% % % plot( zeros(192,1), yAxis, 'or' );


%% SWP struct

swpWinFullWidth = P(1).Sweep.swpWinFullWidth;

SWP = struct();

for iStk = 1:nStack
    
    SWP(iStk).tFrame     = cell2mat([S(iStk).tFrame]);
    SWP(iStk).thetaFrame = cell2mat([S(iStk).thetaFrame]);
    SWP(iStk).tRR        = repelem( cell2mat([S(iStk).tRR]), swpWinFullWidth );
%     SWP(iStk).tRTrigger  = repelem( cell2mat([S(iStk).tRTrigger]), swpWinFullWidth );
    SWP(iStk).tRTrigger  = cell2mat([S(iStk).tRTrigger]);
    
    % Find indices of duplicate frames
    numSwpLoca   = max( P(iStk).Sweep.swpWindows(:) );
    for iSwpLoca = 1:numSwpLoca
        
        iDuplicateFrames = find( P(iStk).Sweep.swpWindows(:) == iSwpLoca );
        
        SWP(iStk).swpLocaIdxArray{iSwpLoca}   = iDuplicateFrames;
        SWP(iStk).nFrameInstances(iSwpLoca)   = numel( iDuplicateFrames );
        SWP(iStk).tFrameSwpLoca{iSwpLoca}     = SWP(iStk).tFrame( iDuplicateFrames );
        SWP(iStk).thetaFrameSwpLoca{iSwpLoca} = SWP(iStk).thetaFrame( iDuplicateFrames );
        SWP(iStk).tRRSwpLoca{iSwpLoca}        = SWP(iStk).tRR( iDuplicateFrames );
    
    end
    
end


%% Examine Parameters of Duplicate Frames
iSwpLoca = 32+32+1;
fprintf('\nSweep Volume Slice Index = %i\n', iSwpLoca );
fprintf('   Sweep Duplicate Indices = %i\n', cell2mat( SWP(iStk).swpLocaIdxArray(iSwpLoca) ) );
fprintf('   Sweep Frame Times = %.3f\n',     cell2mat( SWP(iStk).tFrameSwpLoca(iSwpLoca) ) );
fprintf('   Sweep Phase (theta) = %.3f\n',   cell2mat( SWP(iStk).thetaFrameSwpLoca(iSwpLoca) ) );
fprintf('   Sweep RR interval = %.3f\n',     cell2mat( SWP(iStk).tRRSwpLoca(iSwpLoca) ) );
fprintf('   Sweep Heart Rate = %.0f\n',      60 ./ cell2mat( SWP(iStk).tRRSwpLoca(iSwpLoca) ) );


%% Plot Parameters of Duplicate Frame

% Method below:
% - Frame offset = 


iStk = 4;
numSwpLoca   = max( P(iStk).Sweep.swpWindows(:) );

maxDupl = 3;

x = nan(maxDupl,length(SWP(iStk).tRRSwpLoca));

% This is very ugly...
for iDupl = 1:maxDupl
    
    idx2get = find( cellfun('size',SWP(iStk).tRRSwpLoca,2) >= iDupl);
    
    for iSwpLoca = 1:numel(idx2get)
   
%         x(iDupl,idx2get(iSwpLoca)) = SWP(iStk).tRRSwpLoca{idx2get(iSwpLoca)}(iDupl);
        x(iDupl,idx2get(iSwpLoca)) =  60 ./ SWP(iStk).tRRSwpLoca{idx2get(iSwpLoca)}(iDupl);
        
    end
        
end


xmean = mean(x,1,'omitnan');
xmax  = max(x);
xmin  = min(x);

figure; hold on;
% plot(xmax,'-r');
% plot(xmin,'-r');
plot(xmean,'-r'); % bodge for legend
area(xmax,'FaceColor',[0.9 0.9 0.9]);
area(xmin,'FaceColor','white');
plot(xmean,'-r');
AX = axis;
axis([1 length(SWP(iStk).tRRSwpLoca) round(AX(4)*0.8) AX(4)]);
title('Sweep Windows Heart Rate Variation');
xlabel('Sweep z-Location/Frame Number');
ylabel('Heart Rate [bpm]');
legend('Mean','Min/Max');


% Convert xmean to cardiac phase for comparison with method further down
% TODO:
% - work out why why the plot is smooth. Think to do with calculation step between sweep windows
% - generalise this code if it turns out to be the best way to do it
tF = [];
thetamean = [];

WHY=S(iStk).nLoc+2;

tF = S(iStk).tFrameSwpLoca - S(iStk).tFrameSwpLoca(1);
tF = reshape(tF, [], WHY);

tRRmean = 60 ./ xmean(1:32:end);

nDyn = 32;
frameDuration = 0.073;

for iRR = 1:WHY
    
    nTrigger(iRR) = ceil( nDyn * frameDuration / tRRmean(iRR) );
    tRTriggermean(:,iRR) = tF(1,iRR) + tRRmean(iRR) * (0:nTrigger(iRR));
    
%     % NB: this is slightly different to deltaPhase below
%     % Here, the 33rd frame begins is offset according to RR of previous
%     % window
%     if iRR > 1
%         deltaPhase = cPF(2,iRR-1) - cPF(1,iRR-1);
%         cPF(1,iRR) = cPF(end,iRR-1) + deltaPhase;
%     end
    
    [ ~, cPF(:,iRR) ] = calc_cardiac_timing( tF(:,iRR), tRTriggermean(:,iRR) );
        
end

% offset by phase of preceding window and deltaPhase of new window
deltaPhase = 0;
phaseOffset = 0;

for iRR = 2:WHY
    deltaPhase(iRR) = cPF(2,iRR) - cPF(1,iRR); % change in cardiac phase for current window
    phaseOffset(iRR) = cPF(end,iRR-1) + deltaPhase(iRR); % phase of previous frame
    cPF(:,iRR) = cPF(:,iRR) + phaseOffset(iRR); % offset phases of new window
end
    
thetamean = wrapTo2Pi( 2 * pi * cPF );
thetamean = thetamean(:);

figure; hold on;
plot(thetamean,'ok','MarkerFaceColor','k');
plot(1:32:numSwpLoca,thetamean(1:32:numSwpLoca),'ko','Markersize',15);
xticks(1:32:numSwpLoca);
grid; grid minor;
axis([1 numSwpLoca -1 2*pi*1.1]);
ylabel('Cardiac Phase (theta)');
xlabel('Sweep Volume Location (z-position/frame index)');


%% Plot theta from Consecutive Sweep Windows

% iSwpWin = 1;
% swpWinStride    = P(iStk).Sweep.swpWinStride;
% swpWinFullWidth = P(iStk).Sweep.swpWinFullWidth;
% 
% 
% % figure; hold on;
% % plot(1:96,S(iStk).thetaFrame{1,iSwpWin}(:),'rx');
% % plot(33:96+32,S(iStk).thetaFrame{1,iSwpWin+1}(:),'bo');
% 
% thetaA                     = nan(1,swpWinFullWidth+swpWinStride);
% thetaA(1:swpWinFullWidth)  = S(iStk).thetaFrame{1,iSwpWin}(:);
% thetaB                     = nan(1,swpWinFullWidth+swpWinStride);
% thetaB(swpWinStride+1:end) = S(iStk).thetaFrame{1,iSwpWin+1}(:);
% diffTheta                  = thetaA - thetaB;
% 
% figure; hold on;
% plot(1:swpWinFullWidth,S(iStk).thetaFrame{1,iSwpWin}(:),'rx');
% plot(swpWinStride+1:swpWinFullWidth+swpWinStride,S(iStk).thetaFrame{1,iSwpWin+1}(:),'bo');
% plot(1:(swpWinFullWidth+swpWinStride),diffTheta,'k.','Markersize',10);
% axis([1 (swpWinFullWidth+swpWinStride) min(diffTheta(:))*1.2 2*pi*1.1]);
% title('Difference in Phase Calculated using Sequential Sweep Windows');
% xlabel('Sweep Frame Location');
% ylabel('Theta Difference [radians]');
% grid on; grid minor;











%% Theta Offset - ensure duplicate frames have matching cardiac phases

% for iStk = 1:nStack
%     
%     % Get theta which corresponds to start of next Sweep Window
%     for iSwpWin = 1:P(iStk).Sweep.numSwpWindows
%          
%         SWP(iStk).thetaFrameOffset(iSwpWin) = S(iStk).thetaFrame{1,iSwpWin}(swpWinStride+1);
%     
%     end
%     
% end


%% 
iSwpWin = 5;

win1 = S(iStk).thetaFrame{1,iSwpWin}; 
winOffset1 = win1(swpWinStride+1);

win2 = wrapTo2Pi( S(iStk).thetaFrame{1,iSwpWin+1} + winOffset1 ); 
% winOffset2 = win2(swpWinStride+1); % simply start previous window
winOffset2 = (win1(2*swpWinStride+1) + win2(swpWinStride+1))/2; % average of previous 2 windows

win3 = wrapTo2Pi( S(iStk).thetaFrame{1,iSwpWin+2} + winOffset2 );
% winOffset3 = win3(swpWinStride+1); % simply start previous window
winOffset3 = (win2(2*swpWinStride+1) + win3(swpWinStride+1))/2; % average of previous 2 windows

win4 = wrapTo2Pi( S(iStk).thetaFrame{1,iSwpWin+3} + winOffset3 );

figure; hold on;
plot(1:96,win1,'ro');
plot(33:96+32,win2,'bo'); 
plot(65:96+64,win3,'go'); 
plot(97:96+96,win4,'mo'); 
grid; grid minor;

plot(1,win1(1),'ko','Markersize',10,'MarkerFaceColor','k');
plot(33,win2(1),'bo','Markersize',10,'MarkerFaceColor','b');
plot(65,win3(1),'go','Markersize',10,'MarkerFaceColor','g');
plot(97,win4(1),'mo','Markersize',10,'MarkerFaceColor','m');

% mean Theta at each frame
winMat = nan(4,192); 
winMat(1,1:96) = win1;
winMat(2,1+32:96+32) = win2;
winMat(3,1+64:96+64) = win3;
winMat(4,1+96:96+96) = win4;

winMean = wrapTo2Pi( mean( unwrap(winMat,[],1) ,1,'omitnan') ); % phase unwrap
% winMean = mean( winMat,1,'omitnan' );
plot(winMean,'ok','MarkerFaceColor','k');
% axis([65 96 0 7]);

legend('Window 1','Window 2','Window 3','Window 4');





%% Formalise

% Below: Frame offset = mean of two previous windows

for iStk = 4
% for iStk = 1:nStack
    
    SWw = P(iStk).Sweep.swpWinFullWidth;
    SWs = P(iStk).Sweep.swpWinStride;
    nSW = P(iStk).Sweep.numSwpWindows;
    SWo = zeros(1,nSW);
    nSL = max( P(iStk).Sweep.swpWindows(:) );
    
    % Window Position subfn
    wPos = @(iSW) (iSW-1)*SWs+1:(iSW-1)*SWs+SWw;
    
    thetaMat = nan(nSW,nSL);    
    
    % Initialise Windows 1 & 2
    iSW = 1;
    S(iStk).thetaFrameOffset{1,iSW} = wrapTo2Pi( S(iStk).thetaFrame{1,iSW} + SWo(iSW) );
    thetaMat(iSW,wPos(iSW)) = S(iStk).thetaFrameOffset{1,iSW};
    
    iSW = 2;
    SWo(iSW) = S(iStk).thetaFrameOffset{1,iSW-1}(SWs+1);
    S(iStk).thetaFrameOffset{1,iSW} = wrapTo2Pi( S(iStk).thetaFrame{1,iSW} + SWo(iSW) );
    thetaMat(iSW,wPos(iSW)) = S(iStk).thetaFrameOffset{1,iSW};
    
    for iSW = 3:nSW
        
        % Offset = mean of previous two Windows
        SWo(iSW) = mean( [S(iStk).thetaFrameOffset{1,iSW-2}(2*SWs+1), ...
                          S(iStk).thetaFrameOffset{1,iSW-1}(  SWs+1) ]);  
        
        % Adjust Frame Times by Offset
        S(iStk).thetaFrameOffset{1,iSW} = wrapTo2Pi( S(iStk).thetaFrame{1,iSW} + SWo(iSW) );      
        
        % Collate Duplicate Measurements of Theta at Each Frame Location
        thetaMat(iSW,wPos(iSW)) = S(iStk).thetaFrameOffset{1,iSW};       
        
    end
    
end

% thetaMean = mean(thetaMat,1,'omitnan');
thetaMean = wrapTo2Pi( mean( unwrap(thetaMat,[],1) ,1,'omitnan') ); % phase unwrap
% imtar(thetaMat); axis('square'); % not sure working...
figure; hold on; 
plot(thetaMat','o');
plot(thetaMean,'ok','MarkerFaceColor','k');
plot(1:SWs:nSL,thetaMean(1:SWs:nSL),'ko','Markersize',15);
xticks(1:SWs:nSL);
grid; grid minor;
axis([1 nSL -1 2*pi*1.1]);
ylabel('Cardiac Phase (theta)');
xlabel('Sweep Volume Location (z-position/frame index)');



%% Slightly more simplified plot without mean values (but not as useful)
iSW=3;

figure; hold on;
plot( 0*SWs+1:0*SWs+SWw, S(iStk).thetaFrameOffset{1,iSW+0} , 'ro');
plot( 1*SWs+1:1*SWs+SWw, S(iStk).thetaFrameOffset{1,iSW+1} , 'bo');
plot( 2*SWs+1:2*SWs+SWw, S(iStk).thetaFrameOffset{1,iSW+2} , 'go');
plot( 3*SWs+1:3*SWs+SWw, S(iStk).thetaFrameOffset{1,iSW+3} , 'mo');
grid; grid minor;

plot( 0*SWs+1, S(iStk).thetaFrameOffset{1,iSW+0}(1),'ro','Markersize',10,'MarkerFaceColor','r');
plot( 1*SWs+1, S(iStk).thetaFrameOffset{1,iSW+1}(1),'bo','Markersize',10,'MarkerFaceColor','b');
plot( 2*SWs+1, S(iStk).thetaFrameOffset{1,iSW+2}(1),'go','Markersize',10,'MarkerFaceColor','g');
plot( 3*SWs+1, S(iStk).thetaFrameOffset{1,iSW+3}(1),'mo','Markersize',10,'MarkerFaceColor','m');


%% All windows...
figure; hold on;
for ii = 1:nSW
    plot( (ii-1)*SWs+1:(ii-1)*SWs+SWw, S(iStk).thetaFrameOffset{1,ii}   , 'o');
    plot( (ii-1)*SWs+1,                S(iStk).thetaFrameOffset{1,ii}(1), 'k.', 'Markersize',30);
end














