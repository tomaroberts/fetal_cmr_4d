function S = cardsync_sweep_interstack( S, varargin )
%CARDSYNC_SWEEP_INTERSTACK  align cardiac phases of Sweep volumes.
%
%   S = CARDSYNC_SWEEP_INTERSTACK( S ) 
%
%   CARDSYNC_SWEEP_INTERSTACK( ..., 'name', value ) specifies optional input argument. 
%   See code for name-value pairs.
%
%   See also PREPROC.

%   tar (t.roberts@kcl.ac.uk)


%% Optional Input Argument Default Values

default.reconDir        = pwd;
default.resultsDir      = pwd;
default.nPhases         = 25;       % Number of Phases for determination of Stack-stack alignment
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

add_param_fn(   p, 'recondir', default.reconDir, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'resultsdir', default.resultsDir, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(  p, 'nphases', default.nPhases, ...
        @(x) validateattributes( x, {'numeric'}, {'scalar','nonnegative'}, mfilename ) );

add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, S, varargin{:} );

reconDir        = p.Results.recondir;
resultsDir      = p.Results.resultsdir;
nPhases         = p.Results.nphases;
isVerbose       = p.Results.verbose;


%% Setup

dataDir     = fullfile( reconDir, 'data' );
maskDir     = fullfile( reconDir, 'mask' );
ktreconDir  = fullfile( reconDir, 'ktrecon' );
pngfigDir   = fullfile( resultsDir, 'figs', 'png' );
matfigDir   = fullfile( resultsDir, 'figs', 'fig' );
isVerbose   = false;

nStack      = numel( S );


%% Load Data

for iStk = 1:nStack

    % Load Params
    MAT = matfile( S(iStk).rltParamFile );
    P(iStk) = MAT.PARAM;

    % Load Sweep Data
    R(iStk) = load_untouch_nii( S(iStk).rltAbSwpApodFile );
    
    % Load Masks
    N = load_untouch_nii( S(iStk).maskHeartSwpApodFile );
    M{iStk} = single( N.img );

    clear N MAT
    
end


%% Blood-pool Threshold

for iStk = 1:nStack

    % Apply Mask
    RM{iStk} = R(iStk).img .* M{iStk};
    
    % Blood-pool Threshold
    BP{iStk} = ~(RM{iStk} < prctile( nonzeros( RM{iStk}(:) ), 75 ) );
    
    % Apply Threshold
    RBP{iStk} = R(iStk).img .* BP{iStk};

end

% iS = 530;
% imtar(RM{iStk}(:,:,iS));
% imtar(BP{iStk}(:,:,iS));
% imtar(RBP{iStk}(:,:,iS));


%% Calculate Mean Intensity in Blood-pool

for iStk = 1:nStack
    
    numSwpLoc = size( RBP{iStk}, 3 );
    
    for iLoc = 1:numSwpLoc     
        
        if ~any( M{iStk}(:,:,iLoc), 'all' )
            S(iStk).meanSignalBP(iLoc) = NaN;
        else
            S(iStk).meanSignalBP(iLoc) = mean( nonzeros( RBP{iStk}(:,:,iLoc) ) );
        end
                 
    end    

end

% figure; plot(S(iStk).meanSignalBP,'.-k');
% AX = axis; axis([1 numSwpLoc AX(3) AX(4)]);


%% Bin Cardiac Phases to Estimate Cardiac Volume Curve

for iStk = 1:nStack

    thetaFrameSwpBins        = [];
    thetaFrameSwpBins        = cell2mat( S(iStk).thetaFrameSwpBins );
    S(iStk).thetaFrameSwpLoc = thetaFrameSwpBins(:);
    
    % Apply Apodization
    S(iStk).thetaFrameSwpLoc( S(iStk).isApod == 1) = [];
    
    % Setup Cardiac Phase Bins
    edgesPhase = []; binsPhases = [];
    edgesPhases = linspace(0,2*pi,nPhases+1);
    binsPhases  = discretize( S(iStk).thetaFrameSwpLoc', edgesPhases );

    % Mean Blood-pool Signal over Frames in Each Phase Range   
    for iP = 1:nPhases
        
        framesInCurrentPhaseRange = find( binsPhases == iP );
        
        meanSignalIntensity{iStk}(iP) = mean( nonzeros( RBP{iStk}(:,:,framesInCurrentPhaseRange) ) );
                                            
        clear framesInCurrentPhaseRange
        
    end
       
    % Plot
    figure; hold on;
    plot( edgesPhases(1:end-1), meanSignalIntensity{iStk}, '.--k' );
    plot( edgesPhases(1:end-1), smooth( meanSignalIntensity{iStk}, .2, 'moving' ), '.-b');
    plot( edgesPhases(1:end-1), smooth( meanSignalIntensity{iStk}, .2, 'sgolay' ), '.-r');
    
    [valMax,iMax] = max( meanSignalIntensity{iStk} );
    [valMin,iMin] = min( meanSignalIntensity{iStk} );
    plot( [edgesPhases(iMin), edgesPhases(iMax)], [valMin valMax], 'ok', 'MarkerSize', 10 );
    
    [valMax,iMax] = max( smooth( meanSignalIntensity{iStk}, .2, 'moving' ) );
    [valMin,iMin] = min( smooth( meanSignalIntensity{iStk}, .2, 'moving' ) );
    plot( [edgesPhases(iMin), edgesPhases(iMax)], [valMin valMax], 'ob', 'MarkerSize', 10 );
    
    [valMax,iMax] = max( smooth( meanSignalIntensity{iStk}, .2, 'sgolay' ) );
    [valMin,iMin] = min( smooth( meanSignalIntensity{iStk}, .2, 'sgolay' ) );
    plot( [edgesPhases(iMin), edgesPhases(iMax)], [valMin valMax], 'or', 'MarkerSize', 10 );
    
    grid; grid minor;
    xlabel('Cardiac Phase [radians]');
    ylabel('Mean Signal Intensity [a.u.]');
    title([S(iStk).desc ' -  RR Interval - Intra-stack']);
    legend('Mean','Moving Av.','S-Golay');
    
    % Save Fig
    hFig = gcf;
    hFig.Name = sprintf( 'sweep_cardiac_phase_intrastack_%s', S(iStk).desc );
    save_figs( pngfigDir, hFig, matfigDir )
    
    if ~isVerbose
        close( hFig )
    end
    
end


%% Align Stacks - Cross Correlation
% TODO: implement cross correlation of stack R-R intervals
% for iStk = 1:nStack
%    
%     xcorr ... ?
%     
% end


%% Align End-Diastole of Each Stack

% Alignment Method:
smoothingType = 'moving';
smoothingSpan = 0.2;        % percentage of total points, i.e.: 20%

for iStk = 1:nStack

    [~, edIdx] = max( smooth( meanSignalIntensity{iStk}, smoothingSpan, smoothingType ) );
    
    S(iStk).thetaInterStackOffset = edgesPhases( edIdx );
    
    meanSignalIntensity{iStk} = circshift( meanSignalIntensity{iStk}, -edIdx );
    
    figure; hold on;
    plot( edgesPhases(1:end-1), meanSignalIntensity{iStk}, '.--k' );
    plot( edgesPhases(1:end-1), smooth( meanSignalIntensity{iStk}, .2, 'moving' ), '.-b');
    plot( edgesPhases(1:end-1), smooth( meanSignalIntensity{iStk}, .2, 'sgolay' ), '.-r');
    
    grid; grid minor;
    xlabel('Cardiac Phase [radians]');
    ylabel('Mean Signal Intensity [a.u.]');
    title([S(iStk).desc ' -  RR Interval - Stack-stack aligned']);
    legend('Mean','Moving Av.','S-Golay','Location','SouthEast');
    
    % Save Fig
    hFig = gcf;
    hFig.Name = sprintf( 'sweep_cardiac_phase_interstack_aligned_%s', S(iStk).desc );
    save_figs( pngfigDir, hFig, matfigDir )
    
    if ~isVerbose
        close( hFig )
    end
    
end


%% Update thetaFrame Timings
% TODO: save updated version of sweep_intrastack_cardiac_phases_s*.fig ?
%       i.e.: sweep_interstack_cardiac_phases_s*.fig


for iStk = 1:nStack
    
    S(iStk).thetaFrameSwpLoc = wrapTo2Pi( S(iStk).thetaFrameSwpLoc - S(iStk).thetaInterStackOffset );
    
end


%% Save Results to Text Files

% TODO: finalise which I need

% fid = fopen( fullfile( resultsDir, 'mean_rrinterval.txt' ), 'w' );
% fprintf( fid, '%.6f ', mean( cell2mat( [ S.tRR ] ) ) );
% fclose( fid );
% fid = fopen( fullfile( resultsDir, 'rrintervals.txt' ), 'w' );
% fprintf( fid, '%.6f ', cell2mat( [ S.tRR ] ) );
% fclose( fid );
fid = fopen( fullfile( resultsDir, 'cardphases_interstack_cardsync.txt' ), 'w' );
fprintf( fid, '%.6f ', cat( 1, S.thetaFrameSwpLoc ) );
fclose( fid );


%% Save Results to .mat File

save( fullfile( resultsDir, 'results_cardsync_sweep_interstack.mat' ), 'S', '-v7.3' );


% cardsync_sweep_interstack(...)
end
