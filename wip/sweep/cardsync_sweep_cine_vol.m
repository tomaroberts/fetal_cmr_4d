function S = cardsync_sweep_cine_vol( S, varargin )
%CARDSYNC_SWEEP_CINE_VOL  use cine_vol to perform cardiac synchronisation
%of k-t SWEEP data
%
%   S = CARDSYNC_SWEEP_CINE_VOL( S ) 
%
%   CARDSYNC_SWEEP_CINE_VOL( ..., 'name', value ) specifies optional input argument. 
%   See code for name-value pairs.
%
%   See also .

%   tar (t.roberts@kcl.ac.uk)


%% Optional Input Argument Default Values

default.reconDir        = pwd;
default.resultsDir      = pwd;
default.ktBlockSize     = 32;
default.isVerbose       = false;


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

add_param_fn(  p, 'ktBlockSize', default.ktBlockSize, ...
        @(x) validateattributes( x, {'numeric'}, {'scalar','nonnegative'}, mfilename ) );

add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, S, varargin{:} );

reconDir        = p.Results.recondir;
resultsDir      = p.Results.resultsdir;
ktBlockSize     = p.Results.ktBlockSize;
isVerbose       = p.Results.verbose;


%% Setup

dataDir     = fullfile( reconDir, 'data' );
maskDir     = fullfile( reconDir, 'mask' );
ktreconDir  = fullfile( reconDir, 'ktrecon' );
cineVolDir  = fullfile( reconDir, 'cine_vol' );
pngfigDir   = fullfile( resultsDir, 'figs_swp2cv', 'png' );
matfigDir   = fullfile( resultsDir, 'figs_swp2cv', 'fig' );
isVerbose   = false;

mkdir(pngfigDir);
mkdir(matfigDir);

nStack      = numel( S );
% nStack = 1;


%% Load Data

for iStk = 1:nStack

    % Load Params
    MAT = matfile( S(iStk).rltParamFile );
    P(iStk) = MAT.PARAM;

    % Load Sweep Data
    R(iStk) = load_untouch_nii( fullfile( dataDir, [ S(iStk).desc '_rlt_swp3d_ab.nii.gz' ] ) );
    
    % Load Masks
%     N = load_untouch_nii( fullfile( maskDir, [ S(iStk).desc '_mask_heart_swp3d.nii.gz' ] ) );
%     M(iStk).img = single( N.img );
    N = load_untouch_nii( fullfile( maskDir, [ S(iStk).desc '_mask_heart.nii.gz' ] ) );
    M(iStk).img = single( reshape( permute( repmat( N.img,[1,1,1,P(iStk).Sweep.swpWinFullWidth] ), [1,2,4,3] ), [size(N.img,1),size(N.img,2),P(iStk).Sweep.numSwpLoca] ) );
    
    % Load Cine Vol
    X(iStk) = load_untouch_nii( fullfile( dataDir, [ 'cine_vol_trans2swp_' S(iStk).desc '.nii.gz' ] ) );
    X(iStk).img = single( X(iStk).img );
    
    % Cine Vol background
    B(iStk).img = X(iStk).img(:,:,:,1) > 0;
    
    clear N MAT
    
end


%% Auto Crop to cine_vol x/y dimensions

for iStk = 1:nStack
    for iLoc = 1:size( X(iStk).img,3 )
        
        if ~any( any( X(iStk).img(:,:,iLoc,1) ) )
            isEmptySlice(iLoc) = true;
            
            topRow(iLoc)      = nan;
            bottomRow(iLoc)   = nan;
            leftColumn(iLoc)  = nan;
            rightColumn(iLoc) = nan;
            
        elseif any( any( X(iStk).img(:,:,iLoc,1) ) )
            isEmptySlice(iLoc) = false;
            
            [ nonZeroRows, nonZeroColumns ] = find( X(iStk).img(:,:,iLoc,1) );
            
            topRow(iLoc)      = min( nonZeroRows(:) );
            bottomRow(iLoc)   = max( nonZeroRows(:) );
            leftColumn(iLoc)  = min( nonZeroColumns(:) );
            rightColumn(iLoc) = max( nonZeroColumns(:) );
            
        end
        
    end
    
    cropRegion = [ min(topRow), max(bottomRow), min(leftColumn), max(rightColumn) ];
    
    % Apply Crop
    R(iStk).img = R(iStk).img( cropRegion(1):cropRegion(2), cropRegion(3):cropRegion(4), : );
    M(iStk).img = M(iStk).img( cropRegion(1):cropRegion(2), cropRegion(3):cropRegion(4), : );
    X(iStk).img = X(iStk).img( cropRegion(1):cropRegion(2), cropRegion(3):cropRegion(4), :, : );
    B(iStk).img = B(iStk).img( cropRegion(1):cropRegion(2), cropRegion(3):cropRegion(4), : );
    
end


%% Pad z for consistency with original k-t blocks

% TODO: might need to re-introduce padding when apodization fully implemented

% I.img = padarray( I.img, [0,0,12], 0, 'both' );
% M.img = padarray( M.img, [0,0,12], 0, 'both' );
% X.img = padarray( X.img, [0,0,12,0], 0, 'both' );
% B.img = padarray( B.img, [0,0,12], 0, 'both' );
% 
% thetaFrame = padarray( thetaFrame, 12, NaN, 'both' );


%% View thetaFrame Plots

for iStk = 1:nStack
    
    nFrame(iStk)  = P(iStk).Sweep.numSwpLoca;
    nBlocks(iStk) = nFrame(iStk) / ktBlockSize;
    % ktBlockSize = 32;

    figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    hold on;

    plot( cell2mat( S(iStk).thetaFrame ), 'ok', 'MarkerFaceColor','k' );

    xticks( 1:ktBlockSize:nFrame(iStk) );
    grid; grid minor;
    axis([1 nFrame(iStk)+10 -0.1 2*pi+0.1]);

    title( [ S(iStk).desc ' - cardsync\_interslice'] );
    ylabel('Cardiac Phase (theta)');
    xlabel('Frame Index');

    % Save Fig
    hFig = gcf;
    hFig.Name = sprintf( '%s_thetaFrame_cardsync_interslice', S(iStk).desc );
    save_figs( pngfigDir, hFig, matfigDir )

    if ~isVerbose
        close( hFig )
    end

end


%% View cine vol at k-t Block location
if isVerbose    
    iStk1 = 1;
    
    locStart  = 193;
    locRange  = locStart:locStart+ktBlockSize-1;
    
    implay_RR( [R(iStk1).img(:,:,locRange).*B(iStk1).img(:,:,locRange) ...
        R(iStk1).img(:,:,locRange).*B(iStk1).img(:,:,locRange).*M(iStk1).img(:,:,locRange) ...
        X(iStk1).img(:,:,locRange,1) ...
        X(iStk1).img(:,:,locRange,1).*M(iStk1).img(:,:,locRange)] );
end


%% Stack-Cine Vol Cardiac Synchronisation

% Loop over stacks
for iStk = 1:nStack
% for iStk = 4
    
    disp( ['Performing Stack-Cine Vol synchronisation on stack no.: ' num2str(iStk) ] );

    % Loop over k-t blocks
    thetaFrameModelled = nan( size( cell2mat( S(iStk).thetaFrame ) ) );    
    ktBlockStartIdx    = 1 : ktBlockSize : nFrame(iStk);
    thetaFrame         = cell2mat( S(iStk).thetaFrame );
    tRR                = repelem( cell2mat( S(iStk).tRR ), P(iStk).Sweep.swpWinFullWidth );
    
    for iKT = 1:nBlocks(iStk)
%     for iKT = 3
        
        %% Compare k-t block with corresponding cine_vol
        % ktBlockStart  = 513;
        ktBlockStart  = ktBlockStartIdx( iKT );
        ktBlockRange  = ktBlockStart:ktBlockStart+ktBlockSize-1;
        ktBlockCentre = ktBlockRange( floor(ktBlockSize/2) );
        
        ktBlockthetaFrame = thetaFrame( ktBlockRange );
        ktBlocktRR        = mean( tRR(ktBlockRange) );
        
        if isVerbose
            % cine_vol at centre of k-t block
            implay_RR( X(iStk).img(:,:,ktBlockCentre,:) );
        end
        
        %% Check for Heart Voxels
        if ~any( M(iStk).img(:,:,ktBlockRange) , 'all' )    || ...
           ~any( X(iStk).img(:,:,ktBlockCentre,:) , 'all' ) || ...
            any( isnan( thetaFrame(ktBlockRange) ) )
            
            disp( ['Using original thetaFrame values for k-t block ' num2str(iKT) ' of ' num2str(nBlocks) '...'] );
            thetaFrameModelled( ktBlockRange ) = ktBlockthetaFrame;
            
        else
            
            disp( ['Performing cardsync on k-t block ' num2str( iKT ) ' of ' num2str( nBlocks(iStk) ) '...'] );
            
            %% Interp cardiac frames to same number as in k-t block
            nX = size(X(iStk).img,1); nY = size(X(iStk).img,2); nZ = size(X(iStk).img,3);
            Xinterp.img = zeros( nX, nY, nZ, ktBlockSize );
            interpPoints = linspace( 1,25,ktBlockSize );
            
            for xx = 1:nX
                for yy = 1:nY
                    for tt = 1:nZ
                        sigt = squeeze( X(iStk).img(xx,yy,tt,:) );
                        if any(sigt)
                            Xinterp.img(xx,yy,tt,:) = interp1( sigt, interpPoints );
                        elseif ~any(sigt)
                            Xinterp.img(xx,yy,tt,:) = 0;
                        end
                    end
                end
            end
                      
            
            %% Get Mean Signals from Thresholded Blood Pool
            Rblock.img = R(iStk).img(:,:,ktBlockRange) .* B(iStk).img(:,:,ktBlockRange) .* M(iStk).img(:,:,ktBlockRange);
            BP.img     = ~( Rblock.img < prctile( nonzeros( Rblock.img(:) ), 75 ) );
            RBP.img    = Rblock.img .* BP.img;
            
%             Xblock.img = X(iStk).img(:,:,ktBlockRange) .* B(iStk).img(:,:,ktBlockRange) .* M(iStk).img(:,:,ktBlockRange);
%             % BP.img     = ~( Xblock.img < prctile( nonzeros( Xblock.img(:) ), 75 ) );
%             XBP.img    = Xblock.img .* BP.img;
            Xblock.img = squeeze(Xinterp.img(:,:,ktBlockCentre,:)) .* M(iStk).img(:,:,ktBlockRange);

            % Video
            if isVerbose
                implay_RR( [ Rblock.img, ...
                            RBP.img, ... 
                            squeeze(Xinterp.img(:,:,ktBlockCentre,:)) .* M(iStk).img(:,:,ktBlockRange) ...
                            ] );
            end
            
            % k-t block Signal
            for iF = 1:numel( ktBlockRange )
                Sig_R(iF) = mean( nonzeros( Rblock.img(:,:,iF) ) );
%                 Sig_R(iF) = mean( nonzeros( RBP.img(:,:,iF) ) );
            end
            
            % Cine Vol Signal
            for iP = 1:size(Xinterp.img,4)
                Sig_X(iP) = mean( nonzeros( Xinterp.img(:,:,ktBlockCentre,iP).*M(iStk).img(:,:,ktBlockCentre) ) );
%                 Sig_X(iP) = mean( nonzeros( Xinterp.img(:,:,ktBlockRange(iF),iP).*M(iStk).img(:,:,ktBlockRange(iF)) ) );
%                 Sig_X(iP) = mean( nonzeros( XBP.img(:,:,iP) ) );
            end
            
            % circshift Sig_X so timepoint 1 = ED
            [ ~, idxEndDiastoleFrame ] = max( Sig_X );
            Sig_X = circshift( Sig_X, -idxEndDiastoleFrame+1 );
            
            % Calculate vector of cardiac phases for cine_vol
            thetaFrameCineVol = linspace( 0,2*pi,ktBlockSize ); % Is 0 == 2pi? Should I do ktBlockSize+1?          
            
            %%% PERFORM MINIMISATION
            kt_block_recon_minim;
            close all;
            
            thetaFrameModelled( ktBlockRange ) = ktBlockthetaFrameShifted; % circ shifted only
%             thetaFrameModelled( ktBlockRange ) = ktBlockthetaFramePerturbed; % circshift + perturbed
            
            
        end
    end
    
    % TODO: reformat to same shape/structure as S.thetaFrame
    S(iStk).thetaFrameModelled = thetaFrameModelled;
    S(iStk).thetaFrameModelled = num2cell( reshape( S(iStk).thetaFrameModelled, P(iStk).Sweep.swpWinFullWidth, P(iStk).Sweep.numSwpWindows ) , [1,P(iStk).Sweep.numSwpWindows] );
    
end


%% Plot thetaFrame before/after
if isVerbose
    
    for iStk = 1:nStack
        figure( 'units', 'normalized', 'outerposition', [0.1 0.1 0.8 0.8] );
        hold on;
        
        subplot(2,1,1);
        plot( S(iStk).thetaFrame, 'ok', 'MarkerFaceColor','k' );
        
        xticks( 1:ktBlockSize:nFrame(iStk) );
        grid; grid minor;
        axis([1 nFrame(iStk)+10 -0.1 2*pi+0.1]);
        
        title( [ S(iStk).desc ] );
        ylabel('Cardiac Phase (theta)');
        xlabel('Frame Index');
        
        subplot(2,1,2);
        plot( S(iStk).thetaFrameModelled, 'ok', 'MarkerFaceColor','k' );
        
        xticks( 1:ktBlockSize:nFrame(iStk) );
        grid; grid minor;
        axis([1 nFrame(iStk)+10 -0.1 2*pi+0.1]);
        
        title( [ S(iStk).desc ' - after cardsync\_cine\_vol' ] );
        ylabel('Cardiac Phase (theta)');
        xlabel('Frame Index');
    end
    
end


%% Save Results to Text Files

% TODO: finalise which I need

% fid = fopen( fullfile( resultsDir, 'mean_rrinterval.txt' ), 'w' );
% fprintf( fid, '%.6f ', mean( cell2mat( [ S.tRR ] ) ) );
% fclose( fid );
% fid = fopen( fullfile( resultsDir, 'rrintervals.txt' ), 'w' );
% fprintf( fid, '%.6f ', cell2mat( [ S.tRR ] ) );
% fclose( fid );
fid = fopen( fullfile( resultsDir, 'cardphases_cine_vol_cardsync.txt' ), 'w' );
fprintf( fid, '%.6f ', cell2mat( [S.thetaFrameModelled] ) ); % TODO: check
fclose( fid );


%% Save Results to .mat File

% Tidy S structure
for iStk = 1:nStack
    S(iStk).thetaFrame = S(iStk).thetaFrameModelled;
end
S = rmfield( S, 'thetaFrameModelled' );

save( fullfile( resultsDir, 'results_cardsync_sweep_cine_vol.mat' ), 'S', '-v7.3' );



end % cardsync_sweep_cine_vol(...)