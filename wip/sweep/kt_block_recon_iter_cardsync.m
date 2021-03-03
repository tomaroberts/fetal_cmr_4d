
%% recon iter_cardsync

% TODO:
% - get specific k-t block
% - get cardiac phases
% - compare / do minimisation


%% Admin

swpDir      = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Sweep';
dataDir     = 'data';
ktreconDir  = 'ktrecon';
cardsyncDir = 'cardsync';
cineVolDir  = 'cine_vol';

isVerbose = false;

%% Run: transform_cine_vol_to_stack.bash
% - transform cine_vol to stack space, e.g.:
% - mirtk transform-image cine_vol_3d.nii.gz cine_vol_3d_transf.nii.gz -target s21_rlt_ab_swp3d_apod.nii.gz -interp NN
% - or: mirtk transform-image cine_vol.nii.gz cine_vol_transf.nii.gz -target s21_rlt_ab_swp3d_apod.nii.gz -interp NN


%% Load Data
cd('C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Sweep\c_fcmr334_3d_iter_recon_cardsync');

% Image Frames
I = load_untouch_nii( 's21_rlt_ab_swp3d_apod.nii.gz' );

% Heart Mask
M = load_untouch_nii( 's21_mask_heart.nii.gz' );
M.img = single(M.img);

% 4D cine_vol
% X = load_untouch_nii( 'cine_vol_transf.nii.gz' );
X = load_untouch_nii( 'cine_vol_transf_Bspline.nii.gz' );
X.img = single(X.img);

% cine_vol background
B.img = X.img(:,:,:,1) > 0;

% S
iStk = 1;
load( 'results_cardsync_sweep_resliced.mat' );
thetaFrame = S(iStk).thetaFrameSwpLoc;



%% Auto Crop to cine_vol x/y dimensions

for iLoc = 1:size( X.img,3 )
    
    if ~any( any( X.img(:,:,iLoc,1) ) )
        isEmptySlice(iLoc) = true;
        
        topRow(iLoc)      = nan;
        bottomRow(iLoc)   = nan;
        leftColumn(iLoc)  = nan;
        rightColumn(iLoc) = nan;
        
    elseif any( any( X.img(:,:,iLoc,1) ) )
        isEmptySlice(iLoc) = false;
        
        [ nonZeroRows, nonZeroColumns ] = find( X.img(:,:,iLoc,1) );
        
        topRow(iLoc)      = min( nonZeroRows(:) );
        bottomRow(iLoc)   = max( nonZeroRows(:) );
        leftColumn(iLoc)  = min( nonZeroColumns(:) );
        rightColumn(iLoc) = max( nonZeroColumns(:) );
        
    end
    
end

cropRegion = [ min(topRow), max(bottomRow), min(leftColumn), max(rightColumn) ];

% Apply Crop
I.img = I.img( cropRegion(1):cropRegion(2), cropRegion(3):cropRegion(4), : );
M.img = M.img( cropRegion(1):cropRegion(2), cropRegion(3):cropRegion(4), : );
X.img = X.img( cropRegion(1):cropRegion(2), cropRegion(3):cropRegion(4), : ,: );
B.img = B.img( cropRegion(1):cropRegion(2), cropRegion(3):cropRegion(4), : );


%% Pad z for consistency with original k-t blocks

I.img = padarray( I.img, [0,0,12], 0, 'both' );
M.img = padarray( M.img, [0,0,12], 0, 'both' );
X.img = padarray( X.img, [0,0,12,0], 0, 'both' );
B.img = padarray( B.img, [0,0,12], 0, 'both' );

thetaFrame = padarray( thetaFrame, 12, NaN, 'both' );


%% thetaFrame
nFrame    = numel(thetaFrame);
nBlocks   = 32;
blockSize = 32;

if isVerbose
    figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    hold on;
    
    plot( thetaFrame, 'ok', 'MarkerFaceColor','k' );
    
    xticks( 1:blockSize:nFrame );
    grid; grid minor;
    axis([1 nFrame+10 -0.1 2*pi+0.1]);
    
    title( [ S(iStk).desc ] );
    ylabel('Cardiac Phase (theta)');
    xlabel('Frame Index');
end


%% View

% Full FOV
% implay_RR( [I.img.*B.img I.img.*B.img.*M.img X.img(:,:,:,1) X.img(:,:,:,1).*M.img] );

% k-t block

locStart  = 513;
locRange  = locStart:locStart+blockSize-1;

if isVerbose
    implay_RR( [I.img(:,:,locRange).*B.img(:,:,locRange) ...
        I.img(:,:,locRange).*B.img(:,:,locRange).*M.img(:,:,locRange) ...
        X.img(:,:,locRange,1) ...
        X.img(:,:,locRange,1).*M.img(:,:,locRange)] );
end


%% k-t block loop

thetaFrameModelled = nan( size(thetaFrame) );

ktBlockStartIdx = 1:blockSize:nFrame;


for iKT = 1:nBlocks
    
    %% Compare k-t block with corresponding cine_vol
    
    %     ktBlockStart  = 513;
    ktBlockStart  = ktBlockStartIdx( iKT );
    ktBlockRange  = ktBlockStart:ktBlockStart+blockSize-1;
    ktBlockCentre = ktBlockRange( floor(blockSize/2) );
    
    ktBlockthetaFrame = thetaFrame( ktBlockRange );
    ktBlocktRR        = mean( cell2mat( S(iStk).tRRSwpLoc(ktBlockRange) ) );
    
    if isVerbose
        % cine_vol at centre of k-t block
        implay_RR( X.img(:,:,ktBlockCentre,:) );
    end
    
    %% Check for Heart Voxels
    if ~any( M.img(:,:,ktBlockRange) ,'all' ) || any( isnan( thetaFrame(ktBlockRange) ) )
        
        disp( ['Using original thetaFrame values for k-t block ' num2str(iKT) ' of ' num2str(nBlocks) '...'] );
        thetaFrameModelled( ktBlockRange ) = ktBlockthetaFrame;
        
    else
        
        disp( ['Performing cardsync on k-t block ' num2str(iKT) ' of ' num2str(nBlocks) '...'] );
        
        %% Interp cardiac frames to same number as in k-t block
        nX = size(X.img,1); nY = size(X.img,2); nZ = size(X.img,3);
        Xinterp.img = zeros( nX, nY, nZ, blockSize );
        interpPoints = linspace( 1,25,blockSize );
        
        for xx = 1:nX
            for yy = 1:nY
                for tt = 1:nZ
                    sigt = squeeze( X.img(xx,yy,tt,:) );
                    if any(sigt)
                        Xinterp.img(xx,yy,tt,:) = interp1( sigt, interpPoints );
                    elseif ~any(sigt)
                        Xinterp.img(xx,yy,tt,:) = 0;
                    end
                end
            end
        end
        
        if isVerbose
            implay_RR( [I.img(:,:,ktBlockRange).*B.img(:,:,ktBlockRange).*M.img(:,:,ktBlockRange) ...
                squeeze(Xinterp.img(:,:,ktBlockCentre,:)).*M.img(:,:,ktBlockRange)] );
        end
        
        %% Compare Signal Intensities
        
        % Threshold
        Iblock.img = I.img(:,:,ktBlockRange).*B.img(:,:,ktBlockRange).*M.img(:,:,ktBlockRange);
        BP.img     = ~(Iblock.img < prctile( nonzeros( Iblock.img(:) ), 75 ) );
        IBP.img    = Iblock.img .* BP.img;
        
        for iF = 1:numel( ktBlockRange )
            %     Sig_I(iF) = mean( nonzeros( I.img(:,:,ktBlockRange(iF)).*B.img(:,:,ktBlockRange(iF)).*M.img(:,:,ktBlockRange(iF)) ) );
            Sig_I(iF) = mean( nonzeros( IBP.img(:,:,iF) ) );
        end
        
        for iP = 1:size(Xinterp.img,4)
            Sig_X(iP) = mean( nonzeros( Xinterp.img(:,:,ktBlockCentre,iP).*M.img(:,:,ktBlockCentre) ) );
        end
        
        % circshift Sig_X so timepoint 1 = ED
        [ ~, idxEndDiastoleFrame ] = max( Sig_X );
        Sig_X = circshift( Sig_X, -idxEndDiastoleFrame+1 );
        
        % Calculate vector of cardiac phases for cine_vol
        thetaFrameCineVol = linspace( 0,2*pi,blockSize ); % Is 0 == 2pi? Should I do blockSize+1?
        
        % figure('units','normalized','outerposition',[0.2 0.2 0.4 0.4]);
        % subplot(1,2,1);
        % hold on;
        % plot( Sig_I, '-or' );
        % xlabel('k-t block frame number');
        % ylabel('Signal [a.u.]');
        %
        % subplot(1,2,2);
        % hold on;
        % plot( Sig_X, '-ob' );
        % xlabel('Cardiac phase [frame index]');
        % ylabel('Signal [a.u.]');
        
        
        %%% PERFORM MINIMISATION
        kt_block_recon_minim;
        close all;
        
        thetaFrameModelled( ktBlockRange ) = ktBlockthetaFramePerturbed;
        
    end
    
end

%% Plot thetaFrame before/after
if isVerbose
    figure( 'units', 'normalized', 'outerposition', [0.1 0.1 0.8 0.8] );
    hold on;
    
    subplot(2,1,1);
    plot( thetaFrame, 'ok', 'MarkerFaceColor','k' );
    
    xticks( 1:blockSize:nFrame );
    grid; grid minor;
    axis([1 nFrame+10 -0.1 2*pi+0.1]);
    
    title( [ S(iStk).desc ] );
    ylabel('Cardiac Phase (theta)');
    xlabel('Frame Index');
    
    subplot(2,1,2);
    plot( thetaFrameModelled, 'ok', 'MarkerFaceColor','k' );
    
    xticks( 1:blockSize:nFrame );
    grid; grid minor;
    axis([1 nFrame+10 -0.1 2*pi+0.1]);
    
    title( [ S(iStk).desc ' - after cardsync\_cine\_vol' ] );
    ylabel('Cardiac Phase (theta)');
    xlabel('Frame Index');
end


%% BELOW REPLACED BY: kt_block_recon_minim.m

% %% Minimise SoS
% % nb: something not quite right here - wrapTo2Pi(diff((thetaFrameNew)))
% % gives values which are not identical
%
% KTB      = zeros( blockSize, 3 );
% KTB(:,1) = 1:blockSize;
% KTB(:,2) = Sig_I;
% KTB(:,3) = ktBlockthetaFrame;
%
% % Sort by cardiac phase
% KTBsorted = sortrows( KTB, 3 );
%
% % Perform Minimisation
% Sig_I_sorted = KTBsorted(:,2)';
%
% for iF = 1:numel( ktBlockRange )
%
%     Sig_I_circ = circshift( Sig_I_sorted, iF-1 );
%     Diff       = Sig_I_circ - Sig_X;
%     ssd(iF)    = sum(Diff(:).^2);
%
% end
%
% [ssd_min, ssd_min_idx] = min(ssd);
%
% thetaShift    = 2*pi * ((ssd_min_idx-1)/blockSize); % -1 correct?
% thetaFrameNew = wrapTo2Pi(ktBlockthetaFrame + thetaShift);
%
% KTBsorted(:,3) = wrapTo2Pi( KTBsorted(:,3) + thetaShift );
%
% KTBsorted_shifted = sortrows( KTBsorted, 3 );
%
% % KTBsorted_shifted = circshift( KTBsorted, ssd_min_idx-1, 1 );
% %
% % % Assign new theta values to frames
% % % TODO: check that idxED is the same for every slice in cine_vol
% % KTBsorted_shifted(:,4) = thetaFrameCineVol;
% % Sig_I_New              = KTBsorted_shifted(:,2);
% % KTBsorted_shifted      = sortrows( KTBsorted_shifted, 1 );
% % thetaFrameNew          = KTBsorted_shifted(:,4);
%
%
% %% Perturb R-R interval
%
% deltaTheta    = mean( wrapTo2Pi( diff( ktBlockthetaFrame ) ) );
% Perturbations = deltaTheta .* [0.95 0.97 0.98 0.99 1 1.01 1.02 1.03 1.05];
%
% for iP = 1:numel( Perturbations )
%     thetaFramePerturb(:,iP) = thetaFrameNew(1) : Perturbations(iP) : thetaFrameNew(1)+(Perturbations(iP)*(blockSize-1));
% end
% thetaFramePerturb = wrapTo2Pi(thetaFramePerturb);
%
%
% % Interpolate cine_vol signal curve
% ntheta       = 1000;
% theta        = linspace(0,2*pi,ntheta);
% Sig_X_interp = interp1( thetaFrameCineVol, Sig_X, theta );
%
%
% % Sig_I_New = KTBsorted(:,2);
% Sig_I_New = KTBsorted_shifted(:,2);
%
% for iP = 1:numel( Perturbations )
%
%     thetaIpert = thetaFramePerturb(:,iP)';
%
%     % Find interpolated indices of thetaI
%     idxThetaInterp(iP,:) = dsearchn( theta', thetaIpert' );
%
%     % cine_vol interpolated signal
%     Sig_X_Matrix(iP,:) = Sig_X_interp(idxThetaInterp(iP,:));
%
%     % SSD
%     Diff         = Sig_I_New - Sig_X_Matrix(iP,:);
%     ssd_pert(iP) = sum(Diff(:).^2);
%
% end
%
% [ssd_pert_min, ssd_pert_min_idx] = min(ssd_pert);
%
%
% %% Plots
% figure('units','normalized','outerposition',[0.1 0.2 0.8 0.4]); hold on;
%
% subplot(1,6,1);
% hold on;
% plot( Sig_I, '-or' );
% xlabel('k-t block frame number');
% ylabel('Signal [a.u.]');
% title(' k-t Block Frames ');
%
% subplot(1,6,2);
% plot( KTB(:,3), KTB(:,2), 'xr' );
% xlabel('Cardiac phase [radians]');
% ylabel('Signal [a.u.]');
% title(' k-t Block Before Minimisation ');
%
% subplot(1,6,3);
% plot( Sig_X, '-ob' );
% xlabel('Cardiac phase [frame index]');
% ylabel('Signal [a.u.]');
% title(' cine\_vol model ');
%
% subplot(1,6,4);
% plot( ssd, '-ko' );
% title(' SSD ');
% xlabel('Cardiac phase [frame index]');
%
% subplot(1,6,5); hold on;
% plot( KTBsorted(:,3)', KTBsorted(:,2)' ./ mean(KTBsorted(:,2)'), 'rx' );
% % plot( KTBsorted_shifted(:,4)', KTBsorted_shifted(:,2)' ./ mean(KTBsorted_shifted(:,2)'), 'rx' );
% plot( thetaFrameCineVol, Sig_X ./ mean(Sig_X), '-ob' );
% xlabel('New cardiac phase / frame index');
% ylabel('Normalised Signals [a.u.]');
% title(' k-t Block After Minimisation ');
%
% subplot(1,6,6);
% hold on;
% plot( thetaFrameCineVol, Sig_X ./ mean(Sig_X), '-ob' );
% % plot( KTBsorted_shifted(:,4)', KTBsorted_shifted(:,2)' ./ mean(KTBsorted_shifted(:,2)'), 'rx' );
% plot( thetaFramePerturb(:,ssd_pert_min_idx), KTBsorted_shifted(:,2)' ./ mean(KTBsorted_shifted(:,2)'), 'r^');
% % plot( thetaFramePertMin, KTBsorted_shifted(:,2)' ./ mean(KTBsorted_shifted(:,2)'), '.k');
% % plot( thetaFramePertMax, KTBsorted_shifted(:,2)' ./ mean(KTBsorted_shifted(:,2)'), '.k');
% xlabel('New cardiac phase / frame index');
% ylabel('Normalised Signals [a.u.]');
% title(' k-t Block R-R Perturbation ');
%
%
% %% k-t Block thetaFrame Plots
%
% figure('units','normalized','outerposition',[0.2 0.2 0.4 0.4]);
%
% subplot(1,3,1); hold on;
% plot( locRange, ktBlockthetaFrame, 'ok', 'MarkerFaceColor','k' );
% grid; grid minor;
% axis([locRange(1) locRange(end)+1 -0.1 2*pi+0.1]);
% title( 'Original Theta' );
% ylabel('Cardiac Phase (theta)');
% xlabel('Frame Index');
%
% subplot(1,3,2); hold on;
% plot( locRange, thetaFrameNew, 'ok', 'MarkerFaceColor','k' );
% grid; grid minor;
% axis([locRange(1) locRange(end)+1 -0.1 2*pi+0.1]);
% title( 'Shifted Theta' );
% ylabel('Cardiac Phase (theta)');
% xlabel('Frame Index');
%
% subplot(1,3,3); hold on;
% plot( locRange, thetaFrameNew, 'ok', 'MarkerSize', 9);
% plot( locRange, thetaFramePerturb(:,1), 'ok', 'MarkerFaceColor','k' );
% grid; grid minor;
% axis([locRange(1) locRange(end)+1 -0.1 2*pi+0.1]);
% title( 'Shifted Theta + Perturbations' );
% ylabel('Cardiac Phase (theta)');
% xlabel('Frame Index');
