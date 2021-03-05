% Run: C:\Users\tr17\Documents\CODE_Projects\fetal_cmr_4d\wip\sweep\kt_block_recon_iter_cardsync.m

% Then:

%% Minimisation 1): Circshift
% - align frames with cine_vol model
KTB      = zeros( ktBlockSize, 3 );
KTB(:,1) = 1:ktBlockSize;             % thetaFrame indices
KTB(:,2) = Sig_R;                     % frame Signal intensities
KTB(:,3) = ktBlockthetaFrame;         % thetaFrame

% Sort by cardiac phase
KTB = sortrows( KTB, 3 );
Sig_R_sorted = KTB(:,2)';
ssd = [];

% Perform Minimisation
for iF = 1:numel( ktBlockRange )
    
    Sig_R_circ = circshift( Sig_R_sorted, iF-1 );
    Diff       = Sig_R_circ - Sig_X;
    ssd(iF)    = sum(Diff(:).^2);
    
end

[ssd_min, ssd_min_idx] = min( ssd );

% Shift thetaFrame
thetaShift = 2*pi * ((ssd_min_idx-1)/ktBlockSize); % -1 correct?
KTB        = sortrows( KTB, 1 ); % resort by frame order
KTB(:,3)   = wrapTo2Pi( KTB(:,3) + thetaShift );

ktBlockthetaFrameShifted = KTB(:,3);


%% Minimisation 2): R-R interval Shift
% - apply constant offsets to difference between thetaFrame values

% Setup interpolated cine_vol signal curve
ntheta       = 1000;
theta        = linspace(0,2*pi,ntheta);
Sig_X_interp = interp1( thetaFrameCineVol, Sig_X, theta );

% Set R-R interval Perturbations
deltaTheta    = mean( wrapTo2Pi( diff( ktBlockthetaFrame ) ) );
Perturbations = deltaTheta .* [0.5 0.95 0.97 0.98 0.99 1 1.01 1.02 1.03 1.05 1.5];
% Perturbations = deltaTheta .* [0.5:0.05:1.5];
% Perturbations = deltaTheta .* [0.98:0.005:1.02];

% Extend KTB
KTB(:,1:2,1:numel(Perturbations)) = repmat( KTB(:,1:2), [1,1,numel( Perturbations)] );

% Apply Perturbations to thetaFrame
for iP = 1:numel( Perturbations )
    KTB(:,3,iP) = ktBlockthetaFrameShifted(1) : Perturbations(iP) : ktBlockthetaFrameShifted(1)+(Perturbations(iP)*(ktBlockSize-1));
end
KTB(:,3,:) = wrapTo2Pi( KTB(:,3,:) );

% Plot Perturbed thetaFrame
if isVerbose
    figure;
    for iP = 1:numel( Perturbations )
        subplot(4,3,iP); hold on;
        plot( ktBlockRange, KTB(:,3,iP), 'ok', 'MarkerFaceColor','k' );
        grid; grid minor;
        axis([ktBlockRange(1) ktBlockRange(end)+1 -0.1 2*pi+0.1]);
        title(['Perturbation = ' num2str(deltaTheta - Perturbations(iP))]);
    end
end

% Perform Minimisation
ssd_pert = [];

for iP = 1:numel( Perturbations )
    % for iP = 5
    
    KTBtemp      = sortrows( KTB(:,:,iP), 3 );
    Sig_R_sorted = KTBtemp(:,2)';
    thetaR       = KTBtemp(:,3)';
    
    Sig_R_Matrix(iP,:) = Sig_R_sorted ./ mean( Sig_R_sorted );
    
    % Find interpolated indices of thetaI
    idxThetaInterp(iP,:) = dsearchn( theta', thetaR' );
    
    % cine_vol interpolated signal
    Sig_X_Matrix(iP,:) = Sig_X_interp( idxThetaInterp(iP,:) ) ./ mean( Sig_X_interp( idxThetaInterp(iP,:) ) );
    
    % SSD
    Diff(iP,:)   = Sig_R_Matrix(iP,:) - Sig_X_Matrix(iP,:);
    ssd_pert(iP) = sum(Diff(iP,:).^2);
    
end

[ssd_pert_min, ssd_pert_min_idx] = min( ssd_pert );

% Shift thetaFrame
KTBtemp = sortrows( KTB(:,:,ssd_pert_min_idx), 1 ); % resort by frame order
ktBlockthetaFramePerturbed = KTBtemp(:,3);


%% Plots
% if isVerbose
    
    %% k-t Block step-by-step
    iP0 = find( Perturbations./deltaTheta == 1 ); % no perturbation
    
    figure('units','normalized','outerposition',[0.2 0.05 0.6 0.9]);

    subplot(2,3,1);
    hold on;
    plot( Sig_R, '-or' );
    xlabel('k-t block frame number');
    ylabel('Signal [a.u.]');
    title([' k-t Block Signal (start = ' num2str(ktBlockStart) ')' ]);
    
    subplot(2,3,2);
    plot( Sig_X, '-ob' );
    xlabel('Cardiac phase [frame index]');
    ylabel('Signal [a.u.]');
    title(' cine\_vol Signal ');

    subplot(2,3,3);
    plot( ktBlockthetaFrame, Sig_R, 'xr' );
    xlabel('Cardiac phase [radians]');
    ylabel('Signal [a.u.]');
    title(' k-t Block - Before Minimisation ');

    subplot(2,3,4);
    plot( ssd, '-ko' );
    xlabel('k-t block frame number');
    ylabel('SSD [a.u.]');
    title('SSD - Circshift Frames to Cine Vol');
    
    subplot(2,3,5);
    plot( Perturbations, ssd_pert, '-ko' );
    xlabel('deltaTheta Offset [radians]');
    ylabel('SSD [a.u.]');
    title('SSD - RR Interval Perturbations');
    
    subplot(2,3,6); hold on;
    plot( Sig_X_interp ./ mean(Sig_X_interp ),'b');
    plot( idxThetaInterp(ssd_pert_min_idx,:), Sig_X_Matrix(ssd_pert_min_idx,:), 'bo');
    plot( idxThetaInterp(ssd_pert_min_idx,:), Sig_R_Matrix(ssd_pert_min_idx,:), 'rx');
    xlabel('Cardiac phase [radians]');
    ylabel('Signal [a.u.]');
    title(' k-t Block - After Minimisation ');
    
    % Save Fig
    hFig = gcf;
    hFig.Name = sprintf( '%s_ktblock%s_minim_overview', S(iStk).desc, num2str(ktBlockStart) );
    save_figs( pngfigDir, hFig, matfigDir )
    
    
    %% Signals after Perturbations
    figure('units','normalized','outerposition',[0.2 0.2 0.4 0.7]);
    for iP = 1:numel( Perturbations )
        subplot(4,3,iP); hold on;
        plot( thetaR, Sig_X_Matrix(iP,:), 'b' );
        plot( thetaR, Sig_R_Matrix(iP,:), 'rx' );
        title(['Perturbation = ' num2str(deltaTheta - Perturbations(iP))]);
    end
    
    % Save Fig
    hFig = gcf;
    hFig.Name = sprintf( '%s_ktblock%s_perturbation_signals', S(iStk).desc, num2str(ktBlockStart) );
    save_figs( pngfigDir, hFig, matfigDir )
    
    
    %% Residuals after Perturbations
    figure('units','normalized','outerposition',[0.2 0.2 0.4 0.7]);
    for iP = 1:numel( Perturbations )
        subplot(4,3,iP); hold on;
        plot( thetaR,  Diff(iP,:), 'ko' );
        title(['Perturbation = ' num2str(deltaTheta - Perturbations(iP))]);
        axis([0 2*pi min(Diff(:)) max(Diff(:))])
    end
    
    % Save Fig
    hFig = gcf;
    hFig.Name = sprintf( '%s_ktblock%s_perturbation_residuals', S(iStk).desc, num2str(ktBlockStart) );
    save_figs( pngfigDir, hFig, matfigDir )


    
    %% k-t Block thetaFrame Before/After
    figure('units','normalized','outerposition',[0.2 0.2 0.4 0.4]);
    
    subplot(1,3,1); hold on;
    plot( ktBlockRange, ktBlockthetaFrame, 'ok', 'MarkerFaceColor','k' );
    grid; grid minor;
    axis([ktBlockRange(1) ktBlockRange(end)+1 -0.1 2*pi+0.1]);
    title( 'Original Theta' );
    ylabel('Cardiac Phase (theta)');
    xlabel('Frame Index');
    
    subplot(1,3,2); hold on;
    plot( ktBlockRange, ktBlockthetaFrameShifted, 'ok', 'MarkerFaceColor','k' );
    grid; grid minor;
    axis([ktBlockRange(1) ktBlockRange(end)+1 -0.1 2*pi+0.1]);
    title( 'Shifted Theta' );
    ylabel('Cardiac Phase (theta)');
    xlabel('Frame Index');
    
    subplot(1,3,3); hold on;
    plot( ktBlockRange, ktBlockthetaFrameShifted, 'ok', 'MarkerSize', 9);
    plot( ktBlockRange, ktBlockthetaFramePerturbed, 'ok', 'MarkerFaceColor','k' );
    grid; grid minor;
    axis([ktBlockRange(1) ktBlockRange(end)+1 -0.1 2*pi+0.1]);
    title( 'Shifted Theta + Perturbations' );
    ylabel('Cardiac Phase (theta)');
    xlabel('Frame Index');
    
    % Save Fig
    hFig = gcf;
    hFig.Name = sprintf( '%s_ktblock%s_thetaframe', S(iStk).desc, num2str(ktBlockStart) );
    save_figs( pngfigDir, hFig, matfigDir )
    
    if ~isVerbose
        close all;
    end

% end