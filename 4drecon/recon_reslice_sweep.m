

%% Run recon_binned_sweep.m until apodization is complete

% Then ...


%% Reslice Sweep Volumes into 4D

for iStk = 1:nStack
% for iStk = 1
    
    % Get Cardiac Phases
    thetaFrameSwpBins = [];
    thetaFrameSwpBins = cell2mat( S(iStk).thetaFrameSwpBins );
    thetaFrameSwpBins = thetaFrameSwpBins(:);
    
    
    % Binning Configuration
    nX         = size( R(iStk).img, 1 );
    nY         = size( R(iStk).img, 2 );
    numSwpLoca = max( P(iStk).Sweep.swpWindows(:) ); %TODO: change to size( R(iStk).img, 3 ); so compatible with sweep_window_filter.m ?
    nDyn       = 64;
    nSlices    = numSwpLoca / nDyn;
%    nSlices    = 16;
    
    binWidthSlices = numSwpLoca / nSlices;

    if ~( isreal( binWidthSlices ) && rem( binWidthSlices ,1)==0 )
        error( ['Number of slices does not give integer bin width. Possible bin widths = ' num2str(divisors(numSwpLoca)) ] );
    end
    
    % Init Slice and Cardiac Phase Bins
    edgesSlices = 1:binWidthSlices:numSwpLoca+binWidthSlices;
    binsSlices  = discretize( 1:numSwpLoca, edgesSlices );
    
    % Bin Sweep volume into Slices
%     R_binned{iStk,1} = cell(1,nSlices);
    R_binned(iStk).img = [];
    
    for iS = 1:nSlices

        currentBinRange = edgesSlices(iS):edgesSlices(iS+1)-1;
        
%         R_binned{iStk,1}{1,iS} = R(iStk).img( :,:,currentBinRange );
    R_binned(iStk).img( :,:,iS,: ) = R(iStk).img( :,:,currentBinRange );

    end
    
end


%% Save NIfTI

for iStk = 1:nStack
    
    % Update img
    R(iStk).img = R_binned(iStk).img;
    
    % Update Header
    % TODO: do I need to update more fields? Affine?
    nSlices = size( R_binned(iStk).img, 3 );
    nDyn    = size( R_binned(iStk).img, 4 );
    R(iStk).hdr.dime.dim([4,5]) = [nSlices, nDyn];     
    
    % Save
    S(iStk).rltBinnedAbFile = fullfile( ktreconDir, sprintf( '%s_rlt_ab_swp_sliced.nii.gz', S(iStk).desc ) );
    save_untouch_nii( R(iStk), S(iStk).rltBinnedAbFile );

end

