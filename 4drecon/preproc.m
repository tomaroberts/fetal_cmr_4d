function S = preproc( reconDir, hrRange, acqMethod )
%PREPROC  preprocess data for 4D reconstruction.
%
%   S = PREPROC( reconDir ) reads files from reconDir and returns data structure S.

%   JFPvA (joshua.vanamerom@kcl.ac.uk)
%   TAR   (t.roberts@kcl.ac.uk)


%% arg check
% TODO: make string-value pairs
if nargin < 2
    hrRange   = [105,180];
    acqMethod = 'm2d';
end

if nargin < 3
    acqMethod = 'm2d';    
end


%% Init

dataDir     = fullfile( reconDir, 'data' );
maskDir     = fullfile( reconDir, 'mask' );
ktreconDir  = fullfile( reconDir, 'ktrecon' );


%% Identify Data and Get Parameters

switch acqMethod
    
    % Standard k-t Acquisition
    case 'm2d'

        % Identify Dynamic MR Image Series
        rltFileList       = dir( fullfile( dataDir, '*_rlt_ab.nii.gz' ) );

        % Get Number of Stacks
        nStack = numel(rltFileList);

        % Initialise Stack Struct
        S = struct([]);

        % Read Data Into Stack Struct
        for iStk = 1:nStack

            % Identify Files
            S(iStk).desc          = strrep( rltFileList(iStk).name, '_rlt_ab.nii.gz', '' );
            S(iStk).rltAbFile     = fullfile( rltFileList(iStk).folder,  rltFileList(iStk).name  );
            S(iStk).rltReFile     = fullfile( ktreconDir, strrep( rltFileList(iStk).name, 'ab', 're' ) );
            S(iStk).rltImFile     = fullfile( ktreconDir, strrep( rltFileList(iStk).name, 'ab', 'im' ) );
            S(iStk).rltMatFile    = fullfile( ktreconDir, sprintf( '%s_rlt_recon.mat', S(iStk).desc ) );
            S(iStk).rltParamFile  = fullfile( ktreconDir, sprintf( '%s_rlt_parameters.mat', S(iStk).desc ) );
            S(iStk).dcAbFile      = fullfile( dataDir, sprintf( '%s_dc_ab.nii.gz', S(iStk).desc ) );
            S(iStk).slwAbFile     = fullfile( ktreconDir, sprintf( '%s_slw_ab.nii.gz', S(iStk).desc ) );
            S(iStk).trnAbFile     = fullfile( ktreconDir, sprintf( '%s_trn_ab.nii.gz', S(iStk).desc ) );
            S(iStk).maskHeartFile = fullfile( maskDir, sprintf( '%s_mask_heart.nii.gz', S(iStk).desc ) );

            % Load Parameters
            M = matfile( S(iStk).rltParamFile );
            P = M.PARAM;

            % Extract Parameters
            S(iStk).nLoc             = P.Timing.numLoc;
            S(iStk).sliceThickness   = P.Scan.RecVoxelSize(3);

            % Load NIfTI
            R = load_untouch_nii( S(iStk).rltAbFile );
            S(iStk).niiHdr = R.hdr;
            
            % Update NrDyn
            % TODO: decide how best to implement this
            % - do I want to update PARAM in mrecon_kt.m ?
            if P.Encoding.NrDyn(1) ~= size( R.img,4 )
                warning(['Series ' S(iStk).desc ' - Updating parameter: P.Encoding.NrDyn ... ']);
                P.Encoding.NrDyn = [size( R.img,4 ) size( R.img,4 )];
            end

            % Separate slices
            for iLoc = 1:S(iStk).nLoc

                % Dynamic Image Series
                S(iStk).frameDuration   = P.Timing.frameDuration;
                S(iStk).tFrame{iLoc}    = P.Timing.sliceTime(iLoc) + S(iStk).frameDuration * (0:(P.Encoding.NrDyn(1)-1));
                
            end
            
            % Plot tFrame
            numDyn = P.Encoding.NrDyn(1);
            
            figure; hold on;
            plot( 1:numDyn, S(iStk).tFrame{1,1} );
            for iLoc = 2:S(iStk).nLoc
                xRange = ( (iLoc-1) * numDyn ) + 1:...
                    ( (iLoc) * numDyn );
                plot( xRange, S(iStk).tFrame{1,iLoc});
            end
            title(['Stack ID: ' S(iStk).desc]);
            xlabel('Frame Index');
            ylabel('Time');
            
            hFig = gcf; hFig.Name = sprintf( '%s_kernels_v_time', S(iStk).desc ); % TODO: save in /cardsync
            save_figs( dataDir, gcf, dataDir );

        end
        
	% Sweep k-t Acquisition
    case { 'sweep' , 'swp' }
        
        % Identify Sweep MR Image Series
        rltFileList       = dir( fullfile( ktreconDir, '*_rlt_ab_swp3d.nii.gz' ) ); % originally dataDir, ...

        % Get Number of Stacks
        nStack = numel(rltFileList);

        % Initialise Stack Struct
        S = struct([]);

        % Read Data Into Stack Struct
        for iStk = 1:nStack

            % Identify Files
            S(iStk).desc          = strrep( rltFileList(iStk).name, '_rlt_ab_swp3d.nii.gz', '' );
            S(iStk).rltAbFile     = fullfile( rltFileList(iStk).folder,  rltFileList(iStk).name  );
            S(iStk).rltReFile     = fullfile( ktreconDir, strrep( rltFileList(iStk).name, 'ab', 're' ) );
            S(iStk).rltImFile     = fullfile( ktreconDir, strrep( rltFileList(iStk).name, 'ab', 'im' ) );
            S(iStk).rltMatFile    = fullfile( ktreconDir, sprintf( '%s_rlt_recon.mat', S(iStk).desc ) );
            S(iStk).rltParamFile  = fullfile( ktreconDir, sprintf( '%s_rlt_parameters.mat', S(iStk).desc ) );
            
            % TODO: what to do with below?
            S(iStk).dcAbFile      = fullfile( ktreconDir, sprintf( '%s_dc_ab.nii.gz', S(iStk).desc ) ); % originally dataDir, ...
            S(iStk).slwAbFile     = fullfile( ktreconDir, sprintf( '%s_slw_ab.nii.gz', S(iStk).desc ) );
            S(iStk).trnAbFile     = fullfile( ktreconDir, sprintf( '%s_trn_ab.nii.gz', S(iStk).desc ) );
            S(iStk).maskHeartFile = fullfile( maskDir, sprintf( '%s_mask_heart.nii.gz', S(iStk).desc ) );
            %%%
            
            % Load Parameters
            M = matfile( S(iStk).rltParamFile );
            P = M.PARAM;

            % Load NIfTI
            R = load_untouch_nii( S(iStk).rltAbFile );
            S(iStk).niiHdr = R.hdr;
            
            % Extract Parameters
            S(iStk).nLoc             = size( R.img,3 );
            S(iStk).sliceThickness   = P.Scan.RecVoxelSize(3);

            % Sweep Image Series - Linear Frame Times
            S(iStk).frameDuration   = P.Timing.frameDuration;
            S(iStk).tFrame          = P.Timing.tSeriesOffset + S(iStk).frameDuration * (0:(S(iStk).nLoc-1));
            
%             % View from t = 0
%             figure; plot( S(iStk).tFrame-P.Timing.tSeriesOffset );

        end

end


%% Set Up Additional Data

filePath = fullfile( dataDir, 'tgt_stack_no.txt' );
if ~exist( filePath, 'file' )
  if ~exist( 'iStkTgt', 'var')
    iStkTgt = 1;
  end
    fid = fopen( filePath , 'w' );
    fprintf( fid, num2str(iStkTgt) );
    fclose( fid );
end

filePath = fullfile( dataDir, 'slice_thickness.txt' );
fid = fopen( filePath , 'w' );
fprintf( fid, '%g ', [S.sliceThickness] );
fclose( fid );

filePath = fullfile( dataDir, 'force_exclude_stack.txt' );
if ~exist( filePath, 'file' )
    fid = fopen( filePath , 'w' );
    fprintf( fid, '' );
    fclose( fid );
end

filePath = fullfile( dataDir, 'force_exclude_slice.txt' );
if ~exist( filePath, 'file' )
    fid = fopen( filePath , 'w' );
    fprintf( fid, '' );
    fclose( fid );
end

filePath = fullfile( dataDir, 'force_exclude_frame.txt' );
if ~exist( filePath, 'file' )
    fid = fopen( filePath , 'w' );
    fprintf( fid, '' );
    fclose( fid );
end


end  % preproc(...)

