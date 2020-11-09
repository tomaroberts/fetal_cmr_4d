function fcmr_4d_recon_beastie01( varargin )

%FCMR_4D_RECON_BEASTIE01  single script to perform fetal_cmr_4d
%   reconstruction on beastie01 in Perinatal @ KCL
%
%   This is tailored to beastie01:
%   - Requires user cmo19
%   - Requires MRecon (v557, NOT v515) for ktrecon from .raw
%   - Requires access to .raw scan database
%   - Requires installed fetal_cmr_4d repo
%   - Requires MIRTK / SVRTK installed
%   - Calls MITK Workbench for manual segmentations (horrible work around)
%   - Calls R2016b to run cardsync_interslice (fmincon issue)
%   ---- Hopefully this hasn't affected recon quality at all.
%   - Had to find older functions in R2015 for some scripts
%
%   FCMR_4D_RECON_BEASTIE01( ..., 'scanDate', scanDate, 'patID', patID )
% 
%   FCMR_4D_RECON_BEASTIE01( ..., 'reconDir', reconDir )
%
%

% tar (t.roberts@kcl.ac.uk)
% jfpva (joshua.vanamerom@kcl.ac.uk)


%% Optional Input Argument Default Values

default.scanDate    		= '';  % e.g.: '2011_07_22'
default.patID			    = '';  % e.g.: 'SA_385591'
default.reconDir            = '';  % must end .../<str>%i%i%i - i.e.: .../fcmr202
default.targetVol           = 1;   % TODO: implement option to specify target volume

%% Parse Input

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

add_param_fn(   p, 'scanDate', default.scanDate, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'patID', default.patID, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'reconDir', default.reconDir, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

parse( p, varargin{:} );

scanDate		    = p.Results.scanDate;
patID  				= p.Results.patID;
reconDir       		= p.Results.reconDir;



%% Directories / Dependencies

% Matlab code
% TODO: make these paths relative to files in Git repo download
addpath( genpath( '/home/cmo19/MATLAB/MRecon-3.0.557' ) ); % nb: do not use 3.0.515 - issue with recon oversampling matrix size
addpath( genpath( '/home/cmo19/MATLAB/ktrecon-dev' ) );
addpath( genpath( '/home/cmo19/MATLAB/fetal_cmr_4d-master' ) );
addpath( genpath( '/home/cmo19/MATLAB/TomR_Matlab_Code') );
addpath( genpath( '/home/cmo19/MATLAB/datafun') );

% Path to MITK Workbench
mitkWkbhPath = '/home/cmo19/MITK-2016.11.0-linux64/MitkWorkbench.sh';

% Path to Python3
pyPath = '/home/cmo19/miniconda/bin/python';

% Paths to reconstruction bash scripts
sh.scriptPath   = '/home/cmo19/MATLAB/fetal_cmr_4d-master/4drecon';
sh.RefVol       = fullfile( sh.scriptPath, 'recon_ref_vol.bash' );
sh.DcVol        = fullfile( sh.scriptPath, 'recon_dc_vol.bash' );
sh.SliceCineVol = fullfile( sh.scriptPath, 'recon_slice_cine.bash' );
sh.CineVol      = fullfile( sh.scriptPath, 'recon_cine_vol.bash' );
sh.VelVol       = fullfile( sh.scriptPath, 'recon_vel_vol.bash' );


%% GUI - ask user what to perform
choiceList = {'Part 0) k-t Reconstruction', ...
              'Part 1.1) Partial 3D Reconstruction     (REQUIRES: Part 0 // heart masks)', ...
              'Part 1.2) 4D Magnitude Reconstruction     (REQUIRES: Part 0 & 1.1 // heart and chest masks)', ...
              'Part 2) 4D Flow Reconstruction     (REQUIRES: Part 0 & PART 1 // uterus masks)', ...
              'Part 1.2) and Part 2)', ...
              'FULL Reconstruction     (NOT tested / NOT recommended as lengthy: will prompt for interaction)', ...
              'Run MITK     (for drawing masks)', ...
             };

[dlgChoice,~] = listdlg( ...
    'PromptString', 'What would you like to run?:', ...
    'SelectionMode', 'Single', ...
    'Name', 'FCMR 4D Reconstruction Options', ...
    'ListSize', [500,150], ...
    'ListString', choiceList );

if     dlgChoice == 1
    reconChoice = 'recon_kt';
    fprintf('Part 0) Running k-t reconstruction ...\n');    
elseif dlgChoice == 2
    reconChoice = 'recon_ref_vol';
    fprintf('Part 1.1) Running 3D reference ref_vol reconstrution ...\n');    
elseif dlgChoice == 3
    reconChoice = 'recon_cine_vol';
    fprintf('Part 1.2) Running 4D magnitude cine_vol reconstruction ...\n');
elseif dlgChoice == 4    
    reconChoice = 'recon_vel_vol';
    fprintf('Part 2) Running 4D flow vel_vol reconstruction ...\n');     
elseif dlgChoice == 5    
    reconChoice = 'recon_4d_vols';
    fprintf('Part 1.2 and Part 2) Running 4D magnitude and 4D flow reconstructions ...\n'); 
elseif dlgChoice == 6    
    reconChoice = 'recon_full';
elseif dlgChoice == 7    
    run_mitk;
    return;
elseif isempty(dlgChoice)
    fprintf( 'Exiting ...\n' );
	return;
end


%% Directories

% Ingenia raw folder
ingeniaRawDirPath 	     = '/pnraw01/raw-ingenia';
ingeniaArchiveRawDirPath = '/isi01/archive-rawdata/archive-ingenia';

% reconDir
studyDir = '/scratch/cmo19/Data'; %TODO allow different user

if isempty(reconDir)
    reconDir = uigetdir( studyDir, 'Select fcmr folder:' ); % e.g.: reconDir = '/scratch/cmo19/Data/fcmr202';
    
    % Cancel button
    if reconDir == 0
        fprintf( 'Exiting ...\n' );        
        return;
    end
end

fcmrNo = str2double( reconDir( end-2:end ) );

% Create directory stucture
cd(reconDir);

warning('off');
mkdir ktrecon
mkdir data
mkdir cardsync
mkdir mask
warning('on');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 0 --- Perform kt SENSE reconstruction on acquired MRI data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp( reconChoice, 'recon_kt' ) || strcmp( reconChoice, 'recon_full' )

	%% 0.1) Locate stacks on .raw server
	cd(reconDir);

	% Select patient .raw folder if not specified
    if ~isempty( scanDate ) && ~isempty( patID )
		
		try
			cd( fullfile( ingeniaRawDirPath, scanDate, patID ) );     
		catch
			cd( fullfile( ingeniaArchiveRawDirPath, scanDate, patID ) );
        end
        
        rawDataFileNames = dir('*kt8*.raw');
        rawDataFileNames = {rawDataFileNames.name}; % array for next part
        rawDataDirPath = pwd;
		
	elseif isempty( scanDate ) || isempty( patID )
        
        cd( ingeniaArchiveRawDirPath );
        
        [rawDataFileNames,rawDataDirPath] = uigetfile('*kt8*.raw', ...
            'Select kt8 bffe files:', ...
            'Multiselect','on');
        
        % Cancel button
        if isnumeric(rawDataFileNames)
            fprintf( 'Exiting ...\n' );        
            return;
        end
        
	else    
		error( 'Issue with specified scan data / patient ID.' );
    end

    
    % Auto-detect seriesNos
    cd( rawDataDirPath );
    
    for ii = 1:numel(rawDataFileNames)

        % parse underscores
        usLoc = strfind( rawDataFileNames{ii}, '_' );
        us3 = usLoc(3);
        us4 = usLoc(4);
        seriesNos(ii) = str2double( rawDataFileNames{ii}( us3+1:us4-1 ) );

    end

	% Stacks confirmation
	prompt = sprintf(['Which stacks do you want to reconstruct?\n' ...
		'Auto-detected: %s \n' ... 
		'Enter space-separated numbers: \n'], ...
        num2str(seriesNos) );
	dlgtitle = 'Input';
    definput = {num2str(seriesNos)};
	seriesNos = inputdlg( prompt, dlgtitle, [1 60], definput );
	seriesNos = sscanf(char(seriesNos),'%i')';

    if isempty(seriesNos)
        fprintf( 'Exiting ...\n' );
        return;
    end
    
	% %% Debug - check pnraw data found
	% for seriesNo = seriesNos
	% 	[ rawDataFilePath, coilSurveyFilePath, senseRefFilePath ] = id_pnraw_data( rawDataDirPath, seriesNo );
	% end


	%% 0.2) Perform kt reconstruction
	cd( reconDir );

	maskDirPath        = [];
	cusTrnDirPath      = [];
	makeHarmonicFilter = false;
	isSelfCaliPreProc  = [];
	patchVersion 	   = 'PIH1';
	isGeoCorrn         = true;
	outputDirPath      = fullfile( reconDir, 'ktrecon');

	recon_exam( fcmrNo, seriesNos, rawDataDirPath, patchVersion, isGeoCorrn, maskDirPath, cusTrnDirPath, makeHarmonicFilter, isSelfCaliPreProc, outputDirPath );

	fprintf('\n____________ Part 0) k-t RECONSTRUCTION COMPLETE ____________\n\n');


	%% copy dc / rlt files to data folder
	cd( reconDir );

	copyfile( 'ktrecon/s*dc_ab.nii.gz',  'data' );
	copyfile( 'ktrecon/s*rlt_ab.nii.gz', 'data' );
    
    if strcmp( reconChoice, 'recon_kt' )       
        fprintf('\nNext step: draw heart and uterus masks for each stack.\n');
        return;
    end

end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1 --- 4D Magnitude CINE volume reconstruction:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp( reconChoice, 'recon_ref_vol' )  || ...
   strcmp( reconChoice, 'recon_full' )

	%% 1.1) Manually draw fetal heart masks for each stack in MITK
	%       - use sXX_dc_ab.nii.gz
	%       - save as sXX_mask_heart.nii.gz in /mask


	%% 1.2) Run preproc.m
	cd(reconDir);    
	fprintf('Running preproc.m ...\n');
	S = preproc( reconDir );
	save( fullfile( reconDir, 'data', 'results.mat' ), 'S', '-v7.3' );
	fprintf('Completed preproc.m ...\n');


	%% 1.3) Run recon_ref_vol.bash
    fprintf('Reconstructing ref_vol ...\n');    
	cmdStr = sprintf('''%s'' ''%s'' ''%s''', sh.RefVol, reconDir, 'ref_vol');
	system(cmdStr);
    fprintf('3D reference volume (ref_vol.nii.gz) reconstructed ...\n');
    
    fprintf('\n____________ Part 1.1) 3D ref_vol RECONSTRUCTION COMPLETE ____________\n\n');
        
    % Exit if only reconstructing ref_vol
    if strcmp( reconChoice, 'recon_ref_vol' )        
        fprintf('\nNext step: draw chest mask using ref_vol.nii.gz.\n');
        return;
    end
   
end


if strcmp( reconChoice, 'recon_cine_vol' )  || ...
   strcmp( reconChoice, 'recon_4d_vols' )   || ...
   strcmp( reconChoice, 'recon_full' )
    
	%% 1.4) Draw fetal chest ROI using ref_vol.nii.gz
	% - save as mask_chest.nii.gz in /mask

    if strcmp( reconChoice, 'recon_full' )
    
        refVolNiiFilePath = fullfile( reconDir, 'ref_vol', 'ref_vol.nii.gz' );
        cmdStr = sprintf('''%s'' ''%s''', mitkWkbhPath, refVolNiiFilePath );
        [~,~] = system(cmdStr);
    
    end

	%% 1.5) Run recon_dc_vol.bash
    fprintf('Reconstructing dc_vol ...\n');    
	cmdStr = sprintf('''%s'' ''%s'' ''%s''', sh.DcVol, reconDir, 'dc_vol');
	system(cmdStr);
    fprintf('3D dc volume (dc_vol.nii.gz) reconstructed ...\n');
    

	%% 1.6) Run cardsync_intraslice.m
	cd(reconDir);
    fprintf('Running cardsync_intraslice.m ...\n');  
	dataDir     = fullfile( reconDir, 'data' );
	cardsyncDir = fullfile( reconDir, 'cardsync' );
	M = matfile( fullfile( dataDir, 'results.mat' ) );
	S = cardsync_intraslice( M.S, 'resultsDir', cardsyncDir, 'verbose', true );
	fprintf('Completed cardsync_intraslice.m ...\n');


	%% 1.7) Run recon_slice_cine.bash
	cd(reconDir);
    fprintf('Reconstructing slice cine volumes (slice_cine_vol) ...\n'); 	    
    cmdStr = sprintf('''%s'' ''%s''', sh.SliceCineVol, reconDir);
	system(cmdStr);
    fprintf('Slice cine volumes reconstructed ...\n'); 


	%% 1.8) Run cardsync_interslice.m
	% - this outputs /cardsync/results_interslice_cardsync.mat
	%               /cardsync/cardphases_interslice_cardsync.txt
	%               /cardsync/log_cardsync_interslice.txt

	cd(reconDir);
    
    fprintf('Running cardsync_interslice.m ...\n'); 

% 	% setup
% 	reconDir    = pwd;
% 	dataDir     = fullfile( reconDir, 'data' );
% 	cardsyncDir = fullfile( reconDir, 'cardsync' );
% 	cineDir     = fullfile( reconDir, 'slice_cine_vol' );    
% 	M = matfile( fullfile( cardsyncDir, 'results_cardsync_intraslice.mat' ) );
% 
% 	% target slice
% 	tgtLoc = NaN;
% 	tgtLocFile = fullfile( dataDir, 'tgt_slice_no.txt' );
% 	if exist( tgtLocFile , 'file' )
% 	  fid = fopen( tgtLocFile, 'r' );
% 	  tgtLoc = fscanf( fid, '%f' );
% 	  fclose( fid );
% 	end
% 
% 	% excluded slices
% 	excludeSlice = [];
% 	excludeSliceFile = fullfile( dataDir, 'force_exclude_slice.txt' );
% 	if exist( excludeSliceFile , 'file' )
% 	  fid = fopen( excludeSliceFile, 'r' );
% 	  excludeSlice = fscanf( fid, '%f' ) + 1;  % NOTE: slice locations in input file are zero-indexed
% 	  fclose( fid );
% 	end

	% slice-slice cardiac synchronisation
    cardsync_interslice_R2016b; % Another work around :(
    fprintf('Completed cardsync_interslice.m ...\n');
    
% 	S = cardsync_interslice( M.S, 'recondir', cineDir, 'resultsdir', cardsyncDir, 'tgtloc', tgtLoc, 'excludeloc', excludeSlice );


	%% 1.9) Run recon_cine_vol.bash
    fprintf('Reconstructing 4D magnitude cine volume (cine_vol) ...\n'); 
    cmdStr = sprintf('''%s'' ''%s'' ''%s''', sh.CineVol, reconDir, 'cine_vol');
	system(cmdStr);
    fprintf('4D magnitude cine volume reconstructed ...\n'); 


	%% 1.10) Summarise in MATLAB
	% (non-essential)
	% TODO: fix - no longer works with move to MIRTK

	% S = summarise_recon( [reconDir '/cine_vol'], [reconDir '/cardsync'], 'verbose', true );
	% I = plot_info( [reconDir '/cine_vol/info.tsv'] );
    
    
    %% 1.11) Create FCMR 4D magnitude dicoms
    % Python script: .../lib/dicom/fcmr_4d_make_dicom.py
    % Uses Miniconda install on beastie01 - TODO: make more generalisable
    % --recon_vel 0 option reconstructs 4D magnitude only.
    
    cd(reconDir);
    fprintf('Creating FCMR 4D magnitude dicoms ...\n');   
    
    exportPathCmd = 'export PATH="/usr/bin:/bin:/usr/sbin:/sbin:/home/cmo19/miniconda/bin";'; %nb: Linux PATH not known by Matlab
    pythonCmd     = ['/home/cmo19/miniconda/bin/python3 /home/cmo19/MATLAB/fetal_cmr_4d-master/lib/dicom/fcmr_4d_make_dicom.py -r ' reconDir '/ -f ' num2str(fcmrNo) ' --recon_vel 0'];
    cmdStr        = strcat(exportPathCmd, pythonCmd);
	system(cmdStr);
    
    fprintf('Completed dicom creation ...\n');
    
    
    fprintf('\n____________ Part 1.2) 4D magnitude volume (cine_vol) RECONSTRUCTION COMPLETE ____________\n\n');

end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2 --- 4D Velocity CINE volume reconstruction:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 2 // Full Reconstruction
if strcmp( reconChoice, 'recon_vel_vol' ) || ...
   strcmp( reconChoice, 'recon_4d_vols' ) || ...
   strcmp( reconChoice, 'recon_full' )

	%% 2.1) Manually draw uterus masks for each stack in MITK
	%       - use sXX_dc_ab.nii.gz
	%       - save as sXX_mask_uterus.nii.gz in /mask


	%% 2.2) Run fcmr_4dflow_preprocessing.m
    cd(reconDir);
    fprintf('Running fcmr_4dflow_preprocessing.m ...\n');
    fcmr_4dflow_preprocessing( reconDir );
	fprintf('Completed fcmr_4dflow_preprocessing.m ...\n');


	%% 2.3) Run fcmr_4dflow_get_first_moments.m

	% - NOTE: this script creates gradient first moment .txt files using 
	% data created by reconFrame. Alternatively, the two .txt files required 
	% for velocity volume reconstruction can be calculated manually if you know
	% the parameters of your flow encoding gradients

	cd(reconDir);
    fprintf('Running fcmr_4dflow_get_first_moments.m ...\n');
	fcmr_4dflow_get_first_moments( reconDir );
	fprintf('Completed fcmr_4dflow_get_first_moments.m ...\n');


	%% 2.4) Run recon_vel_vol.bash
    fprintf('Reconstructing 4D velocity cine volume (vel_vol) ...\n'); 
	cmdStr = sprintf('''%s'' ''%s'' ''%s''', sh.VelVol, reconDir, 'vel_vol');
	system(cmdStr);
    fprintf('4D velocity cine volume reconstructed ...\n'); 

    
	%% 2.5) Run fcmr_4dflow_postprocessing.m

	% - TODO: update script with MATLAB niftiread/niftiwrite to avoid having to
	% use reslice_nii.m

	cd(reconDir);
    fprintf('Running fcmr_4dflow_postprocessing.m ...\n');
    
    bloodpoolMask = 'mask_blood_pool'; 
	velMasks      = bloodpoolMask;
	fcmr_4dflow_postprocessing( reconDir, ...
		'useVelDriftCorr', true, ...
		'fileExt', 'polyCorr', ...
		'bloodpoolMask', bloodpoolMask, ...
		'velMasks', bloodpoolMask );
	
    fprintf('Completed fcmr_4dflow_postprocessing ...\n');
    
    
    %% 2.6) Create FCMR 4D flow dicoms
    % Python script: .../lib/dicom/fcmr_4d_make_dicom.py
    % Uses Miniconda install on beastie01 - TODO: make more generalisable
    
    cd(reconDir);
    fprintf('Creating FCMR 4D flow dicoms ...\n');
    
    % Remove existing dcm_4d folder created in Part 1
    rmdir( fullfile( reconDir, 'dcm_4d' ), 's' );
    
    exportPathCmd = 'export PATH="/usr/bin:/bin:/usr/sbin:/sbin:/home/cmo19/miniconda/bin";'; %nb: Linux PATH not known by Matlab
    pythonCmd     = ['/home/cmo19/miniconda/bin/python3 /home/cmo19/MATLAB/fetal_cmr_4d-master/lib/dicom/fcmr_4d_make_dicom.py -r ' reconDir '/ -f ' num2str(fcmrNo)];
    cmdStr        = strcat(exportPathCmd, pythonCmd);
	system(cmdStr);
    
    fprintf('Completed dicom creation ...\n');
    
    
    fprintf('\n____________ Part 2) 4D flow volume (vel_vol) RECONSTRUCTION COMPLETE ____________\n\n');

end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3 --- Whole-heart 4D visualisation in Paraview:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: automate generation of 4D flow Paraview file

%% 3.1) Run fcmr_4dflow_make_vector_vol.py
%
% In python:
% - edit fcmr_4dflow_make_vector_vol.py with correct path information
% Run: fcmr_4dflow_make_vector_vol.py


%% 3.2) Load state file in Paraview.
% /vel_vol_4d/paraview_polyCorr contains .pvsm file
% In Paraview:
% - File > Load State > fcmr***_paraview.pvsm



end % FCMR_4D_RECON_BEASTIE01(...)


function run_mitk( )

%% RUN_MITK
%
% Quick function to call MITK Workbench using system command
%
% Note: MITK can't be called from within Matlab 2015a on beastie01
% Hence, this _horrible_ work around calls Matlab 2016b within the
% run_mitk.bash script.
%
%

% Path to MITK Workbench bash script
mitkWkbhScriptPath = '/home/cmo19/MATLAB/run_mitk.bash';

% Load MITK Workbench
fprintf('\nLoading MITK Workbench ...\n');
cmdStr = sprintf('''%s''', mitkWkbhScriptPath );
[~,~] = system(cmdStr);
fprintf('Exiting MITK Workbench.\n');


% end run_mitk(...)
end


function cardsync_interslice_R2016b( )

%% cardsync_interslice_R2016b
%
% Run cardsync_interslice.m in R2016b (horrible workaround)
%
% Note: can't run easily in R2015a as fmincon changed between versions
%
%

% Path to MITK Workbench bash script
cardsyncInterScriptPath = '/home/cmo19/MATLAB/cardsync_interslice_R2016b.bash';

% Load MITK Workbench
cmdStr = sprintf('''%s''', cardsyncInterScriptPath );
[~,~] = system(cmdStr);


% end run_mitk(...)
end