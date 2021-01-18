%% k-t Training Data Harmonics


%% Load Training Data
studyDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper';
cd(studyDir);

ktreconDir = 'ktrecon';
maskDir    = 'mask';

dimT = 3;
xt2xf = @( xt ) fftshift( fft( xt, [], dimT ), dimT );

caseNos = [189 191 194 197 201 202 206 214 254];

for iC = 1:numel(caseNos)
    
    fprintf( 'Loading fcmr%d ... \n', caseNos(iC)  );
    
    reconDir = [ 'fcmr' num2str(caseNos(iC)) ];
    
    % Load Heart Masks
    cd( fullfile( studyDir, reconDir, maskDir ) );
    maskFiles = dir( 's*_mask_heart.nii.gz' );
    numStacks = numel(maskFiles);
    
    for iM = 1:numStacks
        
        mask = niftiread( maskFiles(iM).name );
        
        S(iM).mask_heart = squeeze( mask );
        
    end
    
    % Load RLT / TRN Data (saved as xtPri in Josh scripts)
    % nb: just single dynamic for RLT to locate heart
    cd( fullfile( studyDir, reconDir, ktreconDir ) );
    priFiles  = dir('s*_pri_recon.mat');
    rltFiles  = dir('s*_rlt_recon.mat');
    
    for iP = 1:numStacks
        
        load( priFiles(iP).name );
        load( rltFiles(iP).name );
        
        S(iP).xtPri    = squeeze( xtPri );
        S(iP).xtRcn    = squeeze( xtRcn(:,:,:,:,1,:,:,:) );
        S(iP).seriesNo = str2double( erase( erase( priFiles(iP).name,'_pri_recon.mat' ), 's') );
        S(iP).xfPri    = xt2xf( squeeze( xtPri ) );
        
    end

    % Save Stacks to Cases Structure
    C(iC).fcmrNum = caseNos(iC);
    C(iC).S = S;
    
    clear S xtPri xtRcn mask maskFiles priFiles rltFiles
    
end

fprintf( 'Collated All Training Datasets. \n\n'  );


%% Identify x-f Line From RLT Data

caseNo   = 3;
seriesNo = 5;

xtRcn = C(caseNo).S(seriesNo).xtRcn;
imtar( abs(xtRcn) );
set(gcf,'position',[20 40 1880 1070]);


%% View x-f of TRN Data

xfLine   = 60;
sliceNum = 9;

xtPri = C(caseNo).S(seriesNo).xtPri;
fetal_cmr_view_xtxf( xtPri, sliceNum, xfLine );


%% Examine Harmonics in Mask Heart Region

mask = C(caseNo).S(seriesNo).mask_heart;

% Find Edge Coordinates of Mask
[mRows, mCols , ~] = find(mask (:,:,sliceNum) );
mRows = unique( mRows ); mCols = unique( mCols );
xMin = min(mCols); xMax = max(mCols);
yMin = min(mRows); yMax = max(mRows);

xfPri = C(caseNo).S(seriesNo).xfPri;
xfPri(:,:,49,:) = 0;

xfHrmRegion = squeeze( xfPri( yMin:yMax,xfLine,:,sliceNum ) );

figure;
subplot(1,2,1);
imagesc( abs(xfHrmRegion), [0, 5 * mean(abs(xfHrmRegion(:)))] ); colormap('parula');
axis('image'); title('x-f Harmonic Region');
subplot(1,2,2);
plot( mean( abs(xfHrmRegion), 1  ) );
title('Mean x-f');
set(gcf,'Position',[900 670 1000 420]);



