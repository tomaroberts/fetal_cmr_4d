%% get data / masks for fcmr segNet

sp1.path = strcat('I:\King', '''', 's College London\Perinatal CMR - Documents\General\Data\fcmr_4dflow_paper');
sp1.prefix = 'fcmr';
sp1.fcmrNums = [189, 191, 194, 197, 201, 202, 206, 213, 214, 230, 254, 255];
sp1.N = numel(sp1.fcmrNums);

sp2.path = strcat('I:\King', '''', 's College London\Perinatal CMR - Documents\General\Data\fcmr_4d_clinical');
sp2.prefix = 'c_fcmr';
sp2.fcmrNums = [240, 242, 246, 250, 252, 253, 256, 266, 267, 271, 273, 285, ...
                286, 287, 288, 289, 291, 293, 295, 302, 305, 308, 311, ...
                316, 319, 327, 329];
sp2.N = numel(sp2.fcmrNums);            


%% Make output dir
imDir  = 'C:\Users\tr17\Dropbox\fcmr_seg\images';
segDir = 'C:\Users\tr17\Dropbox\fcmr_seg\segs';

foldNums = sort([sp1.fcmrNums, sp2.fcmrNums]);

for ii = 1:numel(foldNums)
    cd(imDir);
    mkdir(num2str(foldNums(ii)));
    
    cd(segDir);
    mkdir(num2str(foldNums(ii)));
end


%% Collate 

% sp1
for iN = 1:sp1.N
        
    % images
    cd( fullfile( sp1.path, strcat( sp1.prefix, num2str(sp1.fcmrNums(iN))), 'data' ) );    
    dcFiles = dir('*_dc_ab.nii.gz');
    outDir  = fullfile( imDir, num2str(sp1.fcmrNums(iN)) );
    
    for iDC = 1:numel(dcFiles)
        copyfile( dcFiles(iDC).name, outDir );
    end
    
    % mask segmentations
    cd( fullfile( sp1.path, strcat( sp1.prefix, num2str(sp1.fcmrNums(iN))), 'mask' ) );    
    maskFiles = dir('*_mask_heart.nii.gz');
    outDir  = fullfile( segDir, num2str(sp1.fcmrNums(iN)) );
    
    for iM = 1:numel(maskFiles)
        copyfile( maskFiles(iM).name, outDir );
    end
   
    disp(['Copied number ... ' num2str(sp1.fcmrNums(iN)) ]);
    
end


% sp2
for iN = 1:sp2.N
        
    % images
    cd( fullfile( sp2.path, strcat( sp2.prefix, num2str(sp2.fcmrNums(iN))), 'data' ) );    
    dcFiles = dir('*_dc_ab.nii.gz');
    outDir  = fullfile( imDir, num2str(sp2.fcmrNums(iN)) );
    
    for iDC = 1:numel(dcFiles)
        copyfile( dcFiles(iDC).name, outDir );
    end
    
    % mask segmentations
    cd( fullfile( sp2.path, strcat( sp2.prefix, num2str(sp2.fcmrNums(iN))), 'mask' ) );    
    maskFiles = dir('*_mask_heart.nii.gz');
    outDir  = fullfile( segDir, num2str(sp2.fcmrNums(iN)) );
    
    for iM = 1:numel(maskFiles)
        copyfile( maskFiles(iM).name, outDir );
    end
   
    disp(['Copied number ... ' num2str(sp2.fcmrNums(iN)) ]);
    
end

