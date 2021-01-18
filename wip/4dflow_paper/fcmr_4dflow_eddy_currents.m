
%% Eddy_current investigation
% - Investigating effect of eddy currents for 4D flow paper reviewers
%
% - ran FSL eddy_correct on Jana's computer using:
% - eddy_correct 201805241031_PHILIPSJ2LF9B9_301_WIPfromfileSENSE.nii.gz 201805241031_PHILIPSJ2LF9B9_301_WIPfromfileSENSE_ec.nii.gz 0

eddyDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Eddy_Current_Testing\data_diffusion\eddy_analysis';
cd(eddyDir);


%% load data
dwi    = load_untouch_nii('201805241031_PHILIPSJ2LF9B9_301_WIPfromfileSENSE.nii.gz');
dwi_ec = load_untouch_nii('201805241031_PHILIPSJ2LF9B9_301_WIPfromfileSENSE_ec.nii.gz');

dwi.img      = im2double( mat2gray( single(dwi.img) ) );
dwi_ec.img   = im2double( mat2gray( single(dwi_ec.img) ) );
dwi_diff.img = dwi.img - dwi_ec.img;

sl = 32;
dw = 1;

imtar([ dwi.img(:,:,sl,dw),    ...
        dwi_ec.img(:,:,sl,dw), ...
        dwi_diff.img(:,:,sl,dw) ]);
    
imtar([ dwi.img(:,:,sl,1), dwi.img(:,:,sl,2), dwi.img(:,:,sl,3); ...
        dwi_ec.img(:,:,sl,1), dwi_ec.img(:,:,sl,2), dwi_ec.img(:,:,sl,3) ], 0, 0.25);
    
imtar([ dwi_diff.img(:,:,sl,1), ...
        dwi_diff.img(:,:,sl,2), ...
        dwi_diff.img(:,:,sl,3)]);

%%
b0    = im2double( mat2gray( dwi.img(:,:,:,1) ));
b500  = im2double( mat2gray( dwi.img(:,:,:,2) ));
b1000 = im2double( mat2gray( dwi.img(:,:,:,3) ));

b0_ec    = im2double( mat2gray( dwi_ec.img(:,:,:,1) ));
b500_ec  = im2double( mat2gray( dwi_ec.img(:,:,:,2) ));
b1000_ec = im2double( mat2gray( dwi_ec.img(:,:,:,3) ));

imtar([ b0(:,:,sl),      ...
        b500(:,:,sl),    ...
        b1000(:,:,sl);   ...
        b0_ec(:,:,sl),   ...
        b500_ec(:,:,sl), ...
        b1000_ec(:,:,sl)]);

imtar([ b0(:,:,sl) - b1000(:,:,sl),      ...
        b0_ec(:,:,sl) - b1000_ec(:,:,sl)]);



