

%% ISMRM 2020 Abstract Figures


%% FIGURE 1 --- x-f SPACES COMPARISON

% Run this: E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr194_hrh\fcmr194_hrmflt_hrmsOnly_reg_testing.m
% Then pick out necessary x-f spaces

% fcmr194, s18, slice9, xfLine = 57;

% get standard recon xfPri / xfMask
load('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr194_hrh\ktrecon_hrmflt_hrmsOnly_s18_alpha_testing\fcmr194_std_uniformReg_s18_sl9.mat','xfPri','xfMask');
xfPriUniformReg  = xfPri;
% xfMaskUniformReg = xfMask;
% dimX = 1; fn_rmv_os = @( x ) x((round(size(x,dimX)/4)+1):(size(x,dimX)-round(size(x,dimX)/4)),:,:,:);
% xfMaskUniformReg = fn_rmv_os(xfMaskUniformReg);

xfLine = 57;
xCrop = 20:170;

% get xfMask
load('hrmflt_alpha10.mat','xfMask','xfPri');
xfHarmonicMask = squeeze(xfMask(xCrop,xfLine,:)./max(xfMask(:)));

% xf
stdTrn = squeeze(abs(xfPriJ(xCrop,xfLine,:))); stdTrn(:,49)=0;
stdPri = squeeze(abs(xfPriUniformReg(xCrop,xfLine,:))); stdPri(:,49)=0;
stdPriReg = stdPri./max(stdPri(:)); stdPriReg(stdPriReg<0.0267) = 0;
% stdPriReg = stdPri./max(stdPri(:)); stdPriReg(stdPriReg<0.03) = 0;
stdRcn = squeeze(abs(xfRcnJ(xCrop,xfLine,:))); stdRcn(:,49)=0;
hrhTrn = squeeze(abs(xfRcnL(xCrop,xfLine,:))); hrhTrn(:,49)=0;
hrhPri = squeeze(abs(xfPri(xCrop,xfLine,:)));  hrhPri(:,49)=0;
hrhRcn = squeeze(abs(XFRCN(xCrop,xfLine,:,7))); hrhRcn(:,49)=0;

% maxAlpha100 = max(max(abs(XFRCN(xCrop,xfLine,:,7))));
       
%stdPri ./ max(stdPri(:)), ...       
imtar([stdTrn ./ max(stdTrn(:)), ...
       0.1.*ones(size(xfHarmonicMask)), ...
       stdPriReg, ...
       stdRcn ./ max(stdRcn(:)); ...
       hrhTrn ./ max(hrhTrn(:)), ...
       xfHarmonicMask, ...
       hrhPri ./ max(hrhRcn(:)), ...
       hrhRcn ./ max(hrhRcn(:))],0,.25); colormap('parula'); colorbar('off');
axis('off');

% xt 
yCrop = 25:80;
% yCrop = 1:144;

stdXtPri = xtPriJ(xCrop,yCrop,1);
stdXtRcn = xtRcnJ(xCrop,yCrop,1);
hrhXtPri = xtRcnL(xCrop,yCrop,1);
hrhXtRcn = XTRCN(xCrop,yCrop,1,7);

imtar([stdXtPri ./ max(stdXtPri(:)), ...
       stdXtRcn ./ max(stdXtRcn(:)); ...
       hrhXtPri ./ max(hrhXtPri(:)), ...
       hrhXtRcn ./ max(hrhXtRcn(:))]); colormap('gray'); colorbar('off');
axis('off');








%% FIGURE 2 --- CINE COMPARISON

% fcmr191, s20
fcmrNum = 191;
cd(['E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr' num2str(fcmrNum) '\data']);
s20.std = niftiread('s20_rlt_ab.nii.gz');

cd(['E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr' num2str(fcmrNum) '_hrh_fullRecon\data']);
s20.hrh = niftiread('s20_rlt_ab.nii.gz');

sl = 5;
% implay_RR([s20.std(:,:,sl,:) s20.hrh(:,:,sl,:) s20.std(:,:,sl,:)-s20.hrh(:,:,sl,:)]);

rowCrop = 40:160;
colCrop = 10:100;
% imtar(s20.std(rowCrop,colCrop,sl,1));

%%% CINE COMPARISON MOVIE

stdCrop = squeeze(s20.std(rowCrop,colCrop,sl,:)); stdCrop = stdCrop./max(stdCrop(:));
hrhCrop = squeeze(s20.hrh(rowCrop,colCrop,sl,:)); hrhCrop = hrhCrop./max(hrhCrop(:));
diffCrop = 10.*(stdCrop-hrhCrop);
    
     

MOVIE = [stdCrop hrhCrop diffCrop];
% implay_RR(MOVIE);

saveDir = 'C:\Users\tr17\Documents\Conferences\ISMRM_2020\4DFlow_kt_Hierarchical';
gifName = ['fcmr' num2str(fcmrNum) 'cine_std_v_hrh.gif'];
saveFilePath = fullfile(saveDir, gifName);

fps = 20;
% MOVIE=MOVIE./max(MOVIE(:));
figure; %set(gcf, 'Position', [1, 31, 1920, 1093]);
% str = 'a';
% dim = [.1 .1 .1 .1];
for ii = 1:size(MOVIE,3)
    [movieFrame,cm]=gray2ind(MOVIE(:,:,ii),256);
    imagesc(movieFrame); colormap(cm); axis('image'); axis('off'); 
    title( sprintf('with calibration data') );
%     title('with calibration data');
    drawnow;
%     annotation('textbox',dim,'String',str,'FitBoxToText','on');    
    if ii == 1
        imwrite(movieFrame, cm, saveFilePath, 'gif', 'Delay', 1/fps, 'Loop', inf );
    else
        imwrite(movieFrame, cm, saveFilePath, 'gif', 'Delay', 1/fps, 'WriteMode','append');
    end
end
close;



%% FIGURES 3 & 4 --- 4D Flow Figures
% -fcmr194
% -fcmr202
% --- std recon vs. hrh recon

% Created using:
% C:\Users\tr17\Documents\CODE_Projects\fetal_cmr_4dflow\wip\ismrm2020_abstract\ISMRM2020_Abstract_Figures_3_4.m


%% FIGURE 5 --- Mean Flow bar chart
% Created using:
% C:\Users\tr17\Documents\CODE_Projects\fetal_cmr_4dflow\wip\ismrm2020_abstract\fcmr_4dflow_velocity_analysis_standard_v_HRH.m

% Then ported into Graphpad

% Original values in Excel file:
% ismrm2020_abstract_HRH_Harmonic_notes_analysis.xlsx
%
% Copied manually into here for paired t-test:

load('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\sc_slw\hrh_meanflow_comparison_ttest');

rowOrder = {'AAo','DA','DAo','PA','SVC'}; % from excel file
for ii = 1:5; [h(ii),p(ii)]=ttest(std(ii,:),hrh(ii,:)); end % paired t-test


%% Temporal Dynamics Image
% NOT USED

% fcmr191, s20
fcmrNum = 191;
nF = 96;
cd(['E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr' num2str(fcmrNum) '\data']);
s20.std = niftiread('s20_rlt_ab.nii.gz');

cd(['E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr' num2str(fcmrNum) '_hrh_fullRecon\data']);
s20.hrh = niftiread('s20_rlt_ab.nii.gz');

sl = 5;
implay_RR([s20.std(:,:,sl,:) s20.hrh(:,:,sl,:) s20.std(:,:,sl,:)-s20.hrh(:,:,sl,:)]);


I = s20.std;
imtar(I(:,:,sl,1));
set(gcf, 'Position', [1, 31, 1920, 1093]);
[cX,cY,iVal,iX,iY] = improfile;
close;

% Standard k-t
for ii = 1:nF
    [~,~,iVal(:,ii)] = improfile( I(:,:,sl,ii),iX,iY );
end
s20.std_iVal = iVal; clear iVal;


% HRH k-t
I = s20.hrh;
for ii = 1:nF
    [~,~,iVal(:,ii)] = improfile( I(:,:,sl,ii),iX,iY );
end
s20.hrh_iVal = iVal; clear iVal;

imtar([s20.std_iVal; s20.hrh_iVal]);


%% Registered Volume Subtraction
% NOT USED --- too many subtle differences between datasets, therefore
% difficult to make meaningful comparison easily

fcmrNum = 202;
cd(['E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr' num2str(fcmrNum) '_hrh_fullRecon\4dvol_reg']);

cine_vol.std = niftiread('cine_vol_orig.nii.gz');
mask = ~(cine_vol.std==-1);
mask = imerode(mask, strel('sphere',2));

cine_vol.std = cine_vol.std .* mask;
cine_vol.hrh = mask .* niftiread('cine_vol_hhh_reg.nii.gz');
% cine_vol.diff = cine_vol.std - cine_vol.hrh;

% montage_RR(cine_vol.std(:,:,:,1));
% montage_RR(cine_vol.hrh(:,:,:,1));
% montage_RR(cine_vol.diff(:,:,:,1));

maxVal = max(cine_vol.std(:));
meanOrig = mean(nonzeros(cine_vol.std(:)));
meanHrh = mean(nonzeros(cine_vol.hrh(:)));

sl = 32;
% imtar(100.*[cine_vol.std(:,:,sl,1)./maxVal,...
%             cine_vol.hrh(:,:,sl,1)./maxVal,...
%             cine_vol.diff(:,:,sl,1)./maxVal ],0,100);

% montage_RR(mean(cine_vol.std,4));
% montage_RR(mean(cine_vol.diff,4));

cine_vol.std = cine_vol.std ./ meanOrig;
cine_vol.hrh = cine_vol.hrh ./ meanHrh;
cine_vol.diff = cine_vol.std - cine_vol.hrh;

imtar([cine_vol.std(:,:,sl,1),...
            cine_vol.hrh(:,:,sl,1),...
            cine_vol.diff(:,:,sl,1)]);













