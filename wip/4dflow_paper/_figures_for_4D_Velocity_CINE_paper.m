%% Figures for paper

paperDir = 'C:\Users\tr17\Documents\Papers\4D_Velocity_CINE_Paper';


%% FIGURE: Synthetic phantom

cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\Synthetic_Flow_Phantom\bSSFP_6pipes');

% load('VEL_VOLUME_no_padding.mat');
phantom = load_nii('VEL_VOLUME_srow4_equal0.nii');
vel_vol = load_nii('vel_vol_5stacks/velocity-vector-final.nii.gz');

% % need to transform Vx/Vy/Vz into same space as vel_vol for difference
% % images:
% vv = vel_vol;
% vv.img = [];
% vv.img(:,:,:,1) = Vx; vv.img(:,:,:,2) = Vy; vv.img(:,:,:,3) = Vz;
% vv.hdr.dime.dim(4) = size(Vx,3);
% save_nii(vv,'VEL_VOLUME_no_padding.nii.gz');
% 
% % mirtk transform-image velocity-vector-final.nii.gz velocity-vector-final_transformed.nii.gz -target ../VEL_VOLUME_no_padding.nii.gz

x = 145:272;
y = 141:268;
z = 250:349;

diff.img = phantom.img(y+1,x+1,z+1)-sum(vel_vol.img(:,:,:,:),4);

implay_RR([phantom.img(y+1,x+1,z+1) sum(vel_vol.img(:,:,:,:),4) diff.img(:,:,:,:) ],'jet',[-1.1,1.1]);

sl = 30;
imtar([phantom.img(y+1,x,z(1)+sl-1) sum(vel_vol.img(:,:,sl,:),4) diff.img(:,:,sl)],-1,1); colormap('jet');       
set(gca,'xticklabel',[]); set(gca,'xtick',[]); set(gca,'yticklabel',[]); set(gca,'ytick',[]); set(gca,'visible','off');

% mean differences analysis:
phantom_crop = phantom.img(y+1,x+1,z+1);
% phantom_crop_slice = phantom_crop;
phantom_crop_slice = phantom_crop(:,10:115,30);
idxPipes = find(phantom_crop_slice);
idxBkgrd = find(phantom_crop_slice==0);
diffPipes = abs(mean(diff.img(idxPipes)) * 100)
dirrBkgrd = abs(mean(diff.img(idxBkgrd)) * 100)


%% FIGURE: Physical Flow Phantom

% QFlow:
cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\#Finalised_Data\flow_phantom6_qflow_Tom_recon');
Q = load_nii('3D_QFlow_ffd.nii.gz');
tra = load_nii('tra_MAG.nii.gz'); %magnitude image for figure

% PC-bSSFP (use PC-bSSFP transformed into QFlow grid):
cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\#Finalised_Data\flow_phantom6_bssfp_ffd_registered2qflow');
B = load_nii('vel_vol_bssfp2qflow_3D.nii.gz');

perm = [2,1,3,4];
Q.img = flip(flip(permute(Q.img,perm),1),2);
B.img = flip(flip(permute(B.img,perm),1),2);

% implay_RR([Q.img(:,:,:,1) Q.img(:,:,:,2) Q.img(:,:,:,3); ...
%            B.img(:,:,:,1) B.img(:,:,:,2) B.img(:,:,:,3)] ...
%            ,'jet',[-50,50]);
       
x = 110:150-1;
y = 125:165-1;
sl = 50;

% Q v B:
imtar([Q.img(x,y,sl,1) Q.img(x,y,sl,2) Q.img(x,y,sl,3); ...
       B.img(x,y,sl,1) B.img(x,y,sl,2) B.img(x,y,sl,3)],-50,50); colormap('jet');

%% % loop through individual images
% cd('C:\Users\tr17\Documents\Papers\4D_Velocity_CINE_Paper');
% imStr = {'Vx','Vy','Vz'};
% for ii = 1:3
%     imtar(Q.img(x,y,sl,ii),-50,50); colormap('jet'); set(gca,'xticklabel',[]); set(gca,'xtick',[]); set(gca,'yticklabel',[]); set(gca,'ytick',[]); set(gca,'visible','off'); colorbar('off');
% 
%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset; 
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];
%     
%     saveas(gcf,['FlowPhantom6_QFlow_' imStr{ii} '.png']);
%     
%     imtar(B.img(x,y,sl,ii),-50,50); colormap('jet'); set(gca,'xticklabel',[]); set(gca,'xtick',[]); set(gca,'yticklabel',[]); set(gca,'ytick',[]); set(gca,'visible','off'); colorbar('off');
% 
%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset; 
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];
%     
%     saveas(gcf,['FlowPhantom6_bSSFP_' imStr{ii} '.png']);
%     
%     close all
% end
% 
% % magnitude image:
% tra.img = flip(flip(permute(tra.img,perm),1),2);
% imtar(tra.img(:,:,sl)); set(gca,'xticklabel',[]); set(gca,'xtick',[]); set(gca,'yticklabel',[]); set(gca,'ytick',[]); set(gca,'visible','off'); colorbar('off');
% 
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% 
% cd('C:\Users\tr17\Documents\Papers\4D_Velocity_CINE_Paper');
% saveas(gcf,['FlowPhantom6_MagnitudeImage.png']);


%% get line profiles through phantom tubes:
%-line profile 1
imtar(Q.img(x,y,sl,2),-50,50); colormap('jet');
[Q.cx1,Q.cy1,Q.c1] = improfile;

[B.cx1(:,1),Q.cy1(:,1),Q.c1(:,1)] = improfile(Q.img(x,y,sl,1),Q.cx1(:,1),Q.cy1(:,1));
[B.cx1(:,1),B.cy1(:,1),B.c1(:,1)] = improfile(B.img(x,y,sl,1),Q.cx1(:,1),Q.cy1(:,1));
for ii = [2,3];
    [Q.cx1(:,ii),Q.cy1(:,ii),Q.c1(:,ii)] = improfile(Q.img(x,y,sl,ii),Q.cx1(:,1),Q.cy1(:,1));
    [B.cx1(:,ii),B.cy1(:,ii),B.c1(:,ii)] = improfile(B.img(x,y,sl,ii),Q.cx1(:,1),Q.cy1(:,1));
    close;
end

% plot
figure;
hold on;
for ii = 1:3
    plot(Q.c1(:,ii),'b');
    plot(B.c1(:,ii),'r');
end
plot(1:length(B.c1(:,ii)),zeros(length(B.c1(:,ii))),'k--');

%% plots in Graphpad:
%-copied B.c1 and Q.c1

% Peak velocity magnitudes through P1/P2:
%-P1:
q1 = [24.0540000000000,-24.6972000000000,-12.8788000000000]; norm(q1);
b1 = [19.5092500000000,-24.7532000000000,-13.2243600000000]; norm(b1);

%-P2:
q2 = [4.802399	-10.447600	-6.339600]; norm(q2);
b2 = [2.981817	-9.656078	-6.725181]; norm(b2);


%% 2D ROI version

% QFlow:
cd('C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\#Finalised_Data\flow_phantom6_qflow_Tom_recon');
Q = load_nii('3D_QFlow_ffd.nii.gz');
tra = load_nii('tra_MAG.nii.gz'); %magnitude image for figure

% % 2d ROIs
% Qmask1 = load_nii('Q_mask_P4.nii.gz'); %nb: drawn on 2nd component of 4th dimension
% Qmask2 = load_nii('Q_mask_P2.nii.gz');
% Qmask3 = load_nii('Q_mask_P3.nii.gz');

% 3d ROIs
Qmask1 = load_nii('Q_mask_2d_P1_v3.nii.gz');
Qmask2 = load_nii('Q_mask_2d_P2.nii.gz');
Qmask3 = load_nii('Q_mask_2d_P3.nii.gz');

% PC-bSSFP (use PC-bSSFP transformed into QFlow grid):
cd('C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\#Finalised_Data\flow_phantom6_bssfp_ffd_registered2qflow');
B = load_nii('vel_vol_bssfp2qflow_3D.nii.gz');

% % 2d ROIs
% Bmask1 = load_nii('B_mask_P4.nii.gz');
% Bmask2 = load_nii('B_mask_P2.nii.gz');
% Bmask3 = load_nii('B_mask_P3.nii.gz');

% 3d ROIs
Bmask1 = load_nii('B_mask_2d_P1_v3.nii.gz');
Bmask2 = load_nii('B_mask_2d_P2.nii.gz');
Bmask3 = load_nii('B_mask_2d_P3.nii.gz');

perm = [2,1,3,4];
Q.img = flip(flip(permute(Q.img,perm),1),2);
B.img = flip(flip(permute(B.img,perm),1),2);

Qmask1.img = double(flip(flip(permute(Qmask1.img,perm),1),2)); Qmask1.img(:,:,:,1:3) = cat(4,Qmask1.img(:,:,:,2), Qmask1.img(:,:,:,2), Qmask1.img(:,:,:,2));
Qmask2.img = double(flip(flip(permute(Qmask2.img,perm),1),2)); Qmask2.img(:,:,:,1:3) = cat(4,Qmask2.img(:,:,:,2), Qmask2.img(:,:,:,2), Qmask2.img(:,:,:,2));
Qmask3.img = double(flip(flip(permute(Qmask3.img,perm),1),2)); Qmask3.img(:,:,:,1:3) = cat(4,Qmask3.img(:,:,:,2), Qmask3.img(:,:,:,2), Qmask3.img(:,:,:,2));
Bmask1.img = double(flip(flip(permute(Bmask1.img,perm),1),2)); Bmask1.img(:,:,:,1:3) = cat(4,Bmask1.img(:,:,:,2), Bmask1.img(:,:,:,2), Bmask1.img(:,:,:,2));
Bmask2.img = double(flip(flip(permute(Bmask2.img,perm),1),2)); Bmask2.img(:,:,:,1:3) = cat(4,Bmask2.img(:,:,:,2), Bmask2.img(:,:,:,2), Bmask2.img(:,:,:,2));
Bmask3.img = double(flip(flip(permute(Bmask3.img,perm),1),2)); Bmask3.img(:,:,:,1:3) = cat(4,Bmask3.img(:,:,:,2), Bmask3.img(:,:,:,2), Bmask3.img(:,:,:,2));

% apply masks
for ii = 1:3
    Qmean1(:,ii) = nonzeros(Qmask1.img(:,:,:,ii) .* Q.img(:,:,:,ii));
    Qmean2(:,ii) = nonzeros(Qmask2.img(:,:,:,ii) .* Q.img(:,:,:,ii));
    Qmean3(:,ii) = nonzeros(Qmask3.img(:,:,:,ii) .* Q.img(:,:,:,ii));
    Bmean1(:,ii) = nonzeros(Bmask1.img(:,:,:,ii) .* B.img(:,:,:,ii));
    Bmean2(:,ii) = nonzeros(Bmask2.img(:,:,:,ii) .* B.img(:,:,:,ii));
    Bmean3(:,ii) = nonzeros(Bmask3.img(:,:,:,ii) .* B.img(:,:,:,ii));
end

% V = Vx+Vy+Vz
% Qmean1(:,4) = sum(Qmean1(:,1:3),2);
% Qmean2(:,4) = sum(Qmean2(:,1:3),2);
% Qmean3(:,4) = sum(Qmean3(:,1:3),2);
% Bmean1(:,4) = sum(Bmean1(:,1:3),2);
% Bmean2(:,4) = sum(Bmean2(:,1:3),2);
% Bmean3(:,4) = sum(Bmean3(:,1:3),2);

% magnitude of vector at each point
for ii = 1:size(Qmean1,1); Qmean1(ii,4) = norm(Qmean1(ii,1:3)); end
for ii = 1:size(Qmean2,1); Qmean2(ii,4) = norm(Qmean2(ii,1:3)); end
for ii = 1:size(Qmean3,1); Qmean3(ii,4) = norm(Qmean3(ii,1:3)); end
for ii = 1:size(Bmean1,1); Bmean1(ii,4) = norm(Bmean1(ii,1:3)); end
for ii = 1:size(Bmean2,1); Bmean2(ii,4) = norm(Bmean2(ii,1:3)); end
for ii = 1:size(Bmean3,1); Bmean3(ii,4) = norm(Bmean3(ii,1:3)); end

% std across pipe
Qstd1 = std(Qmean1,1);
Qstd2 = std(Qmean2,1);
Qstd3 = std(Qmean3,1);
Bstd1 = std(Bmean1,1);
Bstd2 = std(Bmean2,1);
Bstd3 = std(Bmean3,1);

% sem across pipe
Qsem1 = Qstd1./sqrt(size(Qmean1,1));
Qsem2 = Qstd2./sqrt(size(Qmean2,1));
Qsem3 = Qstd3./sqrt(size(Qmean3,1));
Bsem1 = Bstd1./sqrt(size(Bmean1,1));
Bsem2 = Bstd2./sqrt(size(Bmean2,1));
Bsem3 = Bstd3./sqrt(size(Bmean3,1));

% mean across pipe
Qmean1 = mean(Qmean1,1);
Qmean2 = mean(Qmean2,1);
Qmean3 = mean(Qmean3,1);
Bmean1 = mean(Bmean1,1);
Bmean2 = mean(Bmean2,1);
Bmean3 = mean(Bmean3,1);

% plot
figure;
titles = {'V_x','V_y','V_z','|V|'};
for nn = 1:4
    cmp = nn; % 1 = Vx, 2 = Vy, 3 = Vz, 4 = Vtot
    subplot(2,2,nn);
    bar([Qmean1(cmp) Bmean1(cmp); Qmean2(cmp) Bmean2(cmp); Qmean3(cmp) Bmean3(cmp)]);
    axis([0.5 3.5 -30 30]);
    title(titles{nn})
end

% for graphpad
cmp=3;
gp_mag = [Qmean1(cmp) Bmean1(cmp); Qmean2(cmp) Bmean2(cmp); Qmean3(cmp) Bmean3(cmp)];
gp_std = [Qstd1(cmp) Bstd1(cmp); Qstd2(cmp) Bstd2(cmp); Qstd3(cmp) Bstd3(cmp)];
gp_sem = [Qsem1(cmp) Bsem1(cmp); Qsem2(cmp) Bsem2(cmp); Qsem3(cmp) Bsem3(cmp)];


%% FIGURE: Phase correction 

% fcmr194:
cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr194');

cd data

s16_ab = load_untouch_nii('s16_rlt_ab.nii.gz');
s16_ph = load_untouch_nii('s16_rlt_ph.nii.gz');
s16_poly = load_untouch_nii('s16_rlt_ph_polynomial_uterus_MITKsave.nii.gz'); %weird... original file wouldn't open so I opened it in MITK and then saved it in there.
s16_ph_corr = load_untouch_nii('s16_rlt_ph_corr_uterus.nii.gz');

% re-orient images into nice view
s16_ab.img = flip(flip(double(s16_ab.img),1),2);
s16_ph.img = flip(flip(double(s16_ph.img),1),2);
s16_poly.img = flip(flip(double(s16_poly.img),1),2);
s16_ph_corr.img = flip(flip(double(s16_ph_corr.img),1),2);

x = 50:110;
y = 70:130;
sl = 5;
t = 18; % nice, strong phase change in DAo

implay_RR(s16_ab.img(:,:,sl,:));

% ab / ph / ph_polynomial / ph_corr
cd(paperDir);
imtar(s16_ab.img(x,y,sl,t),[],[],'paperFig');
saveas(gcf,['fcmr194_rlt_ab.png']);

imtar(s16_ph.img(x,y,sl,t),0.75,1.75,'paperFig'); 
saveas(gcf,['fcmr194_rlt_ph.png']);

% imtar(s16_poly.img(x,y,1,1),-3,-1,'paperFig'); %colorbar;
% TOFIX: think this is legit?: Basically how the polynomial is displayed.
% By .*-1 it becomes: ph_corr = ph - poly. In actual pipeline, it's done
% using complex number maths with exponentials
imtar(-1.*s16_poly.img(x,y,1,1),1,2,'paperFig'); %colorbar;
saveas(gcf,['fcmr194_rlt_ph_poly.png']);

imtar(s16_ph_corr.img(x,y,sl,t),-0.5,0.5,'paperFig'); 
saveas(gcf,['fcmr194_rlt_ph_corr.png']);

close all;


% pi/2 windowing
cd(paperDir);
mkdir pi_over2_plus_1p25
cd pi_over2_plus_1p25

rangeW = pi/3;
offset = 1.25;
lowW   = (-rangeW)+offset;
highW  = ( rangeW)+offset;

imtar(s16_ab.img(:,:,sl,t),[],[],'paperFig');
saveas(gcf,['fcmr194_rlt_ab.png']);

imtar(s16_ph.img(:,:,sl,t),lowW,highW,'paperFig'); 
saveas(gcf,['fcmr194_rlt_ph.png']);

imtar(s16_ab.img(x,y,sl,t),[],[],'paperFig');
saveas(gcf,['fcmr194_rlt_ab_cropped.png']);

imtar(s16_ph.img(x,y,sl,t),lowW,highW,'paperFig'); 
saveas(gcf,['fcmr194_rlt_ph_cropped.png']);

% imtar(s16_poly.img(x,y,1,1),-3,-1,'paperFig'); %colorbar;
% TOFIX: think this is legit?: Basically how the polynomial is displayed.
% By .*-1 it becomes: ph_corr = ph - poly. In actual pipeline, it's done
% using complex number maths with exponentials
imtar(-1.*s16_poly.img(x,y,1,1),lowW,highW,'paperFig'); %colorbar;
saveas(gcf,['fcmr194_rlt_ph_poly.png']);

imtar(-1.*s16_ph_corr.img(x,y,sl,t),-rangeW,rangeW,'paperFig'); 
saveas(gcf,['fcmr194_rlt_ph_corr.png']);

close all;
cd(paperDir);