%% Script for trying to see Concomitant Gradients effect in phantom data
% February 2020

% phantom 8 study
% ax_ph appears to have concomitant phase roll

phantomDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_04_Flow_Phantom8\PC-bSSFP\data';
cd(phantomDir);


% %%% Shen load_nii
% ax.ph = load_nii('ax_ph.nii.gz');
% ax.ph_corr = load_nii('ax_ph_corr.nii.gz');

%% MATLAB niftiread
ax.ph = niftiread('ax_ph.nii.gz');
ax.ph_corr = niftiread('ax_ph_corr.nii.gz');

ax.ph_info = niftiinfo('ax_ph.nii.gz');
ax.ph_corr_info = niftiinfo('ax_ph_corr.nii.gz');

% LR = PA / UD = SI >>> sagittal to bore >>> rows = SI, cols = PA
ax.ph = flip(rot90(permute(ax.ph,[1,3,2])),2); 
ax.ph_corr = flip(rot90(permute(ax.ph_corr,[1,3,2])),2); 

phSlice      = ax.ph(:,:,60);
ph_corrSlice = ax.ph_corr(:,:,60);

imtar(phSlice); axis square;
[cx1,cy1,ax.ph_C,xi1,yi1]      = improfile;
[cx2,cy2,ax.ph_corr_C]         = improfile(ph_corrSlice,xi1,yi1);
close all;

% number of points along improfile line 
numPoints = numel(ax.ph_C);

figure; 
hold on;
plot(ax.ph_C,'ro');
plot(ax.ph_corr_C,'bo');



%% simulation

% approx world location of line profile (read from rview)
point1 = [43.6 -21.1 15.6] * 1e-3;
point2 = [43.6 -21.1 -96.6] * 1e-3;

% max gradient amplitudes taken from GVE (positive lobes only so far)

% TODO: establish which gradients should be used in calculation. I've
% simply included them all, essentially as flat gradients on permanently

% Phantom 8
Gm  = 17.95 * 1e-3;      % [mT/m]
tGm = 1.06e-3;           % [sec]

Gp  = 5.79 * 1e-3;       % [mT/m]
tGp = 1.28e-3;           % [sec]

Gs  = 17.75 * 1e-3;      % [mT/m]
tGs = 1.15e-3;           % [sec]

% fcmr protocol
point1 = [43.6 -21.1 -200] * 1e-3;
point2 = [43.6 -21.1 200] * 1e-3;

Gm  = 11.65 * 1e-3;      % [mT/m]
Gp  = 8.54 * 1e-3;       % [mT/m]
Gs  = 12.68 * 1e-3;      % [mT/m]

G = [Gm; Gp; Gs];

% GAMMA = 42577;   % [Hz/mT] gyromagnetic ratio    %%%%%%%% * 2 * pi * 42577; %Hz/mT
% z     = 20e-2;   % [m] magnet coordinates
% tG    = 2e-3;    % gradient time --- currently estimated
% G     = 40e-3;   % [T/m] gradient strength
B0    = 1.5;     % main magnetic field

% z = linspace(0,20e-2,numPoints);
% 
% phi_c = (GAMMA * z.^2 * G.^2 * tG) / (2 * B0);
% 
% figure;
% plot(z,phi_c,'o');

% term1 = Gx*z - ( (Gz*x)/2 );
% term2 = Gy*z - ( (Gz*y)/2 );
fn_Bc = @(x, y, z, Gx, Gy, Gz) ( 1/(2*B0) * ( ( Gx*z - ( (Gz*x)/2 ) ).^2 + ( Gy*z - ( (Gz*y)/2 ) ).^2 ) );

worldCoords.x = point1(1) * ones(numPoints,1);
worldCoords.y = point1(2) * ones(numPoints,1);
worldCoords.z = linspace(point1(3),point2(3),numPoints)';

Bc = fn_Bc( worldCoords.x, worldCoords.y, worldCoords.z, Gm, Gp, Gs );

% Bc_ppm = Bc * 1e6;

figure;
plot(1e3*worldCoords.z,Bc,'o');
xlabel('Position (mm)');
ylabel('Concomitant Gradient')


%% Total B
% Don't think this is helpful actually...
B = B0 + ...
    ( G(1) .* worldCoords.x ) + ...
    ( G(2) .* worldCoords.y ) + ...
    ( G(3) .* worldCoords.z ) + ...
    Bc;

% B_ppm = B / 1e6;

figure;
plot(1e3*worldCoords.z,B,'o');
xlabel('Position (mm)');
ylabel('Concomitant Gradient')