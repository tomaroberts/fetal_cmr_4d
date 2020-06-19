function fcmr_concomitant_correction( reconDir, sID, varargin )
%FCMR_CONCOMITANT_CORRECTION  concomitant gradient correction of phase data
%
%   FCMR_CONCOMITANT_CORRECTION( reconDir, 'param', val )
%       Concomitant gradient correction within MRecon is only applied to
%       phase contrast flow scans. Not performed for bffe scans. This
%       function calculates concomitant gradients at every point across the
%       imaging volume and subtracts from the phase data
% 
%   Requires:
%       
%   Input:
%       reconDir            - str - local directory of fetal subjects
%       sID                 - int - stack ID
%
%   Optional Parameter-value Pairs:
%       dataDir             fcmr data directory
%       ktreconDir          fcmr ktrecon directory
%       sliceToDisplay      slice to display in comparison figure
%       scaleFactor         phase scale factor [for experimenting]
%
%   Example usage:
%       fcmr_4dflow_postprocessing( 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\fcmr194' )
%
%   TODO:
%       - 
%
%   See also:
%       wrapper_concomitant_correction.m, calculate_concomitant_coeffs.m,
%       synth_make_gradient.m
%       

% Tom Roberts (t.roberts@kcl.ac.uk)


%% Parse Input

default.dataDir            = 'data';
default.ktreconDir         = 'ktrecon';
default.scaleFactor        = 1;
default.sliceToDisplay     = 1;

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(  p, 'studyDir' );

add_param_fn( p, 'dataDir', default.dataDir, ...
        @(x) validateattributes( x, {'char'}, ...
        {}, mfilename ) );

add_param_fn( p, 'ktreconDir', default.ktreconDir, ...
        @(x) validateattributes( x, {'char'}, ...
        {}, mfilename ) );
    
add_param_fn( p, 'scaleFactor', default.scaleFactor, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {}, mfilename ) );
    
add_param_fn( p, 'sliceToDisplay', default.sliceToDisplay, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {}, mfilename ) );

parse( p, reconDir, varargin{:} );

dataDir             = p.Results.dataDir;
ktreconDir          = p.Results.ktreconDir;
scaleFactor         = p.Results.scaleFactor;
sliceToDisplay      = p.Results.sliceToDisplay;


%% Load phase data
cd( fullfile( reconDir, dataDir ) );

ph.im      = niftiread( strcat( 's', num2str(sID), '_rlt_ph.nii.gz' ) );
ph.info    = niftiinfo( strcat( 's', num2str(sID), '_rlt_ph.nii.gz' ) );
ph.im_corr = niftiread( strcat( 's', num2str(sID), '_rlt_ph_corr_uterus.nii.gz' ) ); % 3D polynomial corrected
ab.im      = niftiread( strcat( 's', num2str(sID), '_rlt_ab.nii.gz' ) );

% cd( fullfile( reconDir, ktreconDir ) );
% re.im = niftiread( strcat( 's', num2str(sID), '_rlt_re.nii.gz' ) );
% im.im = niftiread( strcat( 's', num2str(sID), '_rlt_im.nii.gz' ) );
% 
% cx.im = re.im + 1i.*im.im;
% ph.im = angle(cx.im);
% % ph.im = angle(abs(cx.im).*exp(sqrt(-1)*(angle(cx.im)+pi))); %DO I NEED +pi in exp? Not sure

cd( fullfile( reconDir, ktreconDir ) );
load( strcat( 's', num2str(sID), '_rlt_parameters.mat' ), 'PARAM' );


%% Load body coil
load( ['s', num2str(sID), '_csm.mat'], 'imBody' );
% load( ['s', num2str(sID), '_csm.mat'], 'imCoil' );

% Resample body coil image
dimX = 1;
fn_rmv_os = @( x ) x((round(size(x,dimX)/4)+1):(size(x,dimX)-round(size(x,dimX)/4)),:,:,:,:);
im_bdy = fn_rmv_os( imBody );
% ph.im_bdy = angle( imBody ); clear imBody;
ph.im_bdy = angle(abs(im_bdy).*exp(sqrt(-1)*(angle(im_bdy)+pi)));
ph.im_bdy = squeeze( ph.im_bdy );

% imCoil = fn_rmv_os( imCoil );
% % ph.im_bdy = angle( imBody ); clear imBody;
% ph.im_coil = angle(abs(imCoil).*exp(sqrt(-1)*(angle(imCoil)+pi)));
% ph.im_coil = squeeze( ph.im_coil ); ph.im_coil = permute( ph.im_coil, [1 2 4 3] );


%% Determine axes orientation

iLPF = PARAM.Scan.ijk(1:2);
jLPF = PARAM.Scan.ijk(4:5);
kLPF = PARAM.Scan.ijk(7:8);

iXYZ = philips_lpf2xyz(iLPF);
jXYZ = philips_lpf2xyz(jLPF);
kXYZ = philips_lpf2xyz(kLPF);

imtar( ph.im(:,:,5,1) ); 
title([ 'sID = ' num2str(sID) ' --- i = ' iLPF ' / ' iXYZ ...
                                 ', j = ' jLPF ' / ' jXYZ ...
                                 ', k = ' kLPF ' / ' kXYZ ]);
close;


%% Extract voxel > world coordinate affine matrix
% Based on Method 3: https://brainder.org/2012/09/23/the-nifti-file-format/ 

T = ph.info.Transform.T';
L = [T(1:3,1)/2 T(1:3,2)/2 T(1:3,3)/4]'; % image lattice (= mirtk info)

dimX      = ph.info.raw.dim(2); 
dimY      = ph.info.raw.dim(3);
numSlices = ph.info.raw.dim(4);
dimT      = ph.info.raw.dim(5);

voxelRange.i = [0:dimX-1]';
voxelRange.j = [0:dimY-1]';
voxelRange.k = [0:numSlices-1]';

[voxelGrid.X, voxelGrid.Y, voxelGrid.Z] = ...
    meshgrid( voxelRange.i, voxelRange.j, voxelRange.k );

% TODO: vectorise to make faster
for ii = 1:numel(voxelRange.j)
    for jj = 1:numel(voxelRange.i)
        for kk = 1:numel(voxelRange.k)
            
            currVoxelCoord = [voxelGrid.X(ii,jj,kk); ...
                              voxelGrid.Y(ii,jj,kk); ...
                              voxelGrid.Z(ii,jj,kk); ...
                              1];
            
            % apply transformation to voxel location
            currWorldCoord = T * currVoxelCoord;
            
%             % IMPORTANT: i/j swapped in MATLAB
%             worldCoords.X(jj,ii,kk) = currWorldCoord(1) .* 1e-3; % [metres]
%             worldCoords.Y(jj,ii,kk) = currWorldCoord(2) .* 1e-3; % [metres]
%             worldCoords.Z(jj,ii,kk) = currWorldCoord(3) .* 1e-3; % [metres]
            
            worldCoords.X(ii,jj,kk) = currWorldCoord(1) .* 1e-3; % [metres]
            worldCoords.Y(ii,jj,kk) = currWorldCoord(2) .* 1e-3; % [metres]
            worldCoords.Z(ii,jj,kk) = currWorldCoord(3) .* 1e-3; % [metres]
            
        end
    end
end


%% Function to convert voxel > world coordinates

% * 1e3 to put in [mm] for ease of comparison with rview

% Configuration below of getWorldCoords matches rview RADIOLOGICAL orientation
getWorldCoords = @(i,j,k) [worldCoords.X(j,i,k) * 1e3; ...
                           worldCoords.Y(j,i,k) * 1e3; ...
                           worldCoords.Z(j,i,k) * 1e3; ...
                           1];
% i = row
% j = col
% k = slice
% e.g.: for fcmr194, s16:
% ___ MATLAB ____                     = ___ rview ____
% getWorldCoords(row=114,col=40,sl=5) = [-99.76, 41.95, -2.39];
% ph.im(114,40,5,1)                   = 1.67;                


%% Calculate Concomitant Coefficients
[COEFFS, O] = calculate_concomitant_coeffs( reconDir, sID, L );
% close;

A = COEFFS.A; B = COEFFS.B; C = COEFFS.C; D = COEFFS.D;

% [Lnew, O2L] = O_matrix2niiLattice(O);

% Ot = O';
% 
% Ot(:,1) = Ot(:,1) * ph.info.PixelDimensions(1);
% Ot(:,2) = Ot(:,2) * ph.info.PixelDimensions(2);
% Ot(:,3) = Ot(:,3) * ph.info.PixelDimensions(3);


%% Calculate Concomitant Phase Error
fn_phi_c = @( x, y, z ) A .* ( z .^2 ) + ...
                        B .* ( x .^2 + y .^2 ) + ...
                        C .* ( x .* z ) + ...
                        D .* ( y .* z ); % x/y/z here = world coordinates

% THIS IS SO STUPID:
x_axis = philips_determine_order(iXYZ);
y_axis = philips_determine_order(jXYZ);
z_axis = philips_determine_order(kXYZ);
                    
% This long-winded method is definitely matching voxelCoord to worldCoord:
clear phi_c

ctr = 1;
                    
for ii = 1:numel(voxelRange.j)
    for jj = 1:numel(voxelRange.i)
        for kk = 1:numel(voxelRange.k)

            currVoxelCoord = [voxelGrid.X(ii,jj,kk)+1; ...
                              voxelGrid.Y(ii,jj,kk)+1; ...
                              voxelGrid.Z(ii,jj,kk)+1; ...
                              1];
            
            currWorldCoord = 1e-3 .* getWorldCoords(currVoxelCoord(1),... % 1e-3: [mm] > [m]
                                                    currVoxelCoord(2),...
                                                    currVoxelCoord(3));
            
            currWorldCoord(4) = 1;
                                                
            vxRecord(:,ctr) = currVoxelCoord;
            wrRecord(:,ctr) = currWorldCoord;
            ctr = ctr + 1;
            
            % Philips scanner system directions:
            % +x = PA
            % +y = RL
            % +z = FH
            % So... fn_phi_c( -currWorldCoord(3),-currWorldCoord(2),+currWorldCoord(1) )
            %             ===> i = -z // j = -y // k = +x
            
            curr_phi_c = fn_phi_c( sign(x_axis)*currWorldCoord(abs(x_axis)), ...
                                   sign(y_axis)*currWorldCoord(abs(y_axis)), ...
                                   sign(z_axis)*currWorldCoord(abs(z_axis)) );
            
            phi_c(currVoxelCoord(1),currVoxelCoord(2),currVoxelCoord(3)) = curr_phi_c;
                            
        end
    end
end

phi_c = phi_c * scaleFactor;
% ph_com.im = (ph.im+0) - (phi_c);


%% Think I have broken this method because I changed from:
% worldCoords.X(jj,ii,kk) to worldCoords.X(ii,jj,kk) in loop up above

% phi_c = fn_phi_c( worldCoords.X, worldCoords.Y, worldCoords.Z );
% phi_c = fn_phi_c( worldCoords.Y, worldCoords.X, worldCoords.Z );
% % phi_c = fn_phi_c( -worldCoords.Z, worldCoords.Y, worldCoords.X );
% 
% % phi_c = phi_c * (2*pi);
% % ph_com.im = (ph.im - 0.8) - phi_c;
% 
% % phi_c = phi_c;
% phi_c = phi_c * 10;
% ph_com.im = (ph.im+0) - phi_c;


%% Subtract Body Coil phase
ph.im_bdy_corr = ph.im - ph.im_bdy;
ph.im_bdy_com_corr = ph.im - ph.im_bdy - phi_c;


%% View corrected phase data and phase profiles

sl  = sliceToDisplay;
dyn = 1;
% profileRow = round(dimX/2);
profileRow = 120;

% colorbar range
cbR = pi/2;

% normalise absolute image
abNormFactor = max(max(ab.im(:,:,sl,dyn)));


figure('units','normalized','outerposition',[0 0 1 1]);

% subplot1 - phase images
subplot(1,2,1);
imagesc([ (2.*cbR .* (ab.im(:,:,sl,dyn)./abNormFactor) - cbR), ...
           ph.im(:,:,sl,dyn), ...
           phi_c(:,:,sl,dyn); ...
           ph.im_bdy_corr(:,:,sl,dyn), ...
           ph.im_bdy_com_corr(:,:,sl,dyn), ...
           ph.im_corr(:,:,sl,dyn)], ...
           [-cbR,cbR]);
axis('image');
colormap('gray'); colorbar;
title([ 'sID = ' num2str(sID) ' --- i = ' iLPF ' / ' iXYZ ...
                                 ', j = ' jLPF ' / ' jXYZ ...
                                 ', k = ' kLPF ' / ' kXYZ ]);

% subplot2 - phase profile                    
subplot(1,2,2);
hold on
plot( ph.im_bdy_corr(profileRow,:,sl,dyn), '-b' );
plot( ph.im_bdy_com_corr(profileRow,:,sl,dyn), '-r'  );
plot( ph.im_corr(profileRow,:,sl,dyn), '-k' );
axis([1,dimY,-cbR,cbR]);
axis('square');
title(['Profiles across row number: ' num2str(profileRow)]);
legend('Body Coil Corrected','Body Coil + Concomitant Corrected','Polynomial Corrected');


%% View subtraction

sl  = sliceToDisplay;
dyn = 1;
% profileRow = round(dimX/2);
% profileRow = 125;

% colorbar range
cbR = pi;

% normalise absolute image
abNormFactor = max(max(ab.im(:,:,sl,dyn)));


figure('units','normalized','outerposition',[0 0 1 1]);

% subplot1 - phase images
subplot(1,2,1);
imagesc([ (2.*cbR .* (ab.im(:,:,sl,dyn)./abNormFactor) - cbR), ...
           ph.im(:,:,sl,dyn); ...
           ph.im_bdy(:,:,sl,dyn), ...
           ph.im_bdy_corr(:,:,sl,dyn)], ...
           [-cbR,cbR]);
axis('image');
colormap('gray'); colorbar;
title([ 'sID = ' num2str(sID) ' --- i = ' iLPF ' / ' iXYZ ...
                                 ', j = ' jLPF ' / ' jXYZ ...
                                 ', k = ' kLPF ' / ' kXYZ ]);

% subplot2 - phase profile                    
subplot(1,2,2);
hold on
plot( ph.im(profileRow,:,sl,dyn), '-k' );
plot( ph.im_bdy(profileRow,:,sl,dyn), '-r' );
plot( ph.im_bdy_corr(profileRow,:,sl,dyn), '-b' );
axis([1,dimY,-cbR,cbR]);
axis('square');
title(['Profiles across row number: ' num2str(profileRow)]);
legend('Acquisition Phase','Body Coil Phase','Acq - Body');


%% Body coil correction vs. polynomial
% - figure to show Jo
% imtar([ph.im_bdy_corr(:,:,sl,dyn)-.4, ph.im_corr(:,:,sl,dyn)],-.5,.5);

% fn end
end

