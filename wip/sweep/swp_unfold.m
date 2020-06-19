function nii3d = swp_unfold( niiFile, outNiiFileName )

%% TODO:
%
% make sure geometry information correct (i.e.: .nii occupies correct
% volume depending on Sweep parameters)
%
%


%% Quick function to convert a Sweep 4D nifti into a 3D nifti

nii4d = load_untouch_nii( niiFile );

nX = size( nii4d.img, 1 );
nY = size( nii4d.img, 2 );
n3 = size( nii4d.img, 3 ) * size( nii4d.img, 4 );

% create 3d nii
nii3d     = nii4d;
nii3d.img = reshape( permute(nii4d.img,[1,2,4,3]), nX, nY, n3 );

% adjust 4d parameters to 3d parameters
nii3d.hdr.dime.dim(1) = 3;
nii3d.hdr.dime.dim(4) = n3;
nii3d.hdr.dime.dim(5) = 1;

% make pixdim 3d
nii3d.hdr.dime.pixdim(1) = 1;
nii3d.hdr.dime.pixdim(4) = nii4d.hdr.dime.pixdim(5);
nii3d.hdr.dime.pixdim(5) = 1;

% TODO: possibly? update affine matrix
% nii3d.hdr.hist.srow_x/y/z etc.

if nargin == 1
    save_untouch_nii( nii3d, [niiFile(1:end-7) '_swp3d.nii.gz' ] );
else
    save_untouch_nii( nii3d, outNiiFileName );
end


end