function imResized = imresize3_oldMatlab( im, dim, interpMethod )

%% Function to recreate imresize3 for older versions of Matlab, pre-2017a
%
%   INPUT:
%       im                  3D input image
%       dim                 [Nx3] row vector - desired output 3D image size
%       interpMethod        interpolation method, default = linear
%
%   OUTPUT:
%       imResized           Resized 3D image
%
% Tom Roberts (t.roberts@kcl.ac.uk)

if nargin < 3
    interpMethod = 'cubic'; % default in imresize3
end

[y, x, z] = ndgrid( linspace( 1, size(im,1), dim(1) ),...
                    linspace( 1, size(im,2), dim(2) ),...
                    linspace( 1, size(im,3), dim(3) ) );

imResized = interp3( im, x, y, z, interpMethod );


% fn end
end