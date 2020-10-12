%% Sweep Windows Visualiser
%
% Visualise/Tinker with Sweep windowing
% - Code taken from mrecon_ktsweep_window.m
%
%

%% Window Configuration
nZ   = 11;      % slices
nT   = 96;      % dynamics
nZnT = nZ*nT;

swpWinWidths = 96;
swpWinStride = 32;
swpWinOffset = 0;
swpWinLoca = swpWinOffset+swpWinWidths/2 : swpWinStride : nZnT;

swpWinFullWidth = swpWinWidths;


%% Create Matrix of Sweep Window Indices

clear swpWindows ktAcqSwpWin ktTrnSwpWin

% Window Half Width
swpWinHalfWidth = ceil( swpWinFullWidth / 2 );

% Window Locations (centre points)
if swpWinLoca
    swpWinSpacing = unique( diff( swpWinLoca ) );
elseif isempty( swpWinSpacing )
    swpWinSpacing = swpWinFullWidth; % M2D equivalent
end

if isempty( swpWinLoca )
    swpWinLoca = swpWinHalfWidth:swpWinSpacing:nZnT;
end

numSwpWindows = numel( swpWinLoca );

% Sweep Windows Array
for iW = 1:numSwpWindows
    swpWindows(:,iW) = ...
        swpWinLoca(iW)-swpWinHalfWidth+1:swpWinLoca(iW)+swpWinHalfWidth;
end

% Ensure Windows within Acquisition FOV
[~,swpWinOutOfBounds,~] = find(swpWindows > nZnT);
swpWindows( :, unique(swpWinOutOfBounds) ) = [];
numSwpWindows = size( swpWindows, 2 );


%% Visualise Sweep Windows
    
% Sweep Windows
figure; hold on;
plot(swpWindows','.b');

% M2D Frames
m2dWinLoca = 0:nT:nZnT;
for iW = 1:numel(m2dWinLoca)
    plot(1:numSwpWindows, repmat(m2dWinLoca(iW),1,numSwpWindows),'k--');
end

xlabel('Sweep Window No.'); ylabel('Frame Index');
legend('Sweep Windows','Location','NorthWest');
axis([1 numSwpWindows 1 nZnT]);

% % Save
% hFig = gcf; hFig.Name = strcat( outFilePrefix, '_swp_windows' );
% saveas( hFig, [outputDirPath '/' hFig.Name, '.fig' ] );
% saveas( hFig, [outputDirPath '/' hFig.Name, '.png' ] ); clear hFig;
