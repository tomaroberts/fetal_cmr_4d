function G = synth_make_gradient( dur, slope1, slope2, str, tstart, dt, tfinal  )

% TODO: fix colon operator warning to do with integers
warning('off');

np_fullgrad = dur/dt;
np_slope1   = slope1/dt;
np_slope2   = slope2/dt;
np_plateau  = np_fullgrad - np_slope1 - np_slope2;

G( 1:np_slope1 ) = linspace( 0, str, np_slope1 );
G( np_slope1+1 : np_slope1+np_plateau ) = str;
G( np_slope1+np_plateau+1 : np_fullgrad ) = linspace( str, 0, numel(np_slope1+np_plateau+1 : np_fullgrad) );

% Shift gradient to RF
start_point = tstart/dt;
if start_point > 0
    offset = zeros(size(1:start_point));
    G = cat(2, offset, G );
elseif start_point < 0
    G = G(abs(start_point):end);
end

% % Set offset from RF
% offset      = zeros(size(1:start_point));
% G = cat(2, offset, G );

% Full range
np_total = tfinal/dt;
Gfull = zeros(size(1:np_total));
Gfull(1:numel(G)) = G;

% Return full section of sequence
G = Gfull(1:np_total);

warning('on');


% fn end
end