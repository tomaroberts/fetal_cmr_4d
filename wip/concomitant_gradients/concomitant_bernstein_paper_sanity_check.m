
%% Gy causing z^2 phase error
% ie: Gx = Gz = 0
% therefore coefficients: B = C = D = 0

fov   = [-16:16].*1e-2;       % 32 cm
thk   = 10e-3;                % 10 mm
venc  = 5e-2;                 % 5 cm/s
B0    = 1.5;                  % 1.5 T
Gmax  = 22e-3;                % 22 mT/m
slew  = 77;                   % 77 T/m/s
r     = Gmax/slew;            % ramp time, seconds
TE    = 9e-3;                 % 9 ms
GAMMA = 2 * pi * 42.577e6;    % [Hz/T]

% trapezoidal lobe approximation:
F = 5e-3; % plateau duration
M = Gmax.^2 * (F + ((2/3)*r) );
M = ( GAMMA / (2 * B0) ) * M;
M = r2d(M);   % think this is in [degrees / m]
M = M * 1e-2; % now in [degrees / cm] which is paper units
disp(['M = ' num2str(M) ' [degrees / cm]']);

% Figure 4 using given parameters:
perr = (fov*1e2).^2 .* 0.499;
figure; plot(fov*1e2,perr);



%% Recreating Figure 4 using equations
fov   = [-16:16].*1e-2;       % 32 cm
dt = 1e-6;                    % step size for trapz integration
TE = 9e-3;                    % [s]
Gmax  = 22e-3;                % 22 mT/m
plateauDur = 1.4e-3;          % [s]
slopeTime  = 0.3e-3;          % [s] --- basically = r. If I use r, the gradient plot looks a bit aesthetically weird because of dt
% slopeTime  = r;
zeroOffset = TE/4;            % [s]
GAMMA = 2 * pi * 42.577e6;    % [Hz/T]
B0    = 1.5;                  % 1.5 T

Gx = 0;
Gy = synth_make_gradient( plateauDur, slopeTime, slopeTime, Gmax, zeroOffset, dt, TE );
Gz = 0;

figure('units','normalized','outerposition',[0 .2 1 0.6]);
subplot(1,2,1);
plot(Gy); title('Gy Gradient');


% coeffs in radians?
% A = ( GAMMA / (2 * B0) ) * trapz( (Gx.*Gx) + (Gy.*Gy) ) * dt;
A = ( GAMMA / (2 * B0) ) * trapz( ((Gx.*Gx) + (Gy.*Gy)) -  -((Gx.*Gx) + (Gy.*Gy)) ) * dt; % bipolar, as in paper
B = ( GAMMA / (8 * B0) ) * trapz( Gz.*Gz ) * dt;        % nb: NOT yet bipolar as in paper
C = -1 * ( GAMMA / (2 * B0) ) * trapz( Gx .* Gz ) * dt; % nb: NOT yet bipolar as in paper
D = -1 * ( GAMMA / (2 * B0) ) * trapz( Gy .* Gz ) * dt; % nb: NOT yet bipolar as in paper

% coeffs in deg?
A = r2d(A);   % [degrees / m^2]
% A = A * 1e-4; % [degrees / cm^2]

fn_phi_c = @( x, y, z ) (A .* ( z .^2 )) + ...
                        (B .* ( x .^2 + y .^2 )) + ...
                        (C .* ( x .* z )) + ...
                        (D .* ( y .* z )); % x/y/z here = world coordinates
                    
phi_c = fn_phi_c( zeros(size(fov)), zeros(size(fov)), fov );

subplot(1,2,2);
plot(fov,phi_c);
xlabel('z [m]');
ylabel(' Phase Error [degrees] ');

disp(['Coefficient A = ' num2str(A * 1e-4) ' [degrees / cm^2]']);
disp(['Coefficient B = ' num2str(B * 1e-4) ' [degrees / cm^2]']);
disp(['Coefficient C = ' num2str(C * 1e-2) ' [degrees / cm]']);
disp(['Coefficient D = ' num2str(D * 1e-2) ' [degrees / cm]']);



