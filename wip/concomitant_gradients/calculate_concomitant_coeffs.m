function [COEFFS, O] = calculate_concomitant_coeffs( reconDir, sID, L )

% L = lattice matrix from nifti header --- import from nifti header to
% apply world coordinates to gradient directions (?)

%% fcmr_gradient_waveforms

ktreconDir = 'ktrecon';

cd( fullfile( reconDir, ktreconDir ) );

if exist( strcat( 's', num2str(sID), '_goalc.txt' ), 'file' )
    gcFilename = strcat( 's', num2str(sID), '_goalc_TAR.txt' );
else
    gcFilename = strcat( 's', num2str(sID), '_goalc.txt' );
end

gc = get_pcmr_orientation_parameters( gcFilename );

% get M / S gradient strengths, durations and slopes
GRADS.mc0.str  = gc.mc0_str .* 1e-3; % [T/m]
GRADS.m0.str   = gc.m0_str .* 1e-3;
GRADS.m3.str   = gc.m3_str .* 1e-3;
GRADS.py.str   = gc.py_str .* 1e-3;
GRADS.s_ex.str = gc.s_ex_str .* 1e-3;
GRADS.r0.str   = gc.r0_str .* 1e-3;

GRADS.mc0.dur  = gc.mc0_dur .* 1e-3; % [seconds]
GRADS.m0.dur   = gc.m0_dur .* 1e-3;
GRADS.m3.dur   = gc.m3_dur .* 1e-3;
GRADS.py.dur   = gc.py_dur .* 1e-3;
GRADS.s_ex.dur = gc.s_ex_dur .* 1e-3;
GRADS.r0.dur   = gc.r0_dur .* 1e-3;

GRADS.mc0.slope1  = gc.mc0_slope1 .* 1e-3; % [seconds]
GRADS.m0.slope1   = gc.m0_slope1 .* 1e-3;
GRADS.m3.slope1   = gc.m3_slope1 .* 1e-3;
GRADS.py.slope1   = gc.py_slope1 .* 1e-3;
GRADS.s_ex.slope1 = gc.s_ex_slope1 .* 1e-3;
GRADS.r0.slope1   = gc.r0_slope1 .* 1e-3;

GRADS.mc0.slope2  = gc.mc0_slope2 .* 1e-3; % [seconds]
GRADS.m0.slope2   = gc.m0_slope2 .* 1e-3;
GRADS.m3.slope2   = gc.m3_slope2 .* 1e-3;
GRADS.py.slope2   = gc.py_slope2 .* 1e-3;
GRADS.s_ex.slope2 = gc.s_ex_slope2 .* 1e-3;
GRADS.r0.slope2   = gc.r0_slope2 .* 1e-3;

GRADS.mc0.time  = gc.mc0_time .* 1e-3; % [seconds]
GRADS.m0.time   = gc.m0_time .* 1e-3;
GRADS.m3.time   = gc.m3_time .* 1e-3;
GRADS.py.time   = gc.py_time .* 1e-3;
GRADS.s_ex.time = gc.s_ex_time .* 1e-3;
GRADS.r0.time   = gc.r0_time .* 1e-3;

GRADS.mc0.ref  = gc.mc0_ref .* 1e-3; % [seconds]
GRADS.m0.ref   = gc.m0_ref .* 1e-3;
GRADS.m3.ref   = gc.m3_ref .* 1e-3;
GRADS.py.ref   = gc.py_ref .* 1e-3;
GRADS.s_ex.ref = gc.s_ex_ref .* 1e-3;
GRADS.r0.ref   = gc.r0_ref .* 1e-3;

% get orient objects for record
m_orient = gc.m_orient;
p_orient = gc.p_orient;
s_orient = gc.s_orient;
O = [m_orient; p_orient; s_orient];
    


%% Setup

TE = 1.910e-3; % [seconds]
dt = 1e-6;
T  = linspace(0,TE,TE/dt);

figure('units','normalized','outerposition',[0 0 1 1]);


%% MEASUREMENT (Readout)

G_mc0 = synth_make_gradient( GRADS.mc0.dur, ...
                             GRADS.mc0.slope1, ...
                             GRADS.mc0.slope2, ...
                             GRADS.mc0.str, ...
                             GRADS.mc0.time - GRADS.mc0.ref, ...
                             dt, ...
                             TE );
                          
G_m0  = synth_make_gradient( GRADS.m0.dur, ...
                             GRADS.m0.slope1, ...
                             GRADS.m0.slope2, ...
                             GRADS.m0.str, ...
                             GRADS.m0.time - GRADS.m0.ref, ...
                             dt, ...
                             TE );

Gm = G_m0 + G_mc0;
                          
subplot(2,2,1);
hold on; 
title('READOUT');
plot(T,G_mc0);
plot(T,G_m0);
plot(T,Gm,'k--');
xlabel('Time [s]');
ylabel('G [T]');
legend('mc[0]','m[0]','Difference');

                          
%% PHASE

% custom PE gradient strength (measured from GVE, or GR`py:str_step * GR`py:str_factor_max)
% GRADS.py.str = 8.546;
                          
G_py = synth_make_gradient( GRADS.py.dur, ...
                            GRADS.py.slope1, ...
                            GRADS.py.slope2, ...
                            GRADS.py.str, ...
                            GRADS.py.time - GRADS.py.ref, ...
                            dt, ...
                            TE );                      

Gp = G_py;
                         
subplot(2,2,2);
hold on; 
title('PHASE');
plot(T,G_py);
plot(T,Gp,'k--');
xlabel('Time [s]');
ylabel('G [T]');
legend('py','Difference');


%% SLICE SELECT
                          
G_s_ex = synth_make_gradient( GRADS.s_ex.dur, ...
                              GRADS.s_ex.slope1, ...
                              GRADS.s_ex.slope2, ...
                              GRADS.s_ex.str, ...
                              GRADS.s_ex.time - GRADS.s_ex.ref, ...
                              dt, ...
                              TE );
                           
G_r0 = synth_make_gradient( GRADS.r0.dur, ...
                            GRADS.r0.slope1, ...
                            GRADS.r0.slope2, ...
                            GRADS.r0.str, ...
                            GRADS.r0.time - GRADS.r0.ref, ...
                            dt, ...
                            TE );                         

Gs = G_s_ex + G_r0;
                         
subplot(2,2,3);
hold on; 
title('SLICE SELECT');
plot(T,G_s_ex);
plot(T,G_r0);
plot(T,Gs,'k--');
xlabel('Time [s]');
ylabel('G [T]');
legend('s-ex','r[0]','Difference');


%% Convert waveforms MPS > XYZ

Gmps = [Gm; Gp; Gs];

% Gxyz = O' * Gmps; % o_matrix from goalc (ijk to xyz)
Gxyz = L' * Gmps; % lattice from nifti (ie: world)

% Onii = O_matrix2niiLattice(O);
% Gxyz = Onii * Gmps;

Gx = Gxyz(1,:);
Gy = Gxyz(2,:);
Gz = Gxyz(3,:);

%%% NB: might need wf in world coordinates rather than xyz
% Cprime = [0 -1 0; 1 0 0; 0 0 1];
% Gworld = Cprime * Gxyz;
% 
% Gx = Gworld(1,:);
% Gy = Gworld(2,:);
% Gz = Gworld(3,:);

% figure;
subplot(2,2,4);
hold on; 
title('XYZ waveforms');
plot(T,Gx);
plot(T,Gy);
plot(T,Gz);
plot(T,zeros(size(Gx)),'k--');
xlabel('Time [s]');
ylabel('G [T]');
legend('G_x','G_y','G_z');


%% Calculate Concomitant Coefficients

% Gx = Gx * 1e-3; % [T/m]
% Gy = Gy * 1e-3; % [T/m]
% Gz = Gz * 1e-3; % [T/m]

% dt = dt * 1e-3; % [seconds]

GAMMA = 2 * pi * 42.577e6; % [Hz/T]
B0    = 1.5;        % [T]

A = ( GAMMA / (2 * B0) ) * trapz( (Gx.*Gx) + (Gy.*Gy) ) * dt; % * dt tells trapz spacing increment

B = ( GAMMA / (8 * B0) ) * trapz( Gz.*Gz ) * dt;

C = -1 * ( GAMMA / (2 * B0) ) * trapz( Gx .* Gz ) * dt;

D = -1 * ( GAMMA / (2 * B0) ) * trapz( Gy .* Gz ) * dt;

COEFFS.A = A;
COEFFS.B = B;
COEFFS.C = C;
COEFFS.D = D;


% fn end
end




















