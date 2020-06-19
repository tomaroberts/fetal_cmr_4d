%% fcmr_4d_petable
%
% Visualise PE table of fcmr kt acquisition
%
%

% fcmr194
studyDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr194';
stackNum = 14;

% % Adult Mar2019
% studyDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Adult\2019_03_28_Adult_Heart';
% stackNum  = 10;

ktreconDir = 'ktrecon';
paramFile = ['s' num2str(stackNum) '_rlt_parameters.mat'];

cd(studyDir);

%% Param file
load( fullfile( studyDir, ktreconDir, paramFile ) );


%% ktrecon

numCoils   = numel(PARAM.Labels.CoilNrsPerStack{1,1}); % better way of getting this?
ktFactor   = PARAM.Scan.KtFactor;
ky         = PARAM.Labels.Index.ky_label;
maxKy      = max( ky ) + 1;
numKy      = maxKy / ktFactor;

% reduced Ky (without individual coil entries)
redKy      = PARAM.Labels.Index.ky_label(1:numCoils:end);
redKyDiff  = diff( redKy );
startKyTrn = find( redKyDiff==1, 1, 'first' );

% slice location info
offCentre  = PARAM.Labels.Offcentre;


%% Plot
numPE = ktFactor; % number of kspace PE tables to display

kyAcqRange2plot = 1:(numKy * numPE) + 1; % +1 because PE table weirdly starts in centre for first k-line
kyTrnRange2plot = startKyTrn : startKyTrn + (numKy * numPE) - 1;

figuremax;
% Acquisition
subplot(2,1,1);
plot( redKy(kyAcqRange2plot), 'o' );
ax = axis;
title(['Acquisition data: stack no = ' num2str(stackNum) ]);
xlabel('time (TR number)');
ylabel('PE number');
figurefonts;

% Training
subplot(2,1,2);
plot( redKy(kyTrnRange2plot), 'o' );
axis(ax);
title('Training data');
xlabel('time (TR number)');
ylabel('PE number');
figurefonts;



