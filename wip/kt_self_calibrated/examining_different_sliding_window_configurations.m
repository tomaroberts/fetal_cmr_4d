cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr254_tom_recon_hrh\ktrecon');

%% C = [8 4 3 2 1 1]
load('s20_x_C843211_Lucilio_recon.mat');

implay_RR([squeeze(abs(x(:,:,8,1,:))),...
           squeeze(abs(x(:,:,8,2,:))),...
           squeeze(abs(x(:,:,8,3,:))),...
           squeeze(abs(x(:,:,8,4,:))),...
           squeeze(abs(x(:,:,8,5,:))),...
           squeeze(abs(x(:,:,8,6,:)))],...
                       'gray',[0, 1400]);
                   
                   
dimT = 5;
xt2xf = @( xt ) fftshift( fft( xt, [], dimT ), dimT );

xf = xt2xf(x);

xfLine=62;
imtar([squeeze(abs(xf(:,xfLine,8,1,:))),...
       squeeze(abs(xf(:,xfLine,8,2,:))),...
       squeeze(abs(xf(:,xfLine,8,3,:))),...
       squeeze(abs(xf(:,xfLine,8,4,:))),...
       squeeze(abs(xf(:,xfLine,8,5,:))),...
       squeeze(abs(xf(:,xfLine,8,6,:)))],0,5000); colormap('parula');
   
   
   
%% C = [8 7 6 5 4]
load('s20_x_C87654_Lucilio_recon.mat');

implay_RR([squeeze(abs(x(:,:,8,1,:))),...
           squeeze(abs(x(:,:,8,2,:))),...
           squeeze(abs(x(:,:,8,3,:))),...
           squeeze(abs(x(:,:,8,4,:))),...
           squeeze(abs(x(:,:,8,5,:)))],...
                       'gray',[0, 1400]);
                   
                   
dimT = 5;
xt2xf = @( xt ) fftshift( fft( xt, [], dimT ), dimT );

xf = xt2xf(x);

xfLine=60;
imtar([squeeze(abs(xf(:,xfLine,8,1,:))),...
       squeeze(abs(xf(:,xfLine,8,2,:))),...
       squeeze(abs(xf(:,xfLine,8,3,:))),...
       squeeze(abs(xf(:,xfLine,8,4,:))),...
       squeeze(abs(xf(:,xfLine,8,5,:)))],0,5000); colormap('parula');
   
   
%% C = [8 5 3]
load('s20_x_C853_Lucilio_recon.mat');

implay_RR([squeeze(abs(x(:,:,8,1,:))),...
           squeeze(abs(x(:,:,8,2,:))),...
           squeeze(abs(x(:,:,8,3,:))),...
                       ],'gray',[0, 1400]);
                   
                   
dimT = 5;
xt2xf = @( xt ) fftshift( fft( xt, [], dimT ), dimT );

xf = xt2xf(x);

xfLine=62;
imtar([squeeze(abs(xf(:,xfLine,8,1,:))),...
       squeeze(abs(xf(:,xfLine,8,2,:))),...
       squeeze(abs(xf(:,xfLine,8,3,:)))],0,5000); colormap('parula');
   
   
%% C = 4
load('s20_x_C4_Lucilio_recon.mat');

implay_RR([squeeze(abs(x(:,:,8,1,:))),...
                       ],'gray',[0, 1400]);
                   
                   
dimT = 5;
xt2xf = @( xt ) fftshift( fft( xt, [], dimT ), dimT );

xf = xt2xf(x);

xfLine=62;
imtar([squeeze(abs(xf(:,xfLine,8,1,:)))],0,5000); colormap('parula');