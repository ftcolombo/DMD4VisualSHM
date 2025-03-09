clc
clear all
close all

addpath('..\Dados\plastic beam - high speed camera\videos 2022-04-05')
% video records included in the dataset by Garrido et al.
% (https://doi.org/10.1016/j.ymssp.2023.110539)

addpath('../Functions')
addpath('../Data')

vidObj = VideoReader('00_00_00_run1.mp4');

vidWidth = vidObj.Width;
vidHeight = vidObj.Height;

vframes = zeros(vidHeight, vidWidth, 10, 'uint8');

numframes = 0;
vidObj.CurrentTime = 0;
tic
while hasFrame(vidObj)
	numframes = numframes + 1;
    frame = readFrame(vidObj);
    vframes(:,:,numframes) = rgb2gray(frame);
end
toc

% restricting the size of frames to only include the motion of the beam:
% nV is selected by inspection for each video in the dataset
nV = 391:600; %00_00_00_run1.mp4

roif = vframes(nV,:,:);

clearvars -except roif

%% Collecting data for the snaphot matrix

fs = 1000; % frame rate in Hz
T = 2; % analysis of only the first 2 s of the recorded video
ts = 0:1/fs:T;

xy = roif(:,:,1:1:T*fs+1);

clear roif

nheight = size(xy,1);
nwidth = size(xy,2);
nframes = size(xy,3);

npixels = nheight*nwidth;

X = zeros(npixels,nframes);
for i = 1:nframes
    X(:,i) = reshape(xy(:,:,i),[npixels,1]);
end

clear xy npixels nframes

%% Applying optDMD with e_opt_0 as initial guesses for the eigenvalues

% you can run the script 'expvideo_initial_guess4optdmd' to obtain the 
% frequency-stable eigenvalues estimated by DMD for order r = 200, if you 
% have access to the % video '00_00_00_run1.mp4' of the beam by 
% Garrido et al. % (https://doi.org/10.1016/j.ymssp.2023.110539)
load 'initial_guess4optdmd_data_00_00_00_run1I.mat'

addpath('..\optdmd-master\src')
% to run the comand optdmd, the repository https://github.com/duqbo/optdmd
% is required (https://doi.org/10.1137/M1124176)

[w_opt,e_opt2,b_opt] = optdmd(X,ts,length(e_opt_0),2,varpro_opts('ifprint',0),e_opt_0);

%% Comparing the natural frequencies and decay rates computed by DMD with time embeddings and optDMD

% % % Table 4

f_optdmd = abs(e_opt2)/(2*pi);
zeta_optdmd = -real(e_opt2)./abs(e_opt2);
dmd_opt = [f_optdmd zeta_optdmd];

f_dmd = abs(e_opt_0)/(2*pi);
zeta_dmd = -real(e_opt_0)./abs(e_opt_0);
dmd = [f_dmd zeta_dmd];

%% Reconstructing the response with optDMD and DMD

ts = 0:1/fs:0.6;

x_opt_dmd = w_opt*diag(b_opt)*exp(e_opt2*ts);

clearvars -except X x_opt_dmd ts nheight nwidth fs

nx = size(x_opt_dmd,2);
X = X(:,1:nx);

% you can run the script 'expvideo_dmd' to obtain the response x_dmd
% estimated by DMD, if you have access to the video '00_00_00_run1.mp4' 
% of the beam by Garrido et al. % (https://doi.org/10.1016/j.ymssp.2023.110539)

%% Video exp_beam_motion

% vm = VideoWriter('exp_beam_motion.mp4','MPEG-4');
% vm.FrameRate = 30;

a_opt_dmd = double(min(real(x_opt_dmd(:))));
b_opt_dmd = double(max(real(x_opt_dmd(:))-a_opt_dmd));

a_dmd = double(min(real(x_dmd(:))));
b_dmd = double(max(real(x_dmd(:))-a_dmd));

figure, 

% subplot(3,1,1),
% imagesc(reshape(abs(X(:,1)),[nheight nwidth]))
% 
% axis equal
% ax = gca;
% ax.NextPlot = 'replaceChildren';
% ax.YAxis.Visible = 'off'; % remove y-axis
% ax.XAxis.Visible = 'off'; % remove x-axis
% title(['Original Video t = ', num2str(1/fs,'%.2f'),' s'])
% colormap(gray)
% c = colorbar;
% clim('manual')
% clim([0 255])
% c.Ticks = linspace(50,250,5);
% c.TickLabels = num2cell(linspace(50,250,5));

subplot(2,1,1),
xIdmd = reshape(abs(x_opt_dmd(:,1)),[nheight nwidth]);
xIdmd = 255*(xIdmd - a_opt_dmd)./(b_opt_dmd - a_opt_dmd);
imagesc(xIdmd)

axis equal
ax = gca;
ax.NextPlot = 'replaceChildren';
ax.YAxis.Visible = 'off'; % remove y-axis
ax.XAxis.Visible = 'off'; % remove x-axis
title(['OptDMD t = ', num2str(1/fs,'%.2f'),' s'])
colormap(ax,gray(25))
colorbar
clim('manual')

subplot(2,1,2),
xIdmd = reshape(abs(x_dmd(:,1)),[nheight nwidth]);
xIdmd = 255*(xIdmd - a_dmd)./(b_dmd - a_dmd);
imagesc(xIdmd)

axis equal
ax = gca;
ax.NextPlot = 'replaceChildren';
ax.YAxis.Visible = 'off'; % remove y-axis
ax.XAxis.Visible = 'off'; % remove x-axis
title('DMD')
colormap(ax,gray(25))
colorbar
clim('manual')

% open(vm)
for i = 1:500
    % subplot(3,1,1),
    % imagesc(reshape(abs(X(:,i)),[nheight nwidth]))
    % title(['Original Video t = ', num2str(i/fs,'%.2f'),' s'])
    
    subplot(2,1,1),
    xIdmd = reshape(abs(x_opt_dmd(:,i)),[nheight nwidth]);
    xIdmd = 255*(xIdmd - a_opt_dmd)./(b_opt_dmd - a_opt_dmd);
    imagesc(xIdmd)
    title(['OptDMD t = ', num2str(i/fs,'%.2f'),' s'])

    subplot(2,1,2),
    xIdmd = reshape(abs(x_dmd(:,i)),[nheight nwidth]);
    xIdmd = 255*(xIdmd - a_dmd)./(b_dmd - a_dmd);
    imagesc(xIdmd)    
    
    drawnow

    % cdata = print('-RGBImage','-r300','-noui'); 
    % thisFrame = im2frame(cdata);  
    % writeVideo(vm, thisFrame); 

end
% close(vm)

%% Reconstruction of certain points

for i = 1:nx
    xrr(:,:,i) = reshape(X(:,i),[nheight nwidth]);
    
    xI_opt_dmd = reshape(abs(x_opt_dmd(:,i)),[nheight nwidth]);
    xdmdd(:,:,i) = 255*(xI_opt_dmd - a_opt_dmd)./(b_opt_dmd - a_opt_dmd);
    
    xIdmd = reshape(abs(double(x_dmd(:,i))),[nheight nwidth]);
    xdmddk(:,:,i) = 255*(xIdmd - a_dmd)./(b_dmd - a_dmd);
end

ydes = [105 105 50];
xdes = [500 50 20]; % x and y coordinates of the points selected for testing

M = xrr(:,:,1);

% % % Fig 12
figure, 
subplot(2,2,1), 
imshow(mat2gray(M)), colormap(bone)
hold on, 
scatter(xdes(1),ydes(1),'MarkerFaceColor',[0.74 0.16 0.17],'MarkerEdgeColor',[0.74 0.16 0.17])
scatter(xdes(2),ydes(2),'MarkerFaceColor',[0.91 0.78 0.086],'MarkerEdgeColor',[0.91 0.78 0.086])
scatter(xdes(3),ydes(3),'MarkerFaceColor',[0.31 0.68 0.622],'MarkerEdgeColor',[0.31 0.68 0.622])
title('Position of the pixels')

for j = 1:3
    xr = zeros(length(ts),1);
    xd = zeros(length(ts),1);
    xdk = zeros(length(ts),1);
    for i = 1:ts(end)*fs+1
        xr(i) = xrr(ydes(j),xdes(j),i);
        xd(i) = xdmdd(ydes(j),xdes(j),i);
        xdk(i) = xdmddk(ydes(j),xdes(j),i);
    end
    
    subplot(2,2,j+1), hold on,
    plot(ts, xr/255,'LineWidth',0.4,'Color',[0.4 0.4 0.4])
    plot(ts, xd/255,'LineWidth',1,'Color',[0 0 0.64])
    plot(ts, xdk/255,'LineWidth',1,'Color',[0 0.64 0])
    grid on
    xlim([ts(1) ts(end)])
    ylim([0 1])
    xlabel('Time (s)')
    ylabel('Pixel intensity')
    switch j
        case 1
            title('Pixel A','Color',[0.74 0.16 0.17])
            legend('original','optDMD','DMD',Location='southeast')
        case 2
            title('Pixel B','Color',[0.91 0.78 0.086])
        case 3
            title('Pixel C','Color',[0.31 0.68 0.622])
    end

end

clear xrr xdmdd xdmddk M

%% Cosine Similarity

for i = 1:nx
    XI = X(:,i);

    xI_opt_dmd = abs(x_opt_dmd(:,i));
    xI_opt_dmd = 255*(xI_opt_dmd - a_opt_dmd)./(b_opt_dmd - a_opt_dmd);

    xIdmd = abs(double(x_dmd(:,i)));
    xIdmd = 255*(xIdmd - a_dmd)./(b_dmd - a_dmd);

    d_opt_dmd(i) = 1 - pdist2(XI',xI_opt_dmd',"cosine");
    d_dmd(i) = 1 - pdist2(XI',xIdmd',"cosine");
end

% % % Fig 13
figure, plot([ts(1) ts(end)],100*[1 1],'--','LineWidth',1,'Color',[0.65 0.65 0.65])
hold on
plot(ts,100*d_opt_dmd,'LineWidth',1,'Color',[0 0 0.64])
plot(ts,100*d_dmd,'LineWidth',1,'Color',[0 0.64 0])
grid on
legend('reference','OptDMD','DMD',Location='southwest')
title('Cosine Similarity')
xlabel('Time (s)')
ylabel('Similarity %')
ylim(100*[0.97 1.001])

clear XI xI_opt_dmd xIdmd