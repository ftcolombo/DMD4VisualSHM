close all
clear all
clc

addpath('../Data')
addpath('../Functions')

load Jcolormap.mat
load CCcool.mat

% properties of the clamped-free beam made of Polypropylene: 
% lbeam = 0.8; % m, length of the beam
% wbeam = 25.4e-3; % m, width/thickness of the beam
% E = 1300e6; % Pa, Young's modulus
% density = 900; % kg/m^3, density 
% epsilon = 0.00; % damping factor

load displ_beam_impulse025N_damp0p.mat
% y: displacement matrix over time due to an impulse force F = 0.25 % N,
% considering nelem = 50; elements along the beam
% t: time vector
% phi: matrix with the first 5 mode shapes
% wn_Hz: first 5 natural frequencies in Hz

fs = 1/(t(2)-t(1)); % sampling rate used in the simulation of the beam motion

f_exact = wn_Hz';
zeta_exact = epsilon*ones(size(f_exact));
nm = length(f_exact);

%% Collecting data for the snaphot matrix

img_noise = 0;  % flag for adding noise ( 1: noise is added, 0: no noise )
T = 1; % simulation of the motion during 1 s

captured_frames = makevideo(y, 0:1/fs:T, lbeam, wbeam, img_noise);

nheight = size(captured_frames,1);
nwidth = size(captured_frames,2);
nframes = size(captured_frames,3);

npixels = nheight*nwidth;

X = zeros(npixels,nframes,'single');
for i = 1:nframes
    X(:,i) = reshape(captured_frames(:,:,i),[npixels,1]);
end

%% Applying DMD with time embeddings

rb = 400; % maximum number of singular values to be kept in DMD
na = 10; % number of time embeddings considered (m = na - 1)

[dmds, Phis, x_modo, ~, ~] = dmd_free_time_response(X, fs, T, na, rb, T);

%% Comparing some of the eigenvalues estimated by DMD with the analytical model
n_order = [2 8 12 18 25]; % indexes determined by inspection of the temporal 
% patterns that correspond to the natural frequencies and damping factor 
% of the analytical beam

f_dmd = dmds(n_order,1);
zeta_dmd = dmds(n_order,2);
for i = 1:nm
    e_f2(i) = abs(f_dmd(i) - f_exact(i))./f_exact(i);
end

% % % Table 2
ff = [f_exact f_dmd 100*e_f2'];
zetazeta = [zeta_exact zeta_dmd];

%% By stacking past measurements to enrich the dataset, Phis must be restricted to now only consider the first frame of each snapshot

npixels = nheight*nwidth;
Phis = Phis(1:npixels,:);
x_dmd = Phis*x_modo; % data reconstructed by the DMD model

ts = 0:1/fs:(size(x_modo,2)-1)/fs;

%% Comparison between original and reconstructed frames at time k = i
i = 164;

% % % Fig 4
figure, 
sph = subplot(2,1,1); 
imshow(mat2gray(captured_frames(:,:,i)),'Colormap', gray)
colormap(gray)
caxis(sph,[0,1]);
title(['(a) Original frame at t = ' num2str(i/fs,'%.2f')])

sph = subplot(2,1,2); 
xIdmd = reshape(real(x_dmd(:,i)),[nheight nwidth]);
imshow(mat2gray(xIdmd),'Colormap', gray)
caxis(sph,[0,1]);
title(['(b) Reconstructed frame at t = ' num2str(i/fs,'%.2f')])
ax = gca;
c = colorbar(ax,'Position',[0.92 0.12 0.025 0.8]);  % attach colorbar to h
colormap(c,gray)
caxis(ax,[0,1]); 

%% Evaluating the spatial patterns estimated by DMD

nk = 1:21; % considering only the first 21 patterns

IPhi = zeros(nheight,nwidth,8,'single');
j = 1;
for i = 1:length(nk)
    IPhi(:,:,j) = mat2gray(reshape(real(Phis(:,i)),[nheight nwidth]));
    j = j+1;
end
minColorLimit = min(min(min(IPhi))); 
maxColorLimit = max(max(max(IPhi)));

nk = [2 8 12 18]; % indexes determined by inspection that correspond to Modes 1, 2, 3 and 4

% % % Fig 5
figure,
for i = 1:length(nk)
    subplot(length(nk),1,i), imshow(IPhi(:,:,nk(i)),'Colormap', CC)
    clim([minColorLimit,maxColorLimit])
end

nk = [4 6 10];  % indexes determined by inspection that correspond to Modes 1.1, 1.2 and 1.3

% % % Fig 6
figure,
for i = 1:length(nk)
    subplot(length(nk),1,i), imshow(IPhi(:,:,nk(i)),'Colormap', CC)
    clim([minColorLimit,maxColorLimit])
end

nk = 1;  % index determined by inspection that correspond to Mode 0

% % % Fig 7
figure, imshow(IPhi(:,:,nk),'Colormap', CC)
clim([minColorLimit,maxColorLimit])

% colors of the spatial patterns may change according to the normalization
% of the computed eigenvectors

%% Video Beam_4Modes

% indexes determined by inspection:
k1 = [1:7 10 11 14:17 20 21 23 34]; % Mode 1 + harmonics + Mode 0
k2 = [1 8 9 36 37]; % Mode 2 + harmonics + Mode 0
k3 = [1 12 13 89 90]; % Mode 3 + harmonics + Mode 0
k4 = [1 18 19]; % Mode 4 + Mode 0
k5 = [1 25 26]; % Mode 5 + Mode 0

x_dmdk1 = Phis(:,k1)*x_modo(k1,:);
x_dmdk2 = Phis(:,k2)*x_modo(k2,:);
x_dmdk3 = Phis(:,k3)*x_modo(k3,:);
x_dmdk4 = Phis(:,k4)*x_modo(k4,:);
x_dmdk5 = Phis(:,k5)*x_modo(k5,:);


% vm = VideoWriter('beam_4_modes.mp4','MPEG-4');
% vm.FrameRate = 30;

xIdmd1 = reshape(real(x_dmdk1(:,1)),[nheight nwidth]);
xIdmd2 = reshape(real(x_dmdk2(:,1)),[nheight nwidth]);
xIdmd3 = reshape(real(x_dmdk3(:,1)),[nheight nwidth]);
xIdmd4 = reshape(real(x_dmdk4(:,1)),[nheight nwidth]);

a_dmd1 = double(min(real(x_dmdk1(:))));
b_dmd1 = double(max(real(x_dmdk1(:))-a_dmd1));
a_dmd2 = double(min(real(x_dmdk2(:))));
b_dmd2 = double(max(real(x_dmdk2(:))-a_dmd2));
a_dmd3 = double(min(real(x_dmdk3(:))));
b_dmd3 = double(max(real(x_dmdk3(:))-a_dmd3));
a_dmd4 = double(min(real(x_dmdk4(:))));
b_dmd4 = double(max(real(x_dmdk4(:))-a_dmd4));

a_dmd = min([a_dmd1 a_dmd2 a_dmd3 a_dmd4]);
b_dmd = max([b_dmd1 b_dmd2 b_dmd3 b_dmd4]);

a_dmd1 = a_dmd;
b_dmd1 = b_dmd;
a_dmd2 = a_dmd;
b_dmd2 = b_dmd;
a_dmd3 = a_dmd;
b_dmd3 = b_dmd;
a_dmd4 = a_dmd;
b_dmd4 = b_dmd;

amax = 1;%255;

figure, 
subplot(2,2,1),
xIdmd = amax*(xIdmd1 - a_dmd1)./(b_dmd1 - a_dmd1);
imagesc(xIdmd)
colormap(J)
axis equal
ax1 = gca;
ax1.NextPlot = 'replaceChildren';
ax1.YAxis.Visible = 'off'; % remove y-axis
ax1.XAxis.Visible = 'off'; % remove x-axis
title('Mode 1')
clim("manual")         % allow subsequent plots to use the same color limits

subplot(2,2,2),
xIdmd = amax*(xIdmd2 - a_dmd2)./(b_dmd2 - a_dmd2);
imagesc(xIdmd)
colormap(J)
axis equal
ax1 = gca;
ax1.NextPlot = 'replaceChildren';
ax1.YAxis.Visible = 'off'; 
ax1.XAxis.Visible = 'off'; 
title('Mode 2')
clim("manual")

subplot(2,2,3),
xIdmd = amax*(xIdmd3 - a_dmd3)./(b_dmd3 - a_dmd3);
imagesc(xIdmd)
colormap(J)
axis equal
ax1 = gca;
ax1.NextPlot = 'replaceChildren';
ax1.YAxis.Visible = 'off';
ax1.XAxis.Visible = 'off';
title('Mode 3')
clim("manual")

subplot(2,2,4),
xIdmd = amax*(xIdmd4 - a_dmd4)./(b_dmd4 - a_dmd4);
imagesc(xIdmd)
colormap(J)
axis equal
ax1 = gca;
ax1.NextPlot = 'replaceChildren';
ax1.YAxis.Visible = 'off';
ax1.XAxis.Visible = 'off';
title('Mode 4')
clim("manual")
colorbar('south')

% open(vm)
for i = 1:250
    xIdmd1 = reshape(real(x_dmdk1(:,i)),[nheight nwidth]);
    xIdmd2 = reshape(real(x_dmdk2(:,i)),[nheight nwidth]);
    xIdmd3 = reshape(real(x_dmdk3(:,i)),[nheight nwidth]);
    xIdmd4 = reshape(real(x_dmdk4(:,i)),[nheight nwidth]);
    
    subplot(2,2,1),
    xIdmd = amax*(xIdmd1 - a_dmd1)./(b_dmd1 - a_dmd1);
    imagesc(xIdmd)
    subplot(2,2,2),
    xIdmd = amax*(xIdmd2 - a_dmd2)./(b_dmd2 - a_dmd2);
    imagesc(xIdmd)
    subplot(2,2,3),
    xIdmd = amax*(xIdmd3 - a_dmd3)./(b_dmd3 - a_dmd3);
    imagesc(xIdmd)    
    subplot(2,2,4),
    xIdmd = amax*(xIdmd4 - a_dmd4)./(b_dmd4 - a_dmd4);
    imagesc(xIdmd)
    drawnow

    % cdata = print('-RGBImage','-r300','-noui'); 
    % thisFrame = im2frame(cdata);  
    % writeVideo(vm, thisFrame); 
    
    pause(0.05)
end
% close(vm)

clear x_dmdk1 x_dmdk2 x_dmdk3 x_dmdk4 x_dmdk5 xIdmd1 xIdmd2 xIdmd3 xIdmd4

%% Reconstruction of certain points

x = 0:lbeam/size(y,2):lbeam;

f1 = figure;
plot(x, [0 y(1,:)], 'LineWidth',14, 'Color','k')
axis equal
ax_vert = 2*max(abs(y(:)));
axis([x(1) x(end) -ax_vert ax_vert])
ax = gca;
ax.Visible = 'off';
M = getframe;
close(f1)

a_dmd = double(min(real(x_dmd(:))));
b_dmd = double(max(real(x_dmd(:))-a_dmd));

for i = 1:size(x_modo,2)
    xrr(:,:,i) = captured_frames(:,:,i);
    
    xIdmd = reshape(real(x_dmd(:,i)),[nheight nwidth]);
    xdmdd(:,:,i) = (xIdmd - a_dmd)./(b_dmd- a_dmd);
end

ydes = [ceil(nheight/2) ceil(nheight/2) ceil(nheight/2)-30];
xdes = [350 500 520]; % x and y coordinates of the points selected for testing

% % % Fig 8
figure, 
subplot(2,2,1), imshow(M.cdata)
hold on, 
scatter(xdes(1),ydes(1),'MarkerFaceColor',[0.74 0.16 0.17],'MarkerEdgeColor',[0.74 0.16 0.17])
scatter(xdes(2),ydes(2),'MarkerFaceColor',[0.91 0.78 0.086],'MarkerEdgeColor',[0.91 0.78 0.086])
scatter(xdes(3),ydes(3),'MarkerFaceColor',[0.31 0.68 0.622],'MarkerEdgeColor',[0.31 0.68 0.622])
title('Position of the pixels')

for j = 1:3
    xr = zeros(length(ts),1);
    xd = zeros(length(ts),1);
    for i = 1:ts(end)*fs+1
        xr(i) = xrr(ydes(j),xdes(j),i);
        xd(i) = xdmdd(ydes(j),xdes(j),i);
    end
    
    subplot(2,2,j+1), hold on,
    plot(ts, xr/255,'Color',[0.4 0.4 0.4])
    plot(ts, xd,'LineWidth',1,'Color',[0 0 0.64])
    grid on
    xlim([0 0.6])
    xlabel('Time (s)')
    ylabel('Pixel intensity')
    switch j
        case 1
            title('Pixel A','Color',[0.74 0.16 0.17])
        case 2
            title('Pixel B','Color',[0.91 0.78 0.086])
        case 3
            title('Pixel C','Color',[0.31 0.68 0.622])
            legend('original','reconstructed','Location','best')
    end

end

%% Video Beam_undamped_motion

% vm = VideoWriter('beam_undamped_motion.mp4','MPEG-4');
% vm.FrameRate = 30;

xIdmd = reshape(real(x_dmd(:,1)),[nheight nwidth]);

a_dmd = double(min(real(x_dmd(:))));
b_dmd = double(max(real(x_dmd(:))-a_dmd));

amax = 1;%255;

figure, 
subplot(2,1,1),
imagesc(captured_frames(:,:,1)/255)
axis equal
ax1 = gca;
ax1.NextPlot = 'replaceChildren';
ax1.YAxis.Visible = 'off'; % remove y-axis
ax1.XAxis.Visible = 'off'; % remove x-axis
colormap(ax1, gray)
title(['Original Video t = ', num2str(1/fs,'%.2f'),' s'])
clim("manual")         % allow subsequent plots to use the same color limits
clim([0 1])
colorbar

subplot(2,1,2),
xIdmdI = amax*(xIdmd - a_dmd)./(b_dmd - a_dmd);
imagesc(xIdmdI)
axis equal
ax1 = gca;
ax1.NextPlot = 'replaceChildren';
ax1.YAxis.Visible = 'off';
ax1.XAxis.Visible = 'off';
title('Reconstructed Video')
colormap(ax1, gray(3))
clim("manual")
colorbar

% open(vm)
for i = 1:500

    
    subplot(2,1,1),
    imagesc(captured_frames(:,:,i)/255)
    title(['Original Video t = ', num2str(i/fs,'%.2f'),' s'])

    subplot(2,1,2),
    xIdmd = reshape(real(x_dmd(:,i)),[nheight nwidth]);
    xIdmdI = amax*(xIdmd - a_dmd)./(b_dmd - a_dmd);
    imagesc(xIdmdI)
    drawnow

    % cdata = print('-RGBImage','-r300','-noui');
    % thisFrame = im2frame(cdata);
    % writeVideo(vm, thisFrame); 
    
    pause(0.05)
end
% close(vm)