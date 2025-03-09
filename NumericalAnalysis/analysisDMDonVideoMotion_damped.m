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
% epsilon = 0.04; % damping factor

load displ_beam_impulse025N_damp4p.mat
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

[dmds, Phis, x_modo,~,~] = dmd_free_time_response(X, fs, T, na, rb, T);

%% Comparing some of the eigenvalues estimated by DMD with the analytical model
n_order = [9 22 73 108]; % indexes determined by inspection of the temporal 
% patterns that correspond to the natural frequencies and damping factor 
% of the analytical beam

f_dmd = dmds(n_order,1);
zeta_dmd = dmds(n_order,2);
for i = 1:length(n_order)
    e_f2(i) = abs(f_dmd(i) - f_exact(i))./f_exact(i);
    e_zeta2(i) = abs(zeta_dmd(i) - zeta_exact(i))/zeta_exact(i);
end

% % % Table 3
ff = [f_exact(1:length(n_order)) f_dmd 100*e_f2'];
zetazeta = [zeta_exact(1:length(n_order)) zeta_dmd 100*e_zeta2'];

%% By stacking past measurements to enrich the dataset, Phis must be restricted to now only consider the first frame of each snapshot

npixels = nheight*nwidth;
Phis = Phis(1:npixels,:);
x_dmd = Phis*x_modo; % data reconstructed by the DMD model

ts = 0:1/fs:(size(x_modo,2)-1)/fs;

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
xdes = [350 500 520];% x and y coordinates of the points selected for testing

% % % Fig 9

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

%% Video Beam_damped_motion
% vm = VideoWriter('beam_dampedmotion.mp4','MPEG-4');
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
colormap(ax1, gray(3))
title('Reconstructed Video')
clim("manual") 
% clim([0 1])
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
