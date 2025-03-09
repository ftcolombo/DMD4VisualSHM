clc
clear all
close all

addpath('../Data')

%% Analysing the mean eigenvectors for all damage scenarios

np = [4 3 2 1]; % ordering of the modes according to the eigs

load('mean_mode_00_00_00_run1_noise.mat')
modes1 = mean_modes(:,np);
nh1 = nheight;
nw1 = nwidth;

load('mean_mode_05_00_00_run1_noise.mat')
modes2 = mean_modes(:,np);
nh2 = nheight;
nw2 = nwidth;

load('mean_mode_10_00_00_run1_noise.mat')
modes3 = mean_modes(:,np);
nh3 = nheight;
nw3 = nwidth;

load('mean_mode_13_00_00_run1_noise.mat')
modes4 = mean_modes(:,np);
nh4 = nheight;
nw4 = nwidth;

load('mean_mode_13_05_00_run2_noise.mat')
modes5 = mean_modes(:,np);
nh5 = nheight;
nw5 = nwidth;

load('mean_mode_13_05_05_run2_noise.mat')
modes6 = mean_modes(:,np);
nh6 = nheight;
nw6 = nwidth;

load('mean_mode_13_10_05_run2_noise.mat')
modes7 = mean_modes(:,np);
nh7 = nheight;
nw7 = nwidth;

load('mean_mode_13_10_11_run2_noise.mat')
modes8 = mean_modes(:,np);
nh8 = nheight;
nw8 = nwidth;


%%
i = 1;

IPhi01 = mat2gray(reshape(real(modes1(:,i)), nh1,nw1));
IPhi02 = mat2gray(reshape(real(modes2(:,i)), nh2,nw2));
IPhi03 = mat2gray(reshape(real(modes3(:,i)), nh3,nw3));
IPhi04 = mat2gray(reshape(real(modes4(:,i)), nh4,nw4));
IPhi05 = mat2gray(reshape(real(modes5(:,i)), nh5,nw5));
IPhi06 = mat2gray(reshape(real(modes6(:,i)), nh6,nw6));
IPhi07 = mat2gray(reshape(real(modes7(:,i)), nh7,nw7));
IPhi08 = mat2gray(reshape(real(modes8(:,i)), nh8,nw8));

i = 2;

IPhi11 = mat2gray(reshape(real(modes1(:,i)), nh1,nw1));
IPhi12 = mat2gray(reshape(real(modes2(:,i)), nh2,nw2));
IPhi13 = mat2gray(reshape(real(modes3(:,i)), nh3,nw3));
IPhi14 = mat2gray(reshape(real(modes4(:,i)), nh4,nw4));
IPhi15 = mat2gray(reshape(real(modes5(:,i)), nh5,nw5));
IPhi16 = mat2gray(reshape(real(modes6(:,i)), nh6,nw6));
IPhi17 = mat2gray(reshape(real(modes7(:,i)), nh7,nw7));
IPhi18 = mat2gray(reshape(real(modes8(:,i)), nh8,nw8));

A = IPhi11;
B = IPhi01;
IPhi_dif1 = imabsdiff(A, B);

A = IPhi12;
B = IPhi02;
IPhi_dif2 = imabsdiff(imcomplement(A), B); 
% A or imcomplement(A) depends on the sign of the eigenvector 

A = IPhi13;
B = IPhi03;
IPhi_dif3 = imabsdiff(imcomplement(A), B);

A = IPhi14;
B = IPhi04;
IPhi_dif4 = imabsdiff(A, B);

A = IPhi15;
B = IPhi05;
IPhi_dif5 = imabsdiff(A, B);

A = IPhi16;
B = IPhi06;
IPhi_dif6 = imabsdiff(A, B);

A = IPhi17;
B = IPhi07;
IPhi_dif7 = imabsdiff(A, B);

A = IPhi18;
B = IPhi08;
IPhi_dif8 = imabsdiff(A, B);

clearvars -except IPhi_dif1 IPhi_dif2 IPhi_dif3 IPhi_dif4 IPhi_dif5 IPhi_dif6 IPhi_dif7 IPhi_dif8

%%
figure, imshow(IPhi_dif1)
figure, imshow(IPhi_dif2)
figure, imshow(IPhi_dif3)
figure, imshow(IPhi_dif4)
figure, imshow(IPhi_dif5)
figure, imshow(IPhi_dif6)
figure, imshow(IPhi_dif7)
figure, imshow(IPhi_dif8)

%%
close all
A = IPhi_dif1;

BW = imbinarize(A, 'adaptive', 'ForegroundPolarity', 'dark');
figure, imshow(BW)

[B, L] = bwboundaries(BW);
figure, imshow(label2rgb(L,@jet,[.5 .5 .5]))
hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1),'w', 'LineWidth',2)    
end

k = 2;
roi = L == k;
figure, imshow(roi)

boundary = B{k};
maxX = max(boundary(:,2));
indmaxX = find(boundary(:,2) == maxX);
yavg_maxX = mean(boundary(indmaxX,1));

traslX = -50;
ctr = [maxX + traslX yavg_maxX];

%%
close all
A = IPhi_dif8;

BW = imbinarize(A, 'adaptive', 'ForegroundPolarity', 'dark');
figure, imshow(BW)

[B, L] = bwboundaries(BW);
figure, imshow(label2rgb(L,@jet,[.5 .5 .5]))
hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1),'w', 'LineWidth',2)    
end

c_roi = round(ctr);

k = 1;
for i = c_roi(2)-10:1:c_roi(2)+10
    for j = c_roi(1)-10:1:c_roi(1)+10 
        mask(k) = L(i,j);
        k = k +1;
    end
end
M = mode(mask);

roi = L == M;
figure, imshow(roi)

%%
close all
A = IPhi_dif8;

BW = imbinarize(A, 'adaptive', 'ForegroundPolarity', 'dark');

% % % Figure 18
figure, 
subplot(4,1,1), imshow(A)
title('Image dinamic mode')

subplot(4,1,2), imshow(BW)
ax1 = gca;
colormap("gray")
title('Segmentation using adaptive thresholding')

[B, L] = bwboundaries(BW);
maxl = max(L(:));
L(L == 1) = max(L(:))+1;

subplot(4,1,3), imshow(label2rgb(L,'hsv','k','shuffle'))
title('Detection of all objects')

roi = L == 2;
subplot(4,1,4), imshow(roi)
ax1 = gca;
colormap("gray")
title('Identification of the beam')
