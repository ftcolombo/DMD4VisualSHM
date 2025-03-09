clc
clear all
close all

addpath('../Data')

%% Analysing the mean eigenvectors for all damage scenarios

load('mean_mode_00_00_00_run1_noise.mat')
modesA = mean_modes;
nhA = nheight;
nwA = nwidth;

load('mean_mode_05_00_00_run1_noise.mat')
modesB = mean_modes;
nhB = nheight;
nwB = nwidth;

load('mean_mode_10_00_00_run1_noise.mat')
modesC = mean_modes;
nhC = nheight;
nwC = nwidth;

load('mean_mode_13_00_00_run1_noise.mat')
modesD = mean_modes;
nhD = nheight;
nwD = nwidth;

load('mean_mode_13_05_00_run2_noise.mat')
modesE = mean_modes;
nhE = nheight;
nwE = nwidth;

load('mean_mode_13_05_05_run2_noise.mat')
modesF = mean_modes;
nhF = nheight;
nwF = nwidth;

load('mean_mode_13_10_05_run2_noise.mat')
modesG = mean_modes;
nhG = nheight;
nwG = nwidth;

load('mean_mode_13_10_11_run2_noise.mat')
modesH = mean_modes;
nhH = nheight;
nwH = nwidth;


%%
i = 1;

IPhi01 = mat2gray(reshape(real(modesA(:,i)), nhA,nwA));
IPhi02 = mat2gray(reshape(real(modesB(:,i)), nhB,nwB));
IPhi03 = mat2gray(reshape(real(modesC(:,i)), nhC,nwC));
IPhi04 = mat2gray(reshape(real(modesD(:,i)), nhD,nwD));
IPhi05 = mat2gray(reshape(real(modesE(:,i)), nhE,nwE));
IPhi06 = mat2gray(reshape(real(modesF(:,i)), nhF,nwF));
IPhi07 = mat2gray(reshape(real(modesG(:,i)), nhG,nwG));
IPhi08 = mat2gray(reshape(real(modesH(:,i)), nhH,nwH));

i = 2;

IPhi11 = mat2gray(reshape(real(modesA(:,i)), nhA,nwA));
IPhi12 = mat2gray(reshape(real(modesB(:,i)), nhB,nwB));
IPhi13 = mat2gray(reshape(real(modesC(:,i)), nhC,nwC));
IPhi14 = mat2gray(reshape(real(modesD(:,i)), nhD,nwD));
IPhi15 = mat2gray(reshape(real(modesE(:,i)), nhE,nwE));
IPhi16 = mat2gray(reshape(real(modesF(:,i)), nhF,nwF));
IPhi17 = mat2gray(reshape(real(modesG(:,i)), nhG,nwG));
IPhi18 = mat2gray(reshape(real(modesH(:,i)), nhH,nwH));

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
