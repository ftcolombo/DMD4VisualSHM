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
% nV is selected by inspection for all videos in the dataset
nV = 381:610;

%% Corrupting the video, with the goal of analysing DMD with only low-quality images

t_samp = 2; % downsampling the dataset by a factor of 2
frame_orig = vframes(nV,:,1:t_samp:end);

clearvars -except frame_orig t_samp

% rng(...) % the Gaussian noise added depends on the seed of the generator 
% the seed rng(115123) gives similar results

r_size = 0.9; % reducing the resolution to 90% and then transforming back to 100% to create artifacts
var_gauss = 0.005; % adding Gaussian white noise with zero mean and variance of 0.005 

for i = 1:size(frame_orig,3)
    frame_aux = imresize(frame_orig(:,:,i), r_size);
    frame_aux = imresize(frame_aux, 1/r_size);
    frame_noise(:,:,i) = imnoise(frame_aux,'gaussian',0,var_gauss);
end

vframes = frame_noise;

clearvars -except vframes t_samp

%% Collecting data for the snaphot matrix

fs = 1000/t_samp; % 1000 Hz was the original frame rate before downsampling

Ti = 0;
Tf = 2.5;
ts_train = Ti:1/fs:Tf;
frames = vframes(:,:,Ti*fs+1:1:Tf*fs+1);

roif = double(frames);

nheight = size(roif,1);
nwidth = size(roif,2);

xy = roif;

nframes = size(xy,3);
npixels = nheight*nwidth;

X_train = zeros(npixels,nframes);
for i = 1:nframes
    X_train(:,i) = reshape(xy(:,:,i),[npixels,1]);
end

clear xy roif frames vframes 

%% Applying optDMD with e_opt_0 as initial guesses for the eigenvalues

% you can run the script 'initial_modeselection4diagnostics' to obtain the 
% most dominant and frequency-stable eigenvalues estimated by DMD for the 
% first health scenario, if you have access to the % video '00_00_00_run1.mp4' of 
% the beam by Garrido et al. % (https://doi.org/10.1016/j.ymssp.2023.110539)
load 'initial_guess4optdmd_00_00_00_run1_lowqualityf500.mat'

% For subsequent scenarios, use the set of eigenvalues estimated by optdmd
% in the previous scenario

addpath('..\optdmd-master\src')
% to run the comand optdmd, the repository https://github.com/duqbo/optdmd
% is required (https://doi.org/10.1137/M1124176)

[w_opt,e_opt2,b_opt] = optdmd(X_train,ts_train,length(e_opt_0),2,varpro_opts('ifprint',0),e_opt_0);


%% Comparing the natural frequencies and decay rates computed by DMD with time embeddings and optDMD

f_optdmd = abs(e_opt2)/(2*pi);
zeta_optdmd = -real(e_opt2)./abs(e_opt2);
dmd_opt = [f_optdmd zeta_optdmd];

f_dmd = abs(e_opt_0)/(2*pi);
zeta_dmd = -real(e_opt_0)./abs(e_opt_0);
dmd = [f_dmd zeta_dmd];

%% Applying the framework proposed by Sashidhar and Kutz to obtain DMD with uncertainty quantification
% available in the repository https://github.com/dsashid/BOP-DMD
% more information on the paper https://doi.org/10.1098/rsta.2021.0199

num_modes = length(e_opt2);

p = round(0.25*size(X_train,2)); % using 25% of the training data for each model
num_cycles =  50; % number of models to be estimated

for j = 1:num_cycles
    unsorted_ind = randperm(nframes,p); %randomly select p indices from nframes time points
    ind = sort(unsorted_ind);
    xdata_cycle = X_train(:,ind);
    ts_ind = ts_train(ind);
    
    [w_cycle,e1_cycle,b_cycle] = optdmd(xdata_cycle,ts_ind,num_modes,2,varpro_opts('ifprint',0),e_opt2);
    
    [~,imag_ind] = sort(imag(e1_cycle));
    
    b_vec_ensembleDMD(:,j) = b_cycle(imag_ind);
    lambda_vec_ensembleDMD(:,j) = e1_cycle(imag_ind);
    w_vec_ensembleDMD(:,(j-1)*num_modes+1:j*num_modes) = w_cycle(:,imag_ind);
   
end

clear xdata_cycle w_cycle e1_cycle b_cycle ts_ind ind unsorted_ind

%% Computing the histogram of the eigenvalues of all models estimated

sortedLambda_ensembleDMD = sort(lambda_vec_ensembleDMD,1,'ComparisonMethod','abs');
abs_sorted = abs(sortedLambda_ensembleDMD)/(2*pi);

nbins = [];

% % % Figure 16
figure,
j = 1;
for i = 1:2:size(sortedLambda_ensembleDMD,1)
    subplot(4,2,2*j-1)
    h = histfit(abs(sortedLambda_ensembleDMD(i,:)/(2*pi)),nbins);
    h(1).FaceColor = [0 0 0.55];
    h(2).Color = [0.74 0.16 0.17];
    
    mean_value_ensemble(j,1) = mean(abs(sortedLambda_ensembleDMD(i,:))/(2*pi));
    ax = gca;
    hold on
    plot([mean_value_ensemble(j,1) mean_value_ensemble(j,1)], ax.YLim,'--','color',[0.91 0.78 0.086] ,'linewidth',1);
    title([num2str(mean_value_ensemble(j,1),'%.2f'), ' Hz'])
    
    subplot(4,2,2*j)
    h = histfit(-real(sortedLambda_ensembleDMD(i,:))./abs(sortedLambda_ensembleDMD(i,:)),nbins);
    h(1).FaceColor = [0 0 0.55];
    h(2).Color = [0.74 0.16 0.17];
    
    mean_value_ensemble(j,1) = mean(-real(sortedLambda_ensembleDMD(i,:))./abs(sortedLambda_ensembleDMD(i,:)));
    ax = gca;
    hold on
    plot([mean_value_ensemble(j,1) mean_value_ensemble(j,1)], ax.YLim,'--','color',[0.91 0.78 0.086] ,'linewidth',1);
    title([num2str(100*mean_value_ensemble(j,1),'%.2f'), ' %'])
    j = j+1;
end


%% Computing the distribution of the eigenvalues in the s-plane of all models estimated

for jj = 1:num_modes
    
    pd_real = fitdist(real(lambda_vec_ensembleDMD(jj,:)'),'Normal');
    pdf_lambda_vec_real(1,jj) = pd_real.mu;
    pdf_lambda_vec_real(2,jj) = pd_real.sigma;
    
    pd_imag = fitdist(abs(imag(lambda_vec_ensembleDMD(jj,:)')),'Normal');
    pdf_lambda_vec_imag(1,jj) = pd_imag.mu;
    pdf_lambda_vec_imag(2,jj) = pd_imag.sigma;
end

theta = linspace(0,2*pi,1000);


% % % Figure 17 (only scenario A)
figure, hold on
for i = 1:num_modes

    x_eli = pdf_lambda_vec_real(1,i) + 2*pdf_lambda_vec_real(2,i)*sin(theta);
    y_eli = pdf_lambda_vec_imag(1,i) + 2*pdf_lambda_vec_imag(2,i)*cos(theta);
    patch('YData', y_eli,'XData',x_eli, 'LineStyle','none','FaceColor',[0.7 0.7 0.7])

    x_eli = pdf_lambda_vec_real(1,i) + pdf_lambda_vec_real(2,i)*sin(theta);
    y_eli = pdf_lambda_vec_imag(1,i) + pdf_lambda_vec_imag(2,i)*cos(theta);
    patch('YData', y_eli,'XData',x_eli, 'LineStyle','none','FaceColor',[0.2 0.3 0.5])

    scatter(pdf_lambda_vec_real(1,i),pdf_lambda_vec_imag(1,i),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor','none')    
end
grid on
xlabel('Real axis')
ylabel('Imaginary axis')
title('Poles of scenario A')


%% Computing the mean eigenvectors, eigenvalues and b-vector

mean_modes = zeros(npixels,num_modes);

for i = 1:ceil(num_modes/2)
    w_sample = w_vec_ensembleDMD(:,i:num_modes:end);
    mean_modes(:,i) = mean(w_sample,2);
end

mean_lambda = mean(lambda_vec_ensembleDMD,2);
mean_b = mean(b_vec_ensembleDMD,2);

%% Computing the Image Dinamic Mode to be used as a damage-sensitive feature

np = [4 3 2 1]; % ordering of the modes according to the mean eigs
modes1 = mean_modes(:,np);

IPhi01 = mat2gray(reshape(real(modes1(:,1)), nheight,nwidth));
IPhi11 = mat2gray(reshape(real(modes1(:,2)), nheight,nwidth));
IPhi21 = mat2gray(reshape(real(modes1(:,3)), nheight,nwidth));
IPhi31 = mat2gray(reshape(real(modes1(:,4)), nheight,nwidth));

A = IPhi11;
B = IPhi01;

% % % (Image dinamic mode of scenario A)
figure, imshow(imabsdiff(A, B))