clc
clear all
close all

addpath('..\Dados\plastic beam - high speed camera\videos 2022-04-05')
% video records included in the dataset by Garrido et al.
% (https://doi.org/10.1016/j.ymssp.2023.110539)

addpath('../Functions')

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
% nV is selected by inspection for the videos in the dataset
nV = 391:600; %00_00_00_run1.mp4

roif = vframes(nV,:,:);

clearvars -except roif

%% Collecting data for the snaphot matrix

fs = 1000; % frame rate in Hz
T = 0.6; % analysis of only the first 0.6 s of the recorded video

xy = roif(:,:,1:1:T*fs+1);

clear roif

nheight = size(xy,1);
nwidth = size(xy,2);
nframes = size(xy,3);

npixels = nheight*nwidth;

X = zeros(npixels,nframes,'single');
for i = 1:nframes
    X(:,i) = reshape(xy(:,:,i),[npixels,1]);
end

clear xy npixels nframes nheight nwidth

%% Applying DMD with time embeddings for a model order r_test

na = 7; % number of time embeddings considered (m = na - 1)
r_test = 200;

rb = r_test;
[~, ~, ~, s_dmds0, ~] = dmd_free_time_response(X, fs, T, na, rb, T);

rb = r_test + 1;
[~, ~, ~, s_dmds1, ~] = dmd_free_time_response(X, fs, T, na, rb, T);

%% Selecting the frequency-stable eigenvalues estimated by DMD at r_test as initial guesses for optDMD

kk = [r_test r_test+1];

sr = [ [s_dmds0; 0] s_dmds1];

nr = length(kk);

fr = abs(sr)/(2*pi);
dr = 0*fr;
for i = 1:nr
    dr(:,i) = -real(sr(:,i))./abs(sr(:,i));
end

fst = false(size(fr,1),1);
dst = false(size(fr,1),1);

tol = [0.005 0.005];

i = 1;

f0 = fr(1:kk(i), i);
f1 = fr(1:kk(i+1), i+1);
d0 = dr(1:kk(i), i);
d1 = dr(1:kk(i+1), i+1);    
for j = 1:kk(i)
    if abs(f0(j)) < tol(1) && min(abs(f0(j)-f1(:))) < tol(1)
        fst(j) = true;
    end
    if min(abs(f0(j)-f1(:))) < tol(1)*f0(j)
        fst(j) = true;
    end
end

fj = fr(1:r_test,i);
fst = fst(1:r_test);

A = fj.*fst;
A = nonzeros(A);

B = fj; 

[~, idx] = ismembertol(A,B,1e-5);

idxx = zeros(size(idx));
flag_prox = 0;
flag_ant = 0;
for i = 1:length(idx)-1
    if idx(i) == idx(i+1)
        idxx(i) = idx(i);
        idxx(i+1) = idx(i) + 1;
        flag_prox = 1;
    else
        flag_ant = flag_prox;
        flag_prox = 0;
    end

    if ~flag_ant && ~flag_prox
        idxx(i) = idx(i);
    end

    if i == length(idx)-1 && idx(i) ~= idx(i+1)
        idxx(i+1) = idx(i+1);
    end
end
%% Saving initial guess to be used in optDMD 

e_opt = double(s_dmds0);

nk = sort(idxx);
e_opt_0 = e_opt(nk);

f_dmd = abs(e_opt_0)/(2*pi);
zeta_dmd = -real(e_opt_0)./abs(e_opt_0);
dmd = [f_dmd zeta_dmd];

% save('../Data/initial_guess4optdmd_data_00_00_00_run1I.mat','e_opt_0')
