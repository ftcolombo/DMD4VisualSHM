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

roif = single(vframes(nV,:,:));

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

%% Applying DMD for an increasing number of time embeddings
rb = 200;  % maximum number of singular values to be kept in DMD

nai = 1:9; % increasing number of time embeddings (be careful with the amount of memory being used)
s_dmdk = NaN(rb,length(nai));

for i = 1:length(nai)
    i
    na = nai(i); % number of time embeddings considered (m = na - 1)

    [~, ~, ~, s_dmds, ~] = dmd_free_time_response(X, fs, T, na, rb, T);
    neigs = length(s_dmds);
    s_dmdk(1:neigs,i) = s_dmds;
end

% save('../Data/num_time_embed_data_00_00_00_run1I.mat','s_dmdk')
