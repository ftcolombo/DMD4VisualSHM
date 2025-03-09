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

% the seed rng(115123) gives similar results to the one in the paper (but
% not exactly, unfortunately)

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

T = 0.6; % analysis of only the first 0.6 s of the recorded video

xy = single(vframes(:,:,1:1:T*fs+1));

nheight = size(xy,1);
nwidth = size(xy,2);
nframes = size(xy,3);

npixels = nheight*nwidth;

X = zeros(npixels,nframes,'single');
for i = 1:nframes
    X(:,i) = reshape(xy(:,:,i),[npixels,1]);
end

clear xy npixels nframes

%% Applying DMD with time embeddings for a model order r_test

na = 6; % number of time embeddings considered (m = na - 1)
r_test = 200;

rb = r_test;
[~, ~, ~, s_dmds0, Ij] = dmd_free_time_response(X, fs, 1, na, rb, 0.6);

rb = r_test + 1;
[~, ~, ~, s_dmds1, ~] = dmd_free_time_response(X, fs, 1, na, rb, 0.6);

[~,ind] = sort(Ij, 'descend');
[~,indd] = sort(ind);

sr200 = s_dmds0(indd);
sr201 = s_dmds1;

%% Selecting the frequency-stable eigenvalues estimated by DMD at r_test as initial guesses for optDMD
% eigenvalues estimated vary depending on the random seed of the Gaussian noise

% load stabilizationdata_00_00_00_run1I_lowquality_f500 % if you want to see the exact result seen in the paper

kk = [r_test r_test+1];
sr = [ [sr200; 0] sr201];
nr = length(kk);

fr = abs(sr)/(2*pi);
dr = 0*fr;
for i = 1:nr
    dr(:,i) = -real(sr(:,i))./abs(sr(:,i));
end

tol = [0.005 0.005];

r_test = 200;

i = 1;

fj = fr(1:r_test,i);

fst = false(r_test,1);
dst = false(r_test,1);
finstt = false(r_test,1);

fn0 = fr(1:kk(i), i);
fn1 = fr(1:kk(i+1), i+1);
dn0 = dr(1:kk(i), i);
dn1 = dr(1:kk(i+1), i+1);    
for j = 1:kk(i)
    if abs(fn0(j)) < tol(1) && min(abs(fn0(j)-fn1(:))) < tol(1)
        if min(abs(dn0(j)-dn1(:))) < tol(2)*abs(dn0(j))
            dst(j) = true;
        else
            fst(j) = true;
        end
    elseif min(abs(fn0(j)-fn1(:))) < tol(1)*fn0(j)
        if min(abs(dn0(j)-dn1(:))) < tol(2)*dn0(j)
            dst(j) = true;
        else
            fst(j) = true;
        end
    else
        finstt(j) = true;
    end

end

Itotal = sum(Ij,"all");

Ij2 = Ij/Itotal;

Ijf = Ij2'.*fst;
Ijfd = Ij2'.*dst;
Ijinst = Ij2'.*finstt;

ff = fj.*fst;
ffd = fj.*dst;
finst = fj.*finstt;

thresh_Ij = 0.01;

% % % Fig 15
figure,
hold on
stem(finst, 100*Ijinst,'Marker','none','LineWidth',2,'Color',[0.65 0.65 0.65])
stem(ffd, 100*Ijfd,'Marker','none','LineWidth',2,'Color',[0.74 0.16 0.17])
stem(ff, 100*Ijf,'Marker','none','LineWidth',2,'Color', [0 0 0.55])
plot([0 50],[1 1]*100*thresh_Ij,'--','LineWidth',0.7,'Color',[0.25 0.25 0.25])
title('Contribution diagram')
xlabel('Frequency (Hz)')
ylabel('Normalized Mode Contribution')
xlim([0 50])
ylim([0 6.5])
legend('Not Stable', 'Stable in frequency/damping', 'Stable in freq.','')

%% Computing stable and most dominant modes
e_opt = double(sr200);

Ij_thresh_indx = find(Ij2 >= thresh_Ij); % most dominant modes
Ij_thresh = zeros(size(Ij2));
Ij_thresh(Ij_thresh_indx) = 1;

f_stable_indx = find(~finst); % stable modes
f_stable = zeros(size(fj));
f_stable(f_stable_indx) = 1;

idx = f_stable.*Ij_thresh';
nkk = find(idx);

[~, nk_sort] = sort(Ij2(nkk),'descend');
nk_stable_dom = nkk(nk_sort);

e_optnk = e_opt(nk_stable_dom);
f_dmd = abs(e_optnk)/(2*pi);
zeta_dmd = -real(e_optnk)./abs(e_optnk);
dmdnk = [f_dmd zeta_dmd];

e_opt_0 = e_optnk;

% save('../Data/initial_guess4optdmd_00_00_00_run1_lowqualityf500.mat','e_opt_0')
