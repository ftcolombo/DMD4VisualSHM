close all
clear all
clc

addpath('../Data')
addpath('../Functions')

% properties of the clamped-free beam made of Polypropylene: 
% lbeam = 0.8; % m, length of the beam
% wbeam = 25.4e-3; % m, width/thickness of the beam
% E = 1300e6; % Pa, Young's modulus
% density = 900; % kg/m^3, density 
% epsilon = 0.05; % damping factor

load displ_beam_impulse025N_damp5p.mat
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

f_noise = 1; % flag for adding noise ( 1: noise is added, 0: no noise )
SNR = 40; % signal-to-noise ratio (SNR) in decibels to be added if f_noise is true

if f_noise == 0
    yy = y;
else
    yy = 0*y;
    rng(1,'twister')
    for i = 1:size(y,2)
        noise = randn(size(y,1),1)*std(y(:,i))/db2mag(SNR);
        yy(:,i) = y(:,i) + noise;
    end
end

%% Comparing the error of the frequency and damping for an increasing number of time embeddings

f_dmd = zeros(nm,1);
zeta_dmd = zeros(nm,1);

rb = 400; % maximum number of singular values to be kept in DMD
na = 1:21; % number of time embeddings considered (m = na - 1)
T = 1; % dmd will analyze the data during the first T seconds

j = 0;
for nai = na
    [dmds,~,~,~,~] = dmd_free_time_response(yy', fs, T, nai, rb, T);
    j = j + 1;
    for i = 1:nm
        f_dmd(i,j) = dmds(2*i-1,1);
        zeta_dmd(i,j) = dmds(2*i-1,2);
    end
end

e_f = 0*f_dmd;
e_zeta = 0*zeta_dmd;
for i = 1:size(f_dmd,2)
    e_f(:,i) = abs(f_dmd(:,i) - f_exact)./f_exact;
    if ~logical(epsilon)
        e_zeta(:,i) = abs(zeta_dmd(:,i) - zeta_exact);
    else
        e_zeta(:,i) = abs(zeta_dmd(:,i) - zeta_exact)./zeta_exact;
    end
end

% % % Fig 2a
figure,
semilogy(na-1, e_f(1,:),'Color',[0 0 0.64],'LineStyle','--'), hold on
semilogy(na-1, e_f(2,:),'Color',[0.74 0.16 0.17],'LineStyle','-.')
semilogy(na-1, e_f(3,:),'Color',[0.91 0.78 0.086],'Marker','o','MarkerFaceColor',[0.91 0.78 0.08],'MarkerSize',4)
xlabel('Number of time embeddings')
ylabel('Relative Error')
title('Frequency of oscillation')
grid on
legend('1','2','3')

% % % Fig 2b
figure,
semilogy(na-1, e_zeta(1,:),'Color',[0 0 0.64],'LineStyle','--'), hold on
semilogy(na-1, e_zeta(2,:),'Color',[0.74 0.16 0.17],'LineStyle','-.')
semilogy(na-1, e_zeta(3,:),'Color',[0.91 0.78 0.086],'Marker','o','MarkerFaceColor',[0.91 0.78 0.08],'MarkerSize',4)
xlabel('Number of time embeddings')
ylabel('Relative Error')
title('Decay rate')
grid on
legend('1','2','3')

%% Comparing the results for a certain amount of noise and number of time embeddings

na = 8; % number of time embeddings considered (m = na - 1)

[dmds, Phis, x_modo, ~, ~] = dmd_free_time_response(yy', fs, T, na, rb, T);

f_dmd = zeros(nm,1);
zeta_dmd = zeros(nm,1);
for i = 1:nm
    f_dmd(i) = dmds(2*i-1,1);
    zeta_dmd(i) = dmds(2*i-1,2);
    
    e_f2(i) = abs(f_dmd(i) - f_exact(i))./f_exact(i);
    e_zeta2(i) = abs(zeta_dmd(i) - zeta_exact(i))/zeta_exact(i);
end

% % % Table 1
ff = [f_exact f_dmd 100*e_f2'];
zetazeta = [zeta_exact zeta_dmd 100*e_zeta2'];

%% By stacking past measurements to enrich the dataset, Phis must be restricted to now only consider the measurement on each snapshot

nk = 1:2:2*nm; % indexes corresponding to Modes 1, 2, 3, 4 and 5
phi_dmd = Phis(1:size(y,2), nk);

x = 0:lbeam/size(y,2):lbeam;

% % % Fig 3
figure, 
hold on  
for i = 1:nm-1
    c_d = 1/phi(end,i);
    c_dmd = 1/real(phi_dmd(end,i));
    plot(x, [0; phi(:,i)]*c_d,'LineWidth',2,'Color',[0.8 0.8 0.8])
    switch i
        case 1,plot(x, [0; real(phi_dmd(:,i))]*c_dmd,'LineWidth',1,'LineStyle','--','Color',[0.00 0.00 0.64],'DisplayName','1')
        case 2,plot(x, [0; real(phi_dmd(:,i))]*c_dmd,'LineWidth',1,'LineStyle','-.','Color',[0.74 0.16 0.17],'DisplayName','2')
        case 3,plot(x, [0; real(phi_dmd(:,i))]*c_dmd,'LineWidth',1,'LineStyle','--','Color',[0.91 0.78 0.08],'Marker','o','MarkerFaceColor',[0.91 0.78 0.08],'MarkerSize',4,'DisplayName','3')
        case 4,plot(x, [0; real(phi_dmd(:,i))]*c_dmd,'LineWidth',0.8,'LineStyle','-','Color',[0.31 0.68 0.62],'DisplayName','4')    
    end
end
grid on
xlabel('x-direction (m)')
ylabel('Normalized y-displacement')
title('Mode shapes')
f=get(gca,'Children');
legend([f(7),f(5),f(3),f(1)],'Location','best')
ylim([-1.2 1.2])


ts = 0:1/fs:(size(x_modo,2)-1)/fs;
x_dmd = (Phis(1:size(y,2),1:2*(nm-1))*x_modo(1:2*(nm-1),:))';

figure,
subplot(2,1,1), hold on
plot(t, 1000*y(:,1),'LineWidth',2,'Color',[0.8 0.8 0.8])
plot(ts, 1000*real(x_dmd(:,1)),'LineStyle','--','Color',[0.00 0.00 0.64])
grid on
xlim([0 1])
xlabel('Time (s)')
ylabel('y-displacement (mm)')
title('(a) near fixed end')
legend('original','reconstructed','Location','best')
subplot(2,1,2), hold on
plot(t, 1000*y(:,end),'LineWidth',2,'Color',[0.8 0.8 0.8])
plot(ts, 1000*real(x_dmd(:,end)),'LineStyle','--','Color',[0.00 0.00 0.64])
grid on
xlim([0 1])
yticks([-50 -25 0 25 50])
xlabel('Time (s)')
ylabel('y-displacement (mm)')
title('(b) free end')
