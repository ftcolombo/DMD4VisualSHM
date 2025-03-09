clc
clear all
close all

addpath('../Data')

% you can run the script 'dmd_stabilization_data' to obtain the eigenvalues
% estimated by DMD for 0:8 time embeddings, if you have access to the
% video '00_00_00_run1.mp4' of the beam by Garrido et al.
% (https://doi.org/10.1016/j.ymssp.2023.110539)

% Otherwise, you can load the data below that includes the eigenvalues computed
load('num_time_embed_data_00_00_00_run1I.mat')

%% Evaluating the location of the eigenvalues computed by DMD for an increased number of time embeddings

p = s_dmdk(2:3,:);

% % % Fig 10
figure, hold on
plot(p(:,1),'MarkerSize',13,'Marker','x','LineWidth',2,'LineStyle','none','DisplayName','0')
plot(p(:,2),'MarkerSize',13,'Marker','x','LineWidth',2,'LineStyle','none','DisplayName','1')
plot(p(:,3),'MarkerSize',13,'Marker','x','LineWidth',2,'LineStyle','none','DisplayName','2')
plot(p(:,4),'MarkerSize',13,'Marker','x','LineWidth',2,'LineStyle','none','DisplayName','3')
plot(p(:,5),'MarkerSize',13,'Marker','x','LineWidth',2,'LineStyle','none','DisplayName','4')
plot(p(:,6),'MarkerSize',13,'Marker','x','LineWidth',2,'LineStyle','none','DisplayName','5')
plot(p(:,7),'MarkerSize',13,'Marker','x','LineWidth',2,'LineStyle','none','DisplayName','6')
plot(p(:,8),'MarkerSize',13,'Marker','x','LineWidth',2,'LineStyle','none','DisplayName','7')

r1 = 43; r2 = 45;
zeta = 0.03;
theta1 = acos(-zeta);
plot([r1 r2]*cos(theta1),[r1 r2]*sin(theta1),'LineStyle','--','Color',[0.5 0.5 0.5],'DisplayName',num2str(zeta))
zeta = 0.05;
theta1 = acos(-zeta);
plot([r1 r2]*cos(theta1),[r1 r2]*sin(theta1),'LineStyle','--','Color',[0.5 0.5 0.5],'DisplayName',num2str(zeta))
zeta = 0.07;
theta1 = acos(-zeta);
plot([r1 r2]*cos(theta1),[r1 r2]*sin(theta1),'LineStyle','--','Color',[0.5 0.5 0.5],'DisplayName',num2str(zeta))
zeta = 0.09;
theta1 = acos(-zeta);
plot([r1 r2]*cos(theta1),[r1 r2]*sin(theta1),'LineStyle','--','Color',[0.5 0.5 0.5],'DisplayName',num2str(zeta))
zeta = 0.04;
theta1 = acos(-zeta);
plot([r1 r2]*cos(theta1),[r1 r2]*sin(theta1),'LineStyle','--','Color',[0.5 0.5 0.5])
zeta = 0.06;
theta1 = acos(-zeta);
plot([r1 r2]*cos(theta1),[r1 r2]*sin(theta1),'LineStyle','--','Color',[0.5 0.5 0.5])
zeta = 0.08;
theta1 = acos(-zeta);
plot([r1 r2]*cos(theta1),[r1 r2]*sin(theta1),'LineStyle','--','Color',[0.5 0.5 0.5])

ang = pi/2+(2*pi/180 : 0.01 : 7*pi/180); 
f = 6.9; r = f*2*pi;
plot(r*cos(ang),r*sin(ang),'LineStyle','--','Color',[0.5 0.5 0.5],'DisplayName',num2str(f));
f = 6.95; r = f*2*pi;
plot(r*cos(ang),r*sin(ang),'LineStyle','--','Color',[0.5 0.5 0.5],'DisplayName',num2str(f));
f = 7; r = f*2*pi;
plot(r*cos(ang),r*sin(ang),'LineStyle','--','Color',[0.5 0.5 0.5],'DisplayName',num2str(f));
f = 7.05; r = f*2*pi;
plot(r*cos(ang),r*sin(ang),'LineStyle','--','Color',[0.5 0.5 0.5],'DisplayName',num2str(f));

xlim([-4 -1])
ylim([43.25 44.25])
xlabel('Real Axis')
ylabel('Imaginary Axis')
