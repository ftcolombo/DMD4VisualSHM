clc
clear all
close all

addpath('../Data')

%% Computing the distribution of the eigenvalues in the s-plane for all damage scenarios

np = [4 3 2 1]; % ordering of the modes according to the natural frequencies
nscen = 8; 

pdf_imag = zeros(2*nscen,length(np));
pdf_real = zeros(2*nscen,length(np));

load('pdf_lambda_00_00_00_run1_noise')
pdf_imag(1:2,:) = pdf_lambda_vec_imag(:,np);
pdf_real(1:2,:) = pdf_lambda_vec_real(:,np);

load('pdf_lambda_05_00_00_run1_noise.mat')
pdf_imag(3:4,:) = pdf_lambda_vec_imag(:,np);
pdf_real(3:4,:) = pdf_lambda_vec_real(:,np);

load('pdf_lambda_10_00_00_run1_noise.mat')
pdf_imag(5:6,:) = pdf_lambda_vec_imag(:,np);
pdf_real(5:6,:) = pdf_lambda_vec_real(:,np);

load('pdf_lambda_13_00_00_run1_noise.mat')
pdf_imag(7:8,:) = pdf_lambda_vec_imag(:,np);
pdf_real(7:8,:) = pdf_lambda_vec_real(:,np);

load('pdf_lambda_13_05_00_run2_noise.mat')
pdf_imag(9:10,:) = pdf_lambda_vec_imag(:,np);
pdf_real(9:10,:) = pdf_lambda_vec_real(:,np);

load('pdf_lambda_13_05_05_run2_noise.mat')
pdf_imag(11:12,:) = pdf_lambda_vec_imag(:,np);
pdf_real(11:12,:) = pdf_lambda_vec_real(:,np);

load('pdf_lambda_13_10_05_run2_noise.mat')
pdf_imag(13:14,:) = pdf_lambda_vec_imag(:,np);
pdf_real(13:14,:) = pdf_lambda_vec_real(:,np);

load('pdf_lambda_13_10_11_run2_noise.mat')
pdf_imag(15:16,:) = pdf_lambda_vec_imag(:,np);
pdf_real(15:16,:) = pdf_lambda_vec_real(:,np);

theta_eli = linspace(0,2*pi,100);

f0 = sqrt(pdf_real(1:2:end,1).^2 + pdf_imag(1:2:end,1).^2)/(2*pi);
f1 = sqrt(pdf_real(1:2:end,2).^2 + pdf_imag(1:2:end,2).^2)/(2*pi);
f2 = sqrt(pdf_real(1:2:end,3).^2 + pdf_imag(1:2:end,3).^2)/(2*pi);
f3 = sqrt(pdf_real(1:2:end,4).^2 + pdf_imag(1:2:end,4).^2)/(2*pi);

zeta0 = -pdf_real(1:2:end,1)./sqrt(pdf_real(1:2:end,1).^2 + pdf_imag(1:2:end,1).^2);
zeta1 = -pdf_real(1:2:end,2)./sqrt(pdf_real(1:2:end,2).^2 + pdf_imag(1:2:end,2).^2);
zeta2 = -pdf_real(1:2:end,3)./sqrt(pdf_real(1:2:end,3).^2 + pdf_imag(1:2:end,3).^2);
zeta3 = -pdf_real(1:2:end,4)./sqrt(pdf_real(1:2:end,4).^2 + pdf_imag(1:2:end,4).^2);

wn = [f0 f1 f2 f3]*2*pi;
zeta = [zeta0 zeta1 zeta2 zeta3];

% % % Figure 17 
figure, subplot(1,3,1), hold on

j = 2;
minz = min(zeta(:,j));
maxz = max(zeta(:,j));
theta_circ = linspace(acos(-0.9*minz),acos(-1.1*maxz),10);

nwn = 4;
wn_circ = linspace(1.01*wn(1,j),0.99*wn(8,j),nwn);
for i = 1:nwn
    x_circ = wn_circ(i)*cos(theta_circ);
    y_circ = wn_circ(i)*sin(theta_circ);
    plot(x_circ, y_circ,'--','Color',[0.8 0.8 0.8])
    text(x_circ(1),mean(y_circ),[num2str(wn_circ(i)/(2*pi),'%.1f'),' Hz'],'BackgroundColor','w')
end

minw = min(wn(:,j));
maxw = max(wn(:,j));

ntn = 4;
theta_circ = linspace(acos(-0.95*minz),acos(-1.05*maxz),ntn);
for i = 1:ntn
    x_line = [1.03*maxw*cos(theta_circ(i)) 0.99*minw*cos(theta_circ(i))];
    y_line = [1.03*maxw*sin(theta_circ(i)) 0.99*minw*sin(theta_circ(i))];
    plot(x_line, y_line,'--','Color',[0.8 0.8 0.8])
    text(x_line(1),y_line(1),[num2str(-cos(theta_circ(i))*100,'%.2f'),'%'],'BackgroundColor','w')
end

for i = 1:nscen
    x_eli = pdf_real(2*i-1,j) + 2*pdf_real(2*i,j)*sin(theta_eli);
    y_eli = pdf_imag(2*i-1,j) + 2*pdf_imag(2*i,j)*cos(theta_eli);
    p1 = patch('YData', y_eli,'XData',x_eli, 'LineStyle','none','FaceColor',[0.65 0.65 0.65]);
    
    p1.FaceVertexAlphaData = 0.95;
    p1.FaceAlpha = 'flat';
    
    x_eli = pdf_real(2*i-1,j) + pdf_real(2*i,j)*sin(theta_eli);
    y_eli = pdf_imag(2*i-1,j) + pdf_imag(2*i,j)*cos(theta_eli);
    p2 = patch('YData', y_eli,'XData',x_eli, 'LineStyle','none','FaceColor',[0.2 0.3 0.75]);
 
    p2.FaceVertexAlphaData = 0.95;
    p2.FaceAlpha = 'flat';  
    
    scatter(pdf_real(2*i-1,j),pdf_imag(2*i-1,j),'MarkerEdgeColor',[1 0 0])%,'MarkerFaceColor',[i/10 0 0],'MarkerEdgeColor','none')
    switch i
        case 1
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'A')
        case 2
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'B')
        case 3
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'C')    
        case 4
            text(1.014*pdf_real(2*i-1,j),1.014*pdf_imag(2*i-1,j),'D')
        case 5
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'E')
        case 6
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'F') 
        case 7
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'G')
        case 8
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'H')
    end         
end

xlim([-0.82 -0.35])
ylim([33 46])
xlabel('Real axis')
ylabel('Imaginary axis')
title('Mode 1')

subplot(1,3,2), hold on

j = 3;
minz = min(zeta(:,j));
maxz = max(zeta(:,j));
theta_circ = linspace(acos(-0.9*minz),acos(-1.1*maxz),10);

nwn = 4;
wn_circ = linspace(1.01*wn(1,j),0.99*wn(8,j),nwn);
for i = 1:nwn
    x_circ = wn_circ(i)*cos(theta_circ);
    y_circ = wn_circ(i)*sin(theta_circ);
    plot(x_circ, y_circ,'--','Color',[0.8 0.8 0.8])
    text(x_circ(1),mean(y_circ),[num2str(wn_circ(i)/(2*pi),'%.1f'),' Hz'],'BackgroundColor','w')
end

minw = min(wn(:,j));
maxw = max(wn(:,j));

ntn = 4;
theta_circ = linspace(acos(-0.95*minz),acos(-1.05*maxz),ntn);
for i = 1:ntn
    x_line = [1.03*maxw*cos(theta_circ(i)) 0.99*minw*cos(theta_circ(i))];
    y_line = [1.03*maxw*sin(theta_circ(i)) 0.99*minw*sin(theta_circ(i))];
    plot(x_line, y_line,'--','Color',[0.8 0.8 0.8])
    text(x_line(1),y_line(1),[num2str(-cos(theta_circ(i))*100,'%.2f'),'%'],'BackgroundColor','w')
end

for i = 1:nscen
    x_eli = pdf_real(2*i-1,j) + 2*pdf_real(2*i,j)*sin(theta_eli);
    y_eli = pdf_imag(2*i-1,j) + 2*pdf_imag(2*i,j)*cos(theta_eli);
    p1 = patch('YData', y_eli,'XData',x_eli, 'LineStyle','none','FaceColor',[0.65 0.65 0.65]);
    
    p1.FaceVertexAlphaData = 0.95;
    p1.FaceAlpha = 'flat';
    
    x_eli = pdf_real(2*i-1,j) + pdf_real(2*i,j)*sin(theta_eli);
    y_eli = pdf_imag(2*i-1,j) + pdf_imag(2*i,j)*cos(theta_eli);
    p2 = patch('YData', y_eli,'XData',x_eli, 'LineStyle','none','FaceColor',[0.2 0.3 0.75]);
 
    p2.FaceVertexAlphaData = 0.95;
    p2.FaceAlpha = 'flat';    

    scatter(pdf_real(2*i-1,j),pdf_imag(2*i-1,j),'MarkerEdgeColor',[1 0 0])%,'MarkerFaceColor',[i/10 0 0],'MarkerEdgeColor','none')
    switch i
        case 1
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'A')
        case 2
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'B')
        case 3
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'C')    
        case 4
            text(1.014*pdf_real(2*i-1,j),1.014*pdf_imag(2*i-1,j),'D')
        case 5
            text(1.01*pdf_real(2*i-1,j),0.99*pdf_imag(2*i-1,j),'E')
        case 6
            text(0.992*pdf_real(2*i-1,j),1.015*pdf_imag(2*i-1,j),'F') 
        case 7
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'G')
        case 8
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'H')
    end         
end

xlim([-1.9 -0.7])
ylim([65 91])
xlabel('Real axis')
% ylabel('Imaginary axis')
title('Mode 1 - 1st H')

subplot(1,3,3), hold on

j = 4;
minz = min(zeta(:,j));
maxz = max(zeta(:,j));
theta_circ = linspace(acos(-0.85*minz),acos(-1.15*maxz),10);

nwn = 4;
wn_circ = linspace(1.012*wn(1,j),0.988*wn(8,j),nwn);
for i = 1:nwn
    x_circ = wn_circ(i)*cos(theta_circ);
    y_circ = wn_circ(i)*sin(theta_circ);
    plot(x_circ, y_circ,'--','Color',[0.8 0.8 0.8])
    text(x_circ(1),mean(y_circ),[num2str(wn_circ(i)/(2*pi),'%.1f'),' Hz'],'BackgroundColor','w')
end

minw = min(wn(:,j));
maxw = max(wn(:,j));

ntn = 4;
theta_circ = linspace(acos(-0.85*minz),acos(-1.15*maxz),ntn);
for i = 1:ntn
    x_line = [1.03*maxw*cos(theta_circ(i)) 0.99*minw*cos(theta_circ(i))];
    y_line = [1.03*maxw*sin(theta_circ(i)) 0.99*minw*sin(theta_circ(i))];
    plot(x_line, y_line,'--','Color',[0.8 0.8 0.8])
    text(x_line(1),y_line(1),[num2str(-cos(theta_circ(i))*100,'%.2f'),'%'],'BackgroundColor','w')
end

for i = 1:nscen
    x_eli = pdf_real(2*i-1,j) + 2*pdf_real(2*i,j)*sin(theta_eli);
    y_eli = pdf_imag(2*i-1,j) + 2*pdf_imag(2*i,j)*cos(theta_eli);
    p1 = patch('YData', y_eli,'XData',x_eli, 'LineStyle','none','FaceColor',[0.65 0.65 0.65]);
    
    p1.FaceVertexAlphaData = 1;
    p1.FaceAlpha = 'flat';
    
    x_eli = pdf_real(2*i-1,j) + pdf_real(2*i,j)*sin(theta_eli);
    y_eli = pdf_imag(2*i-1,j) + pdf_imag(2*i,j)*cos(theta_eli);
    p2 = patch('YData', y_eli,'XData',x_eli, 'LineStyle','none','FaceColor',[0.2 0.3 0.75]);
 
    p2.FaceVertexAlphaData = 1;
    p2.FaceAlpha = 'flat';      

    scatter(pdf_real(2*i-1,j),pdf_imag(2*i-1,j),'MarkerEdgeColor',[1 0 0])%,'MarkerFaceColor',[i/10 0 0],'MarkerEdgeColor','none')
         
    switch i
        case 1
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'A')
        case 2
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'B')
        case 3
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'C')    
        case 4
            text(1.014*pdf_real(2*i-1,j),1.014*pdf_imag(2*i-1,j),'D')
        case 5
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'E')
        case 6
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'F') 
        case 7
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'G')
        case 8
            text(0.992*pdf_real(2*i-1,j),0.992*pdf_imag(2*i-1,j),'H')
    end    
end

xlim([-3.1 -0.8])
ylim([98 137])
xlabel('Real axis')
% ylabel('Imaginary axis')
title('Mode 1 - 2nd H')
