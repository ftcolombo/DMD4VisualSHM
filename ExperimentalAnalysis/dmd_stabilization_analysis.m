clc
clear all
close all

addpath('../Data')

% you can run the script 'dmd_stabilization_data' to obtain the eigenvalues
% estimated by DMD for model orders r = 175:1:250, if you have access to the
% video '00_00_00_run1.mp4' of the beam by Garrido et al.
% (https://doi.org/10.1016/j.ymssp.2023.110539)

% Otherwise, you can load the data below that includes the eigenvalues computed
load('stabilizationdata_00_00_00_run1I.mat')

%% Evaluating the stability of the eigenvalues computed by DMD for an increased order r of the SVD

o_min = 180; % analysing only the models with order >= o_min
f_max = 50; % restricting the frequency range to [ 0 ; f_max ]

nr = length(kk);

fr = abs(sr)/(2*pi);
dr = 0*fr;
for i = 1:nr
    dr(:,i) = -real(sr(:,i))./abs(sr(:,i));
end

fst = false(size(fr));
dst = false(size(fr));

tol = [0.005 0.005];
for i = 1:nr-1
    f0 = fr(1:kk(i), i);
    f1 = fr(1:kk(i+1), i+1);
    d0 = dr(1:kk(i), i);
    d1 = dr(1:kk(i+1), i+1);    
    for j = 1:kk(i)
        if abs(f0(j)) < tol(1) && min(abs(f0(j)-f1(:))) < tol(1)
            fst(j,i) = true;
        end
        if min(abs(f0(j)-f1(:))) < tol(1)*f0(j)
            fst(j,i) = true;
        end
        if min(abs(d0(j)-d1(:))) < tol(2)*abs(d0(j))
            dst(j,i) = true;
        end
    end
end


i_fd = 0;
i_f = 0;
i_ns = 0;
for i = 1:nr
    for j = 1:kk(i)
        if fst(j,i) == true
            if dst(j,i) == true
                i_fd = i_fd+1;
                v_fd(1,i_fd) = fr(j,i);
                v_fd(2,i_fd) = kk(i);
            else
                i_f = i_f+1;
                v_f(1,i_f) = fr(j,i);
                v_f(2,i_f) = kk(i);
            end
        else
            i_ns = i_ns+1;
            v_ns(1,i_ns) = fr(j,i);
            v_ns(2,i_ns) = kk(i);
        end
    end
end

a_ns = find(v_ns(2,:) >= o_min,1);
a_fd = find(v_fd(2,:) >= o_min,1);
a_f = find(v_f(2,:) >= o_min,1);

v_ns2 = v_ns(:,a_ns:end);
v_fd2 = v_fd(:,a_fd:end);
v_f2 = v_f(:,a_f:end);

j = 1;
for i = 1:size(v_ns2,2)
    if v_ns2(1,i) <= f_max
        v_ns3(:,j) = v_ns2(:,i);
        j = j+1;
    end
end

j = 1;
for i = 1:size(v_fd2,2)
    if v_fd2(1,i) <= f_max
        v_fd3(:,j) = v_fd2(:,i);
        j = j+1;
    end
end

j = 1;
for i = 1:size(v_f2,2)
    if v_f2(1,i) <= f_max
        v_f3(:,j) = v_f2(:,i);
        j = j+1;
    end
end

% % % Fig 11
figure, hold on
plot(v_ns3(1,:),v_ns3(2,:),'.', 'MarkerSize', 6, 'Color', [0.7 0.7 0.7])
plot(v_fd3(1,:),v_fd3(2,:),'x', 'MarkerSize', 5, 'Color', [0.74 0.16 0.17],'MarkerFaceColor', [0.74 0.16 0.17]) 
plot(v_f3(1,:),v_f3(2,:),'o', 'MarkerSize', 2, 'Color', [0 0 0.55])           
legend('Not Stable', 'Stable in frequency/damping', 'Stable in freq.')
legend('Orientation','horizontal','Location','best')
grid on
xlabel('Frequency (Hz)')
ylabel('Model Order')
title('Stability Diagram')
ylim([180 259])

