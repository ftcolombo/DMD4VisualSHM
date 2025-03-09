function [dmds, Phis, x_modo, s_dmds, Ij] = dmd_free_time_response(X, fs, Td, na, rb, Tj)
% function dmd_free_time_response applies DMD to the time-resolved data

% INPUTS:
% X:        rows are assumed to be measurements, columns are assumed to be time points
% fs:       sampling rate of the measurements
% Td:       ( may be different than t(end) ) maximum time instant analyzed by the time embeddings
% na:       number of time embeddings considered (m = na - 1)
% rb:       user specified order of the reduced order model approximated by DMD
% Tj:       ( may be different than t(end) ) maximum time instant analyzed by the contribution of each mode

% OUTPUTS:
% dmds:     [f_dmds zeta_dmds], frequencies of oscillation in Hz and decay/growth rates estimated by DMD
% Phis:     modes estimated by DMD
% x_modo:   displacement over time of the DMD modes
% s_dmds:   eigenvalues estimated by DMD
% Ij:       contribution of each mode to the free response

% F. T. Colombo, September 2024

[X1, X2, ~] = time_embedding(X, na, 0:1/fs:Td);

clear X

[U,S,V] = svd(X1,'econ');

x0 = X1(:,1);

clear X1

r = min([rb size(S,1)]);

Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);

clear U S V

Vs = Vr/Sr;
Atilde = Ur'*X2*Vs;
[W,D] = eig(Atilde); 

clear Atilde

lambda = diag(D);

dt = 1/fs;

s_dmd = log(lambda)/dt;
f_dmd = abs(s_dmd)/(2*pi);
zeta_dmd = -real(s_dmd)./abs(s_dmd);

dmd = [f_dmd zeta_dmd];

% Sorting the DMD modes according to their contribution to the free response

alpha = (W\Ur')*x0; % Mode amplitude vector
clear Ur Sr Vr x0

Phi = X2*(Vs*W);
clear X2 Vs W

k = size(Phi,2);
Ij = zeros(1,k);
for j =  1:k
    for i = 1:fs*Tj
        Ij(j) = Ij(j) + abs(alpha(j)*(lambda(j))^(i-1));
    end

    Ij(j) = Ij(j)*norm(Phi(:,j),'fro'); % mode selection criterion
end

[~,ind] = sort(Ij, 'descend');

Ds = D(ind,ind);
Phis = Phi(:,ind);

clear Phi D

lambdas = diag(Ds);

clear Ds

s_dmds = log(lambdas)/dt;

f_dmds = abs(s_dmds)/(2*pi); 
zeta_dmds = -real(s_dmds)./abs(s_dmds);

dmds = [f_dmds zeta_dmds];

b = alpha(ind);

t = 0:1/fs:Tj;
x_modo = zeros(size(b,1),length(t));
for i = 1:length(t)
    x_modo(:,i) = b.*exp(s_dmds*t(i));
end
