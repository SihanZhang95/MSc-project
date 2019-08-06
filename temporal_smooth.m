function [smoothed_scalo] = temporal_smooth(cwt_x,cwt_y,a,OMEGA,delt,d)
% perform smoothing on desired section
N = size(cwt_x);
N = N(2);
a_0 = a./delt;
%nu = ceil(a*40);
%     % calculate smoothing parameter delta   
M = ceil(OMEGA*abs(a_0));
nu = ceil(3*abs(a_0)*d);

smoothed_scalo = NaN(size(cwt_x));

for jj = 1:length(M)


    b_start = nu(jj) + M(jj);
    %b_start = b_start*delt;
    b_end = N - 1 - nu(jj) - M(jj);
    %b_end = b_end*delt;
    for b = b_start:b_end;
        tau = [b-M(jj):b+M(jj)];
        smoothed_scalo(jj,b) = mean(cwt_x(jj,tau).*conj(cwt_y(jj,tau)));
    end
end