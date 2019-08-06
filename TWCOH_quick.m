function [coherence,dof] = TWCOH_quick(X,Y,delt,kappa,d,aorf)
% aorf indicates whether wavelet coherence is to be calculated with
% frequency of scale as the key parameter of interest. 1 is scale, 0 is
% frequency.

% determine number of data points
N = length(X);
NY = length(Y);
% error message if different lengths
if N ~= NY
    error('Time series X and Y are different lengths. Abort')
end

% some wavelet parameters
fd = 1.8;
fn = 1./(2*delt);

amin = fd/fn;
amax = delt*(N-5)/(6*d + 2*kappa);

if aorf
a = linspace(amin,amax,100);
else
fmin = 1/amax;
fmax = 1/amin;
f = linspace(fmax,fmin,100);
a = 1./f;
end

[cwt_X] = calc_cwt(X,delt,a,d);
[cwt_Y] = calc_cwt(Y,delt,a,d);

[smoothed_scaloXX] = temporal_smooth(cwt_X,cwt_X,a,kappa,delt,d);
[smoothed_scaloYY] = temporal_smooth(cwt_Y,cwt_Y,a,kappa,delt,d);
[smoothed_scaloXY] = temporal_smooth(cwt_X,cwt_Y,a,kappa,delt,d);
        
coherence = (abs(smoothed_scaloXY).^2)./(smoothed_scaloXX.*smoothed_scaloYY);

a0=a/delt;
NB=2*ceil(abs(a0)*kappa)+1;
NS=2*ceil(3*abs(a0)*d)+1;
NP=NB+NS-1;
phi=zeros(1,length(a0));
for i = 1:length(a0)
g=arrayfun(@(n) pi^(-1/4)*(abs(a0(i))*d)^(-1/2)*exp(-((n-ceil(3*abs(a0(i))*d))/a0(i)*d)^2/2),0:(NS(i)-1));
B=zeros(NP(i),NB(i));
  for j = 0:(NB(i)-1)
    a=zeros(1,j);
    b=g;
    c=zeros(1,(NB(i)-1-j));
    B(:,(j+1))=[a,b,c];
  end
  B=B/(sqrt(NB(i)));
  opm = B * transpose(B);
  evs=eig(opm);
  phi(i)=1/sum(evs(evs>0).^2);
end

  count=sum(~isnan(coherence));
  dof=zeros(1,length(sum(count)));
  jj=1;
  for ii = 1:length(count)
  if count(ii)~=0
     dof(jj:(jj+(count(ii)-1)))=phi(1:count(ii));
     jj=jj+count(ii);
  end
  end
end


