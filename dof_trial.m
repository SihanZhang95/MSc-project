function [dof] = dof_trial(fd,N,delt,kappa,d)
fn = 1/(2*delt);

amin = fd/fn;
amax = delt*(N-5)/(6*d + 2*kappa);

a = linspace(amin,amax,100);

a0=a./delt;
NB=2.*ceil(abs(a0).*kappa)+1;
NS=2.*ceil(3.*abs(a0).*d)+1;
NP=NB+NS-1;
dof=zeros(1,length(a0));
for i = 1:length(a0)
g = arrayfun(@(n) exp(-((n-ceil(3*abs(a0(i))*d))/(a0(i)*d))^2/2),0:(NS(i)-1));
% sum of squared values
 C = sum(g.*g);
% normalised
 g = g/sqrt(C);

B=zeros(NP(i),NB(i));
  for j = 0:(NB(i)-1)
    a=zeros(1,j);
    b=g;
    c=zeros(1,(NB(i)-1-j));
    B(:,(j+1))=[a,b,c];
  end
  B=B/(sqrt(NB(i)));
  OP = B * transpose(B);
  
  % find largest NB eigenvalues and eigenvectors
  options.disp=0; % suppresses output
  [u,D,flag]=eigs(OP,NB(i),'la',options);
   if flag ~= 0
    error(' failure of eigenvector routine eigs');
   end
% the vectors have norm unity:
%
% now sort out polarity of eigenvectors

% eigenvalues
  lambda=diag(D);
%%% calculate degrees of freedom
 dof(i) = 1/sum(lambda.^2);
end
end