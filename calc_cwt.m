function [cwt_x] = calc_cwt(x,delt,a,d)

N = length(x);

% Fourier transfrom signal x
x_hat = (fft(x));

% Calculate signal nyquist frequency
Fn = 1/delt;

% calc frequency vector
l = [0:N-1];
f = l./N;

% calc scale coefficient a_{0}
a_0 = a./delt;

% How many scales
nscales = length(a);

% calc f_{0} vector with which to calculate wavelet transform frequency responses
% corresponding to scale coefficient a_{0}
f_pos = f(1:N/2);
f_neg = (f(N/2+1:end)-1);

f_0 = [f_pos f_neg];
f_0= repmat(a_0',1,N).*repmat(f_0,nscales,1);

% Initiate the wavelet transform matrices
cwt_x = zeros(nscales,N);


    % Calculate frequency response of analytic wavelet at frequencies in f_{0}
    % Do in matrix form, each row represents
    Psi = (pi^(1/4))*((2*d)^(1/2))*exp(-2*(d*pi*(f_0-1)).^2);
    
%     for ii = 1:10;
%         Psi(ii,:) = ifft(Psi(ii,:));
%         Psi(ii,(ceil(4*a_0(ii))+1):(end-ceil(3*a_0(ii))-1))=0;
%         Psi(ii,:) = fft(Psi(ii,:));
%     end
    
    % Calculate the products
    x_hat_Psi     = repmat(x_hat,length(a_0),1).*Psi;

    % For each scale
    for jj = 1:length(a_0);

        % calculate the cwts for each type of wavelet
        cwt_x(jj,:)       = (sqrt(abs(a_0(jj)))/(N*sqrt(delt))) * ifft(x_hat_Psi(jj,:));

    end
    
%     testcwt_x = zeros(nscales,N);
%     WW = zeros(1,N);
%     b_0 = [0:N-1];
%     for ii = 1:length(a_0);
%         for jj = 1:length(b_0);
%             for nn = 0:N-1;
%                 WW(nn+1) = sum(Psi(ii,:).*exp(i*2*pi*[0:N-1]*(b_0(jj) - nn)/N));
%             end
%             testcwt_x(ii,jj) = ((abs(a)^1/2)/N)*sum(x.*WW); 
%         end
%     end
    