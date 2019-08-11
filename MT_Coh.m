% % Multitaper AR(4) and AR(2) coherence % %

AR4model = arima('Constant',0,'AR',{2.7607,-3.8106,2.6535,-0.9238},'Variance',1);

AR2model = arima('Constant', 0, 'AR', {5/6,-1/6}, 'Variance', 1);

% number of realisations

N = 1024;

% sampling interval

Delta = 1;

% bandwidth

NW = 8;



% frequency in [0, 0.5]

f = linspace(0,1/2,N/2+1);

% dpss up to order 2NW

dps = dpss(N, NW);



% obeservations in AR(4) process

X = simulate(AR4model,N);

% obeservations in AR(2) process

Y = simulate(AR2model,N);



% taper time series

hX = dps.*X;

hY = dps.*Y;



dseX = [];

dseY = [];

csdeXY = [];



for i = 1:2*NW

    dseX(i,:) = J(hX(:,i),Delta,f).*conj(J(hX(:,i),Delta,f));

    dseY(i,:) = J(hY(:,i),Delta,f).*conj(J(hY(:,i),Delta,f));

    csdeXY(i,:) = J(hX(:,i),Delta,f).*conj(J(hY(:,i),Delta,f));

end



% apply multitaper

dseXmt = transpose(1./(1:2*NW)).*cumsum(dseX);

dseYmt = transpose(1./(1:2*NW)).*cumsum(dseY);

csdeXYmt = transpose(1./(1:2*NW)).*cumsum(csdeXY);



% find coherence

figure()

gammaHat = [];

for i = 1:NW

    gammaHat(i,:) = abs(csdeXYmt((2*i-1),:)).^2./(dseXmt((2*i-1),:).*dseYmt((2*i-1),:));

    subplot(NW/2,2,i);
    
    plot(f,gammaHat(i,:),'color','green');
  
    xlim([0 0.5]);

    ylim([0 1]);

    xlabel('f');

    ylabel(sprintf('K=%d',(2*i-1))) 
end



% obeservations in combined AR(4) and AR(2) processes

V = X+Y;

W = X+normrnd(0,1,[N,1]);



% (N*1) vector of sdf for frequency in [0, 0.5]

SX = ARsdf(4,[2.7607;-3.8106;2.6535;-0.9238],transpose(f));

SY = ARsdf(2,[5/6;-1/6],transpose(f));



% true coherence

gamma = abs(SX).^2./((SX+SY).*(SX+1));



% taper time series

hW = dps.*V;

hZ = dps.*W;



dseW = [];

dseZ = [];

csdeWZ = [];



for i = (1:2*NW)

    dseW(i,:) = J(hW(:,i),Delta,f).*conj(J(hW(:,i),Delta,f));

    dseZ(i,:) = J(hZ(:,i),Delta,f).*conj(J(hZ(:,i),Delta,f));

    csdeWZ(i,:) = J(hW(:,i),Delta,f).*conj(J(hZ(:,i),Delta,f));

end



% apply multitaper

dseWmt = transpose(1./(1:2*NW)).*cumsum(dseW);

dseZmt = transpose(1./(1:2*NW)).*cumsum(dseZ);

csdeWZmt = transpose(1./(1:2*NW)).*cumsum(csdeWZ);



% find coherence:,i

figure()

gammaHat = [];

for i = (1:NW)

    gammaHat(i,:) = abs(csdeWZmt((2*i-1),:)).^2./(dseWmt((2*i-1),:).*dseZmt((2*i-1),:));

    subplot(NW/2,2,i)

    plot(f,gammaHat(i,:),'k');

    hold on

    plot(f,gamma,'r');

    p(2).Marker = '*';
    
    hold off

    xlim([0 0.5]);

    ylim([0 1]);

    xlabel('f');

    ylabel(sprintf('K=%d',(2*i-1)))

end