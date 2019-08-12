alpha = 0.05;
Delta = 1;
kappa= 24;
C=0.2902;
n=100; % number of repetitions
m=6; % number of trials we have

v = [0.5;0.5];
B = [0.9, 0.8; 0.8, 0.9];
lambda0 = [0.5,0.5];

power = zeros(1,5);

for i = 1:5
    A = [0.5,0.2*i;0.2*i,0.5];
    string=zeros(1,n);
    
for j = 1:n
    sumTWCOH=0;
    for jj = 1:m
        [N,X,E,M,E_source,lambda,spectra,cross_spectra,coherences] = MvHawkesSimulation(v,A,B,lambda0);

        Xt = transpose(X(1,:));
        Yt = transpose(X(2,:));
        
        [x,y] = TWCOH_quick(Xt',Yt',Delta,kappa,1,1);
        sumTWCOH=sumTWCOH+x;
    end
        average=sumTWCOH/m;
        xx=(average(~isnan(average)));
        f=0;
           for jj = 1:length(xx)
               percentile=fzero(@(x) estcoh(x,y(jj),2,1-alpha,0), [0.000001,0.999999]);
                if (xx(jj)>percentile || xx(jj)==percentile)
                      f=f+1;
                end
           end
        string(j)=(f/length(xx))>C;
        [i,j]
end
    power(i)=mean(string);
end

