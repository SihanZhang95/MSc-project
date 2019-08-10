alpha = 0.05;
Delta = 1;
C=0.2902;
v = [0.5;0.5];
B = [0.9, 0.8; 0.8, 0.9];
lambda0 = [0.5,0.5];

power = zeros(1,10);

for i = 1:5
    A = [0.5,0.2*i;0.2*i,0.5];
    string=zeros(1,100);
    for j = 1:100
        [N,X,E,M,E_source,lambda,spectra,cross_spectra,coherences] = MvHawkesSimulation(v,A,B,lambda0);
        Xt = transpose(X(1,:));
        Yt = transpose(X(2,:));
        
        [x,y] = TWCOH_quick(Xt',Yt',Delta,24,1,1);
        xx=(x(~isnan(x)));
        f=0;
           for jj = 1:length(xx)
               percentile=fzero(@(x) estcoh(x,y(jj),2,1-alpha,0), [0.000001,0.999999]);
                if (xx(jj)>percentile || xx(jj)==percentile)
                      f=f+1;
                end
           end
        string(j)=f/length(xx)>C;
        [i,j]
    end
    power(i)=mean(string);
end

plot(0.2*[1:5],power(1:5),'-o','color','red')
ylim([0 1])
ylim([0 1])
xlabel('Mutual excitation parameter')
ylabel('Power (%)')
title('Power-a plot of Multivariate Hawkes process using New Global Test')