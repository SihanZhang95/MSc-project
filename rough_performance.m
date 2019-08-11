v = [0.5;0.5;0.5;0.5;0.5];

A = [0.1, 0, 0, 0, 0.1;0, 0.1, 0, 0, 0;0, 0, 0.1, 0, 0;0, 0, 0, 0.1, 0;0.1, 0, 0, 0, 0.1];

B = [0.5, 0.5, 0.5, 0.5, 0.5;0.5, 0.5, 0.5, 0.5, 0.5;0.5, 0.5, 0.5, 0.5, 0.5;0.5, 0.5, 0.5, 0.5, 0.5;0.5, 0.5, 0.5, 0.5, 0.5];

lambda0 = [0.5,0.5,0.5,0.5,0.5];

rate = zeros(5,5);

[N,X,E,M,E_source,lambda,spectra,cross_spectra,coherences] = MvHawkesSimulation(v,A,B,lambda0);
data = X;
for n = 1:3
        [N,X,E,M,E_source,lambda,spectra,cross_spectra,coherences] = MvHawkesSimulation(v,A,B,lambda0);
        data(:,:,(n+1)) = X;
end

  for m = 1:5
      for n = m+1:5  
          string=zeros(1,4);
          for i = 1:4
          X=data(m,:,i);
          Y=data(n,:,i);
          [x,y]=TWCOH_quick(X,Y,1,24,1,1);
                xx=(x(~isnan(x)));
                f=0;
                for ii = 1:length(xx)
                percentile=fzero(@(x) estcoh(x,y(ii),2,0.95,0), [0.000001,0.999999]);
                if (xx(ii)>percentile || xx(ii)==percentile)
                    f=f+1;
                end
                end
                string(i)=f/length(xx)>0.2902 ;
           end
            rate(m,n) = mean(string);
            [m n]
        end
  end
  
 gamma_mat = rate;

    for i = 1:5
        for j = 1:i
            if i==j
                gamma_mat(i,j)=1;
            else
                gamma_mat(i,j) = gamma_mat(j,i);
            end
        end
    end
    
    figure()
    heatmap(gamma_mat)