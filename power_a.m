function [power,sig] = power_a(n,m,N,C,ka,d,delta)
  z=linspace(0,1,m);
  power_sig=zeros(m,1);
  for i = 1:m
    coef=z(i);
    string=zeros(n,1);
    for j = 1:n
    V = sim_ar2(5/6,-1/6,N);
    W = sim_ar2(3/4,-1/8,N);
    Z = sim_ar2(7/12,-1/12,N);
    X = coef*V + (1-coef)*Z;
    Y = coef*V + (1-coef)*W;
    [x,y]=TWCOH_quick(X,Y,delta,ka,d,1);
    xx=(x(~isnan(x)));
    f=0;
      for jj = 1:length(xx)
          percentile=fzero(@(x) estcoh(x,y(jj),2,0.95,0), [0.000001,0.999999]);
          if (xx(jj)>percentile || xx(jj)==percentile)
                    f=f+1;
          end
      end
      string(j)=f/length(xx)>C;
    end
    power_sig(i)=mean(string);
  i;
  end
 sig=power_sig(1);
 power=power_sig(2:m);
end