function [proportion] = C(n,N,ka,d,delta)
  proportion=zeros(1,n);
  for i = 1:n
    X = normrnd(0,1,[1,N]);
    Y = normrnd(0,sqrt(10),[1,N]);
    [x,y]=TWCOH_quick(X,Y,delta,ka,d,1);
    xx=(x(~isnan(x)));
    f=0;
      for jj = 1:length(xx)
          percentile=fzero(@(x) estcoh(x,y(jj),2,0.95,0), [0.000001,0.999999]);
          if (xx(jj)>percentile || xx(jj)==percentile)
                    f=f+1;
          end
      end
      proportion(i)=f/length(xx);
  end
end