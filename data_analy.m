Delta = 0.28;
alpha = 0.05;
kappa=3;
C=0.0627;

for sfi = 7
    
    rate = zeros(185,185);
    
    for m = 1:185
        for n = m+1:185
            string=zeros(1,4);
            for i = 1:4
                sfc1 = get_tf(cells, m,i+4*(sfi-1));
                proc1_freq1 = histcounts(sfc1,1:Delta:15);
                sf1 = transpose(proc1_freq1-mean(proc1_freq1));

                sfc2 = get_tf(cells,n,i+4*(sfi-1));
                proc1_freq2 = histcounts(sfc2,1:Delta:15);
                sf2 = transpose(proc1_freq2-mean(proc1_freq2));
                [x,y]=TWCOH_quick(sf1.',sf2.',1,3,1,1);
                xx=(x(~isnan(x)));
                f=0;
                for ii = 1:length(xx)
                percentile=fzero(@(x) estcoh(x,y(ii),2,0.95,0), [0.000001,0.999999]);
                if (xx(i)>percentile || xx(i)==percentile)
                    f=f+1;
                end
                end
                string(i)=f/length(xx)>C ;
            end
            % rejection rate when goes through the stepdown procedure
            rate(m,n) = mean(string);
            [m n]
        end
    end

    gamma_mat = rate;

    for i = 1:185
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
    
end