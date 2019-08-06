function theoretical_values = Goodman_QQ_Plots(true_value,observed_values,kk)

L = length(observed_values);

pp=(1:L)./(L+1);

theoretical_values = zeros(1,L);

for j=1:L
    prob=pp(j);
    theoretical_values(j) = fzero(@(x) estcoh(x,kk,2,prob,true_value),[0.000001,0.999999]);
end