Delta = 0.025;
kappa=12;
m=6;

sumTWCOH=0;
for i = 1:m
sfc1 = get_sf(cells,108,i+6*(3-1));
sfc2 = get_sf(cells,117,i+6*(3-1));
proc1_freq1 = histcounts(sfc1,1:Delta:8);
proc1_freq2 = histcounts(sfc2,1:Delta:8);
[x,y] = TWCOH_quick(proc1_freq1,proc1_freq2,1,kappa,1,1);
sumTWCOH=sumTWCOH+x;
end

average=sumTWCOH/m;
contourf(linspace(1,8,280), linspace(3.6732,9.1667,100),average)
colormap(hot)
title('Average TWCOH','FontSize',16)
xlabel('Time (seconds)','FontSize',16) 
ylabel('Scale','FontSize',16) 
