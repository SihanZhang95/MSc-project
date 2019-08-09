a=linspace(0.1111,1,9);
plot(a,k3,'-o','MarkerIndices',1:1:length(a));
hold on 
plot(a,k6,'-o','MarkerIndices',1:1:length(a));
hold on 
plot(a,k12,'-o','MarkerIndices',1:1:length(a));
hold on 
plot(a,k24,'-o','MarkerIndices',1:1:length(a));
legend({'kappa=3','kappa=6','kappa=12','kappa=24'},'FontSize',12,'Location','southeast')
hold off