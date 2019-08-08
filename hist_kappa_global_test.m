subplot(2,2,1);
edges=linspace(0,1,101);
histogram(C_k3,'BinEdges',edges);
grid on;
xlim([0, 1]);
xline(0.0581,'--r');
xlabel('kappa=3', 'FontSize', 14);
ylabel('Frequency', 'FontSize', 14);
text(0.5,450,'950th percentile=0.0581','Color','red','FontSize',14);

subplot(2,2,2);
edges=linspace(0,1,101);
histogram(C_k6,'BinEdges',edges);
grid on;
xlim([0, 1]);
ylim([0, 500]);
xline(0.1150,'--r');
xlabel('kappa=6', 'FontSize', 14);
ylabel('Frequency', 'FontSize', 14);
text(0.5,450,'950th percentile=0.1150','Color','red','FontSize',14);

subplot(2,2,3);
edges=linspace(0,1,101);
histogram(C_k12,'BinEdges',edges);
xlabel('kappa=12');
grid on;
xlim([0, 1]);
xline(0.1927,'--r');
xlabel('kappa=12', 'FontSize', 14);
ylabel('Frequency', 'FontSize', 14);
text(0.5,450,'950th percentile=0.1927','Color','red','FontSize',14);


subplot(2,2,4);
edges=linspace(0,1,101);
histogram(C_k24,'BinEdges',edges);
grid on;
xlim([0, 1]);
xline(0.2902,'--r');
xlabel('kappa=24', 'FontSize', 14);
ylabel('Frequency', 'FontSize', 14);
text(0.5,900,'950th percentile=0.2902','Color','red','FontSize',14);
