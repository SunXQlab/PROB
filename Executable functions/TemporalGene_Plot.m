function TemporalGene_Plot(geneID,TimeSampled,Data_smooth)
 % geneID =  the ID of genes you want to plot.
 % TimeSampled = standardized time-points sampled for Data_smooth.
 % Data_smooth = ordered and smoothed gene expression data along with pseudotemporal progression trajectory. 
figure,
for i=geneID
hold on, plot(TimeSampled,Data_smooth(i,:),'-','color',[rand rand rand],'LineWidth',3)
end
set(gca,'FontSize',15)
xlabel('Inferred pseudo-progression','FontSize',20);
ylabel('Gene expression','FontSize',20)
Leg=legend(num2str(geneID'),'Location','NorthEastOutside','Box','off')
set(Leg,'Box','off','FontSize',12)
% axis([0 1 -1 3])