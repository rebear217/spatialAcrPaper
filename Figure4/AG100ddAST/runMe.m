clear all
close all
clc

fileName{1} = 'allDataKO';
fileName{2} = 'allDataWT';

for k = 1:2
    close all
    
    figure(1)
    set(1,'pos',[742   939   736   407]);

    load([fileName{k},'.mat'])
    spatialDoseResponsedata = y1;
    clear y1;

    v = mean(spatialDoseResponsedata);
    [M,N] = size(spatialDoseResponsedata);
    
    mv = max(v);
    sd = std(spatialDoseResponsedata,[],1)/sqrt(M-1);

    X = (0:N-1)/(N-1);    

    plots = [];
    
    plots(2) = plotshaded(X,[v-1.96*sd ; v + 1.96*sd],'r');
    hold on
    set(plots(2),'facecolor',[0.7 0.7 0.7]);
    plots(1) = plot(X,v,'-k','linewidth',1);
    box on
    legend(plots,{'mean OD','+/- 95% CI (n = 30)'},'location','northwest');
    legend('boxoff');
    axis tight
    set(gca,'Xtick',[0 1])
    set(gca,'Xticklabel',{'centre','edge'},'fontsize',14)
    drawnow
    
    export_fig(['rawODProfile',fileName{k},'.pdf'])
end

%%

colours = {[0.7 0.7 0.7],[1 0.7 0.7]};
linecolours = {[0.2 0.2 0.2],[1 0.2 0.2]};
plots = [];

for k = 1:2
    
    figure(1)
    set(1,'pos',[742   939   736   407]);

    load([fileName{k},'.mat'])
    spatialDoseResponsedata = y1;
    clear y1;

    v = mean(spatialDoseResponsedata);
    [M,N] = size(spatialDoseResponsedata);
    
    mv = max(v);
    sd = std(spatialDoseResponsedata,[],1)/sqrt(M-1);

    X = (0:N-1)/(N-1);    
    
    p = plotshaded(X,[v-1.96*sd ; v + 1.96*sd],'r');
    plots = [plots p];
    hold on
    set(p,'facecolor',colours{k});
    p = plot(X,v,'-k','linewidth',1,'color',linecolours{k});
    plots = [plots p];
	box on
    axis tight
    set(gca,'Xtick',[0 1])
    set(gca,'Xticklabel',{'centre','edge'},'fontsize',14)
        
end

legend(plots,{'AG100A +/- 95% CI (n = 30)','AG100A mean OD','AG100 +/- 95% CI (n = 30)','AG100 mean OD'},'location','northwest');
legend('boxoff');


export_fig('rawODProfileOverlay.pdf')
