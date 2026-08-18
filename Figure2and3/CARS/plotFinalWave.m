function plotFinalWave(mySSurface,myRSurface,myASurface,myCSurface,...
    parameters,label)

    if nargin < 4
        label = '';
    end

    clc
    X = -(parameters.X);
    [N,~] = size(mySSurface); 
	thickness = 4;

    figure();
    set(gcf,'Position',[1200         813         569         393])

    j = N;
    S = mySSurface(j,:);
    R = myRSurface(j,:);
    C = myCSurface(j,:);
    C = C/max(C);
    A = myASurface(j,:);
    A = A/max(A);
    plot(X,S,'-','linewidth',thickness);
    hold on
    plot(X,R,'-','linewidth',thickness);
    plot(X,A,'--k','linewidth',3);
    plot(X,C,'-','color',[0 0.75 0],'linewidth',2);

    %area(X,S,'FaceColor',[0,1,0]);
    %alpha(0.2)

    set(gca,'Xtick',[-1 -0.5 0])
    set(gca,'Xticklabel',{'drug source','spatial coord.','least drug'})
    legend({'final S density','final R density','final relative antibiotic concentration','final relative nutrient concentration (99.9% consumed)'})
	legend('location','northeast')
    axis tight
    ax = axis;

    axis([min(X) max(X) 0 ax(4)])
    ylabel('bacterial density')

    export_fig(['figures/finalSRWave-',date,'-',label,'.pdf'])            

end