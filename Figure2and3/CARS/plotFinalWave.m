function plotFinalWave(mySSurface,myRSurface,parameters,label)

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
    plot(X,S,'-','linewidth',thickness);
    hold on
    plot(X,R,'-','linewidth',thickness);
    %area(X,S,'FaceColor',[0,1,0]);
    %alpha(0.2)
    set(gca,'Xtick',[-1 -0.5 0])
    set(gca,'Xticklabel',{'drug source','spatial coord.','least drug'})
    legend({'final S density','final R density'})
	legend('location','southeast')
    axis tight
    ax = axis;

    axis([min(X) max(X) 0 ax(4)])
    ylabel('bacterial density')

    export_fig(['figures/finalSRWave-',date,'-',label,'.pdf'])            

end