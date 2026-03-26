function outputs = summaryPlot(epoch,state,totalCells)

colormap(cool)

subplot(1,10,1)
    semilogy([0 ; (1:length(totalCells))'],[state.initialBacteria ; totalCells],'-b');
    xlabel(['time (epoch ',num2str(epoch),')'])
    ylabel('log #cells')
    axis tight

subplot(1,10,2)
    hist(state.EffluxPumps);
    ylabel('#cells')
    xlabel('#flux pumps')

subplot(1,10,3)
    surf(state.sugar,'EdgeColor','none','FaceColor','interp');
    xlim([1 state.spaceSize]);
    ylim([1 state.spaceSize]);
    title('sugar')
    zlim([0 state.initialSugar])
    ylabel('r')
    xlabel('theta')    

subplot(1,10,4)
    surf(state.drug,'EdgeColor','none','FaceColor','interp');
    xlim([1 state.spaceSize]);
    ylim([1 state.spaceSize]);
    title('drug')
    zlim([0 state.initialDrug])
    ylabel('r')
    xlabel('theta')

subplot(1,10,5)
    [spaceBugMatrix,spaceDrugMatrix,spaceSugarMatrix,pumpMatrix] = makeSpatialBacterialDistribution(state);
    surf(spaceBugMatrix,'EdgeColor','none','FaceColor','interp');
    xlim([1 state.spaceSize]);
    ylim([1 state.spaceSize]);
    view(2)
    title('cell density')
    ylabel('r')
    xlabel('theta')

W = 0;
subplot(1,10,6)
    surf((W+spaceDrugMatrix)./(W+spaceBugMatrix),'EdgeColor','none','FaceColor','interp');
    title('drug/cell')
    xlim([1 state.spaceSize]);
    ylim([1 state.spaceSize]);
    zlim([0 12])
    view(2);
    grid off
    box on
    ylabel('r')
    xlabel('theta')
    colorbar

subplot(1,10,7)
    surf((W+spaceSugarMatrix)./(W+spaceBugMatrix),'EdgeColor','none','FaceColor','interp');
    title('sugar/cell')
    xlim([1 state.spaceSize]);
    ylim([1 state.spaceSize]);
    view(2);
    grid off
    box on
    ylabel('r')
    xlabel('theta')

subplot(1,10,8)
    bugs = sum(spaceBugMatrix,2);
    X = state.spaceSize:-1:1;
    plot(X,bugs,'Color','b','Linewidth',2);
    xlim([1 state.spaceSize])
    ylabel('bacterial density')
    xlabel('distance to ab (r)')
    axis tight

subplot(1,10,9)
    plot(X,sum(spaceDrugMatrix,2),'Color','k','Linewidth',2);
    xlim([1 state.spaceSize])
    ylabel('intracell drug')
    xlabel('distance to ab (r)')
    axis tight

subplot(1,10,10)
    PumpsPerCell = sum(pumpMatrix,2)./bugs;
    plot(X,PumpsPerCell,'Color','r','Linewidth',2,'DisplayName','#pumps/cell')
    ylabel('pumps/cell')
    xlabel('distance to ab (r)')
    xlim([1 state.spaceSize])    
    axis tight

drawnow

outputs.X = state.spaceSize:-1:1;
outputs.Y = bugs;

end