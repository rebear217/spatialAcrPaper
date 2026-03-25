function [state,plotStates] = simulation(Tmax)

    close all
    clc

    nCells = 2000;
    nCellsMax = 1000000;
    timeMax = 20*Tmax;
    
    state = defineNewState(nCells,nCellsMax);    
    totalCells = zeros(timeMax,1);
    
    figure(1);
    set(1,'position',[64         110        1374         823])
    set(1,'color','white')
    
    figure(2);
    set(2,'position',[147        1003         951         366])
    set(2,'color','white')
    
    figure(3);
    set(3,'position',[1157        1027        1147         342])
    set(3,'color','white')

    figure(4);
    set(4,'position',[1445         109         595         824])
    set(4,'color','white')
    
    plotStates = zeros(timeMax/20,100);
    count = 1;
    for j = 1:timeMax
        state = uptakeSugarAndDrug(state);
        state = makeATP(state);
        state = diffuseSugar(state);
        state = diffuseDrug(state);
        state = killCells(state);
        state = divideCells(state);
        state = degradeDrug(state);
        %state = replenishDrug(state);
        
        totalCells(j) = state.aliveBacteria;
        if mod(j,20) == 0
            outState = plotState(state,totalCells(1:j));
            plotStates(count,:) = outState.Y;
            plotHisto(state);
            count = count + 1;
        end

    end

    figure(4)
    subplot(3,1,3)
    [spaceBugMatrix,~,~,pumpMatrix] = makeSpatialBacterialDistribution(state);
    bugs = sum(spaceBugMatrix,2);
    X = state.spaceSize:-1:1;
    PumpsPerCell = sum(pumpMatrix,2)./bugs;
    LM = fitlm(X,PumpsPerCell);
    hold on
    plot(X,LM.feval(X),'-k','DisplayName','linear regression');
    legend('location','northwest')

end

function plotHisto(state)
    figure(2);
    subplot(1,2,1);
    hist(state.drugBindingHalfSat);
    ylabel('numbers of cells')
    xlabel('drug-binding half saturation constant')
    subplot(1,2,2);
    hist(state.EffluxPumps);
    ylabel('numbers of cells')
    xlabel('no of efflux pumps')
    drawnow;
end



function outputs = plotState(state,totalCells)

    figure(1);
    
    subplot(2,2,1)
    semilogy([0 ; (1:length(totalCells))'],[state.initialBacteria ; totalCells],'-b');
    xlabel('time')
    ylabel('total cell count (log scale)')
    
    subplot(2,2,2)
    surf(state.sugar,'EdgeColor','none');
    xlim([1 state.spaceSize]);
    ylim([1 state.spaceSize]);
    %view(2);
    title('spatial sugar distribution')
    zlim([0 state.initialSugar])
    ylabel('spatial coord (drug source at top)')
    xlabel('spatial coord')    

    subplot(2,2,3)
    surf(state.drug,'EdgeColor','none');
    xlim([1 state.spaceSize]);
    ylim([1 state.spaceSize]);
    %view(2);
    title('spatial drug distribution')
    zlim([0 state.initialDrug])
    ylabel('spatial coord (drug source at top)')
    xlabel('spatial coord')

    [spaceBugMatrix,spaceDrugMatrix,spaceSugarMatrix,pumpMatrix] = makeSpatialBacterialDistribution(state);
    subplot(2,2,4)
    surf(spaceBugMatrix,'EdgeColor','none');
    xlim([1 state.spaceSize]);
    ylim([1 state.spaceSize]);
    view(2);
    title('spatial cell distribution')
    %colorbar
    ylabel('spatial coord (drug source at top)')
    xlabel('spatial coord')
    
    figure(3);
    subplot(1,2,1)
    
    surf(-1+(1+spaceDrugMatrix)./(1+spaceBugMatrix),'EdgeColor','none');
    title('drug molecules per cell')
    xlim([1 state.spaceSize]);
    ylim([1 state.spaceSize]);
    view(2);
    grid off
    box on
    ylabel('spatial coord (drug source at top)')
    xlabel('spatial coord')

    subplot(1,2,2)
    surf(-1+(1+spaceSugarMatrix)./(1+spaceBugMatrix),'EdgeColor','none');
    title('sugar molecules per cell')
    xlim([1 state.spaceSize]);
    ylim([1 state.spaceSize]);
    view(2);
    grid off
    box on
    ylabel('spatial coord (drug source at top)')
    xlabel('spatial coord')
    
    figure(4)
    bugs = sum(spaceBugMatrix,2);
    subplot(3,1,1)
    X = state.spaceSize:-1:1;
    plot(X,bugs,'Color','b','Linewidth',2);
    xlim([1 state.spaceSize])
    ylabel('total bacterial density')
    outputs.X = state.spaceSize:-1:1;
    outputs.Y = bugs;
    xlabel('spatial location: distance from highest drug concentration')

    subplot(3,1,2)
    plot(X,sum(spaceDrugMatrix,2),'Color','k','Linewidth',2);
    xlim([1 state.spaceSize])
    ylabel('total drug (within-cell)')
    xlabel('spatial location: distance from highest drug concentration')

    subplot(3,1,3)
    PumpsPerCell = sum(pumpMatrix,2)./bugs;
    plot(X,PumpsPerCell,'Color','r','Linewidth',2,'DisplayName','pumps per cell')
    ylabel('efflux pumps per cell')
    xlabel('spatial location: distance from highest drug concentration')
    xlim([1 state.spaceSize])
   
    drawnow;
end

