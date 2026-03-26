function [state,plotStates] = simulation(Tmax,PLOTSAVEON)

    if nargin < 2
        PLOTSAVEON = 0;
    end

    close all
    clc

    nCells = 2000;
    nCellsMax = 1000000;
    timeMax = 20*Tmax;
    
    state = defineNewState(nCells,nCellsMax);    
    totalCells = zeros(timeMax,1);
    
    figure(1);
    set(1,'position',[59         906        2502         258])
    set(1,'color','white')
    
    figure(2);
    set(2,'position',[62   355   560   420])
    set(2,'color','white')    

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
        epoch = j/20;
        if mod(j,20) == 0
            figure(1)
            outState = summaryPlot(epoch,state,totalCells(1:j));
            plotStates(count,:) = outState.Y;

            figure(2)
            plotHisto(state);
            count = count + 1;
        end
        if PLOTSAVEON && mod(j,100) == 0
            figure(1)
            exportgraphics(gcf,['./figures/stateEpoch',num2str(epoch),'.pdf'],'ContentType','vector')
            %exportgraphics(gcf,['./figures/stateEpoch',num2str(epoch),'.jpg'])
        end

    end

    figure(1)
    subplot(1,10,10)
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
    subplot(1,2,1);
    hist(state.drugBindingHalfSat);
    ylabel('numbers of cells')
    xlabel('drug-binding half saturation constant')
    subplot(1,2,2);
    hist(state.EffluxPumps);
    ylabel('numbers of cells')
    xlabel('no of efflux pumps')
    drawnow
end


