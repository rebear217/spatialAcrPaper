function [spaceMatrix,spaceDrugMatrix,spaceSugarMatrix,pumpMatrix] = makeSpatialBacterialDistribution(state)
    spaceMatrix = zeros(state.spaceSize);
    spaceDrugMatrix = zeros(state.spaceSize);
    spaceSugarMatrix = zeros(state.spaceSize);
    pumpMatrix = zeros(state.spaceSize);
    
    for cell = 1:state.aliveBacteria
        X = state.bacteriaX(cell);
        Y = state.bacteriaY(cell);
        spaceMatrix(X,Y) = spaceMatrix(X,Y) + 1;
        spaceSugarMatrix(X,Y) = spaceSugarMatrix(X,Y) + state.internalSugar(cell);
        spaceDrugMatrix(X,Y) = spaceDrugMatrix(X,Y) + state.internalDrug(cell);
        pumpMatrix(X,Y) = pumpMatrix(X,Y) + state.EffluxPumps(cell);
    end
end
