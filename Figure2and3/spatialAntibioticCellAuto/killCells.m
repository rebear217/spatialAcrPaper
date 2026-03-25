function outState = killCells(inState)
    N = inState.aliveBacteria;

    outState = inState;
    notDying = inState.ATP > inState.antibioticDeathThreshold;
	notDying = notDying(1:N);
    if not(N == sum(notDying))
        outState.aliveBacteria = sum(notDying);

        outState.bacteriaX = outState.bacteriaX(notDying);
        outState.bacteriaY = outState.bacteriaY(notDying);
        outState.ATP = outState.ATP(notDying);

        outState.internalSugar = outState.internalSugar(notDying);
        outState.internalDrug = outState.internalDrug(notDying);

        outState.freeTransporters = outState.freeTransporters(notDying);
        outState.drugBoundTransporters = outState.drugBoundTransporters(notDying);

        outState.drugBindingHalfSat = outState.drugBindingHalfSat(notDying);
        
        outState.aliveBacteria = sum(notDying);
        
        outState.EffluxPumps = outState.EffluxPumps(notDying);
        outState.boundEffluxPumps = outState.boundEffluxPumps(notDying);
    end
end
