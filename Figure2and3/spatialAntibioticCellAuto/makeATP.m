function outState = makeATP(inState)
    outState = inState;
    N = inState.aliveBacteria;
    
    metaboliseAmount = inState.metabolicRate*rand(N,1);
    sugarLoss = outState.internalSugar(1:N).*metaboliseAmount(1:N);
    
    outState.internalSugar(1:N) = outState.internalSugar(1:N) - sugarLoss;
    outState.ATP(1:N) = outState.ATP(1:N) + sugarLoss;
end