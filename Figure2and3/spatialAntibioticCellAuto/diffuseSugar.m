function outState = diffuseSugar(inState)
    outState = inState;
    
    DiffusionMatrix = inState.DiffusionMatrix;        
    outSideSugar = inState.sugar;
    outSideSugar = outSideSugar(:);
    spaceSize = inState.spaceSize;
    
    outSideSugar = DiffusionMatrix*outSideSugar;
    outState.sugar = reshape(outSideSugar,spaceSize,spaceSize);
end
