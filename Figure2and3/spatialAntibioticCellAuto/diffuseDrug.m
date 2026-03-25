function outState = diffuseDrug(inState)
    outState = inState;
    
    DiffusionMatrix = inState.DiffusionMatrix;        
    outSideDrug = inState.drug;
    outSideDrug = outSideDrug(:);
    spaceSize = inState.spaceSize;
    
    outSideDrug = DiffusionMatrix*outSideDrug;
    outState.drug = reshape(outSideDrug,spaceSize,spaceSize);
end
