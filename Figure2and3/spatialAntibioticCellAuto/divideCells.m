function outState = divideCells(inState)
    outState = inState;
    dividingNow = inState.ATP > inState.divisionATP;
    NdividingNow = sum(dividingNow);
    if NdividingNow > 0
        parentcell = dividingNow;
        daughtercell = outState.aliveBacteria + (1:NdividingNow);

        Xshift = floor(3*rand(NdividingNow,1))-1;
        Yshift = floor(3*rand(NdividingNow,1))-1;

        outState.bacteriaX(daughtercell) = outState.bacteriaX(parentcell)+Xshift;
        outState.bacteriaY(daughtercell) = outState.bacteriaY(parentcell)+Yshift;

        outState.ATP(daughtercell) = outState.ATP(parentcell)/2;
        outState.ATP(parentcell) = outState.ATP(parentcell)/2;

        outState.internalSugar(daughtercell) = outState.internalSugar(parentcell)/2;
        outState.internalSugar(parentcell) = outState.internalSugar(parentcell)/2;

        outState.freeTransporters(daughtercell) = outState.freeTransporters(parentcell);
        outState.drugBoundTransporters(daughtercell) = outState.drugBoundTransporters(parentcell);
        
        outState.EffluxPumps(daughtercell) = outState.EffluxPumps(parentcell);
        outState.boundEffluxPumps(daughtercell) = outState.boundEffluxPumps(parentcell);
        
        if outState.mutationSize > 0
            %drugBindingHalfSat IS NOT kA
            %drugBindingHalfSat IS kAT
            
            pH = (1+outState.mutationSize*(rand(NdividingNow,1)-0.5)).*outState.drugBindingHalfSat(parentcell);
            pH = 0*(pH < 0) + pH.*(pH >= 0).*(pH <= outState.maxDrugBindingHalfSat) + ...
                outState.maxDrugBindingHalfSat*(pH > outState.maxDrugBindingHalfSat);
            outState.drugBindingHalfSat(daughtercell) = pH;
        else
            outState.drugBindingHalfSat(daughtercell) = outState.drugBindingHalfSat(parentcell);
        end

        outState.internalDrug(daughtercell) = outState.internalDrug(parentcell)/2;
        outState.internalDrug(parentcell) = outState.internalDrug(parentcell)/2;
        
        outState.aliveBacteria = outState.aliveBacteria + sum(dividingNow);

        bX = outState.bacteriaX;
        bX = 1*(bX < 1) + bX.*(bX >= 1).*(bX <= outState.spaceSize) + outState.spaceSize*(bX > outState.spaceSize);
        outState.bacteriaX = bX;

        bY = outState.bacteriaY;
        bY = 1*(bY < 1) + bY.*(bY >= 1).*(bY <= outState.spaceSize) + outState.spaceSize*(bY > outState.spaceSize);
        outState.bacteriaY = bY;
    end
    
end
