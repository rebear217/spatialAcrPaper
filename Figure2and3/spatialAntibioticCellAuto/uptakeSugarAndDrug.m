function outState = uptakeSugarAndDrug(inState)

    N = inState.aliveBacteria;
    reorder = randperm(N);
    
    resourceBindingHalfSat = inState.resourceBindingHalfSat;
    drugDiffusionRate  = inState.drugDiffusionRate;
    
    SugarReleaseProbability = rand(N,1);
    SugarBindProbability = rand(N,1);
    
    TransporterReleaseProbability = inState.drugReleaseRate.*rand(N,1);
    TransporterBindProbability = rand(N,1);
    
    pumpShift = floor(3*rand(N,1))-1;
    bindingPumpsCoeff = rand(N,1);
    pumpchangeCoeff = rand(N,1);

    for j = 1:N
        cell = reorder(j);
        
        X = inState.bacteriaX(cell);
        Y = inState.bacteriaY(cell);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %take up sugar
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        externalSugar = inState.sugar(X,Y);
        internalSugar = inState.internalSugar(cell);
        
        freeTransporters = inState.freeTransporters(cell);
        boundTransporters = inState.maxTransporters - freeTransporters - inState.drugBoundTransporters(cell);
        
        transportingNow = round(SugarReleaseProbability(cell)*boundTransporters);
        internalSugar = internalSugar + transportingNow;
        
        freeTransporters = freeTransporters + transportingNow;
        %boundTransporters = boundTransporters - transportingNow;
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        %resourceBindingHalfSat = kS
        bindingProbability = externalSugar/(resourceBindingHalfSat + externalSugar);
        bindingTransportersNow = round(SugarBindProbability(cell) * freeTransporters * bindingProbability);
        
        if bindingTransportersNow > 0
            if externalSugar > bindingTransportersNow
                freeTransporters = freeTransporters - bindingTransportersNow;
                %boundTransporters = boundTransporters + bindingTransportersNow;

                externalSugar = externalSugar - bindingTransportersNow;
            else
                freeTransporters = freeTransporters - externalSugar;
                %boundTransporters = boundTransporters + externalSugar;
                externalSugar = 0;            
            end
        end
        
        inState.sugar(X,Y) = externalSugar;
        inState.internalSugar(cell) = internalSugar;
        inState.freeTransporters(cell) = freeTransporters;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %take up drug
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        externalDrug = inState.drug(X,Y);
        internalDrug = inState.internalDrug(cell);
        
        deltaDrug = drugDiffusionRate * (internalDrug - externalDrug);
        internalDrug = internalDrug - deltaDrug;
        externalDrug = externalDrug + deltaDrug;
        
        inState.drug(X,Y) = externalDrug;
        inState.internalDrug(cell) = internalDrug;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %bind drug to transporters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        drugBindingHalfSat = inState.drugBindingHalfSat(cell);
        drug = inState.internalDrug(cell);
        freeTransporters = inState.freeTransporters(cell);
        drugBoundTransporters = inState.drugBoundTransporters(cell);
        
        %kA IS NOT drugBindingHalfSat 
        bindingProbability = drug/(drugBindingHalfSat + drug);
        bindingTransportersNow = round(TransporterBindProbability(cell) * freeTransporters * bindingProbability);
        
        if bindingTransportersNow > 0
            if freeTransporters > bindingTransportersNow
                freeTransporters = freeTransporters - bindingTransportersNow;
                drugBoundTransporters = drugBoundTransporters + bindingTransportersNow;
            else
                drugBoundTransporters = drugBoundTransporters + freeTransporters;
                freeTransporters = 0;
            end
            drug = drug - bindingTransportersNow;
            if drug < 0
                drug = 0;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %the following line leads to problems if state.drugReleaseRate*state.maxTransporters < 0.5
        %as then drug is never released from a drug-bound sugar transporter
        releaseBoundTransportersNow = round(drugBoundTransporters*TransporterReleaseProbability(cell));
        %the following is one fix:
        %%%%%%%%%%%%%%%%%%%%%%%%
        %releaseBoundTransportersNow = ceil(drugBoundTransporters*TransporterReleaseProbability(cell));
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        drugBoundTransporters = drugBoundTransporters - releaseBoundTransportersNow;
        freeTransporters = freeTransporters + releaseBoundTransportersNow;
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        inState.internalDrug(cell) = drug;
        inState.freeTransporters(cell) = freeTransporters;
        inState.drugBoundTransporters(cell) = drugBoundTransporters;
        mEP = inState.maxEffluxPumps;
        ATPCost = inState.maxEffluxPumpsATPCost;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %efflux drug
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        drug = inState.internalDrug(cell);
        externalDrug = inState.drug(X,Y);
        
        pumps = inState.EffluxPumps(cell);
        boundPumps = inState.boundEffluxPumps(cell);
        unboundPumps = pumps - boundPumps;
        pumpRegulationRate = inState.pumpRegulationRate;
        
        %kA = state.EffluxDrugHalfSat
        bindingPumpsNow = round(bindingPumpsCoeff(cell)*unboundPumps*drug/(inState.EffluxDrugHalfSat + drug));

        %pumpNow = round(pumpoutCoeff(cell)*boundPumps);
        pumpNow = boundPumps;
        
        %pumpOut and pumpNow are the same, unless there is less drug than
        %can be pumped by one pump:
        pumpOut = (pumpNow > 0)*(pumpNow*(pumpNow <= drug) + drug*(pumpNow > drug));

        externalDrug = externalDrug + pumpOut;
        drug = drug - pumpOut;

        boundPumps = boundPumps - pumpNow + bindingPumpsNow;
        
        if pumpchangeCoeff(cell) < pumpRegulationRate
            pumps = pumps + pumpShift(cell);
            if (pumps < boundPumps)
                pumps = boundPumps;
            end
            if (pumps < 0)
                pumps = 0;
            end
            if (pumps > mEP)
                pumps = mEP;
            end
            ATP = inState.ATP(cell);
            if (pumpShift(cell) > 0) && (pumps < mEP)
                ATP = ATP - ATPCost;
            end
            if ATP < 0
                ATP = 0;
            end
            inState.ATP(cell) = ATP;
        end
        
        inState.drug(X,Y) = externalDrug;
        inState.internalDrug(cell) = drug;
        inState.EffluxPumps(cell) = pumps;
        inState.boundEffluxPumps(cell) = boundPumps;
        
    end
    
    outState = inState;
    
end