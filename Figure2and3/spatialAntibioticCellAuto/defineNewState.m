function state = defineNewState(m,n)
    if m < n
        state.initialBacteria = m;
        state.maxBacteria = n;
        state.spaceSize = 100;

        state.initialSugar = 250;
        state.initialDrug = 30;
        state.drugDegradation = 0.001;
        state.mutationSize = 1;

        state.drugDiffusionRate = 0.1;
        
        state.drugReleaseRate = 0.185;        %should be less than 1
        state.maxTransporters = 5;
        if state.drugReleaseRate*state.maxTransporters <= 0.5
            error('Can''t use these parameters or else drug never unbinds from transporters')
        end
        
        state.maxDrugBindingHalfSat = 6;
        
        %should be less than 1
        state.metabolicRate = 0.2;
        
        state.divisionATP = 10;
        state.antibioticDeathThreshold = 5.5;
        
        state.resourceBindingHalfSat = 10;

        state.bacteriaX = ceil(state.spaceSize*rand(n,1));
        state.bacteriaY = ceil(state.spaceSize*rand(n,1));

        state.internalSugar = ones(n,1);
        state.internalDrug = zeros(n,1);

        state.freeTransporters = state.maxTransporters*ones(n,1);    
        state.drugBoundTransporters = zeros(n,1);
        
        state.drugBindingHalfSat = ones(n,1)/10;

        state.ATP = (3/2)*state.antibioticDeathThreshold*ones(n,1);

        state.sugar = state.initialSugar*ones(state.spaceSize,state.spaceSize);
        
        state.drug = zeros(state.spaceSize,state.spaceSize);        
        state.drugLocationsX = floor(state.spaceSize/3):state.spaceSize;
        state.drugLocationsY = 1:state.spaceSize;
        state.drug(state.drugLocationsX,state.drugLocationsY) = state.initialDrug;

        [~,~,L]=laplacian([state.spaceSize,state.spaceSize],{'NN','NN'});

        state.DiffusionCoefficient = 0.3;
        state.DiffusionMatrix = myExp(-L*state.DiffusionCoefficient,5);

        state.aliveBacteria = m;
        
        state.maxEffluxPumps = 40;
        state.maxEffluxPumpsATPCost = 0.000;
        state.EffluxPumps = round(state.maxEffluxPumps*rand(n,1));
        state.boundEffluxPumps = zeros(n,1);
        state.EffluxDrugHalfSat = 0.01;
        %a number less than one: closer to zero means pump no. is more stable
        state.pumpRegulationRate = 0.5;
    else
        error('More bacteria are being used than can be accounted for')
    end
end

function ExpM = myExp(A,N)
    I = sparse(eye(size(A)));
    A = sparse(A);
    Exp = I;
    for i = N:-1:1 %a for loop with a **descending** argument
        Exp = I + Exp*A/i;
    end
    ExpM = Exp;
end