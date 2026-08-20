function [iC,globalParameters] = setInitialDataR(globalParameters,option)

    if nargin < 2
        option = 1;
    end

    N = globalParameters.N;
    antibioticPos = globalParameters.antibioticPos;

    C = ones(N,1)*globalParameters.C0;
    A = zeros(N,1);
    R = zeros(N,1);
    R1 = zeros(N,1);
    S = ones(N,1);
    
    M = floor(antibioticPos*N);

    A(M:end) = globalParameters.A0;
    
    switch option
        case 1
            S = R;
            S(1:10) = globalParameters.S0/10;
            R(1:10) = S(1:10)/1000;
            globalParameters.abDecay = 0.01;
        case 2
            S = globalParameters.S0*S;
        case 3
            S = globalParameters.S0*S;
            R = S;
        case 4
            S = R;
            S(1:floor(9*N/10)) = globalParameters.S0/10;
            globalParameters.abDecay = 0.02;          
            %globalParameters.mu = 0;
            %globalParameters.diffusionRate = globalParameters.diffusionRate * 2;
        case 5
            S = R;
            S(1:floor(9*N/10)) = globalParameters.S0/20;
            globalParameters.abDecay = 0.04;     
            globalParameters.mu = 0;
        case 6
            %globalParameters.S0 = globalParameters.S0 * 2;
            S = R;
            S(1:floor(90*N/100)) = globalParameters.S0/100;
            globalParameters.abDecay = 0.03;
            globalParameters.killing = globalParameters.killing/3;
            globalParameters.mu = globalParameters.mu / 10;
            globalParameters.diffusionRate = globalParameters.diffusionRate / 20;
        case 7
            %globalParameters.S0 = globalParameters.S0 * 2;
            S = R;
            S(1:floor(90*N/100)) = globalParameters.S0/100;
            globalParameters.abDecay = 0.03;
            globalParameters.killing = globalParameters.killing/3;
            globalParameters.mu = 0;
            globalParameters.diffusionRate = globalParameters.diffusionRate / 20;
        case 8
            %globalParameters.S0 = globalParameters.S0 * 2;
            S = R;
            S(1:floor(90*N/100)) = globalParameters.S0/100;
            globalParameters.abDecay = 0.01;
            globalParameters.killing = globalParameters.killing/5;
            globalParameters.mu = 0;
            %globalParameters.diffusionRate = globalParameters.diffusionRate / 20;
        otherwise
            S = globalParameters.S0*S;
            A = 0*A;
            R = 0*R;
            S(M:end) = 0;
    end
    
    iC = [S ; R ; R1 ; A ; C];
    
end