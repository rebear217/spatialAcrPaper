function v = fR(u,globalParameters)

    N = globalParameters.N;
    mu = globalParameters.mu;
    V = globalParameters.V;
    K = globalParameters.K;
    r = globalParameters.r;
    ab = globalParameters.abBinding;
    abD = globalParameters.abDecay;
    rK = globalParameters.reducedKill;
    killing = globalParameters.killing;
    
    diffusionRate = globalParameters.diffusionRate;
    ABdiffusionRate = globalParameters.ABdiffusionRate;
    yield = globalParameters.yield;

    S = u(1:N);
    R = u(N+1:2*N);
    R1 = u(2*N+1:3*N);
    A = u(3*N+1:4*N);
    C = u(4*N+1:5*N);
    
    uptakeRate = V*C./(C+K);
    gRate = yield*uptakeRate;
    
    dC = ABdiffusionRate*D(C) - uptakeRate.*(S+R);
    dA = ABdiffusionRate*D(A) - A*abD - (S+R).*A*ab;
    dR = r*R.*gRate - rK*killing*A.*R + mu*S + diffusionRate*D(R) - mu*R;
    dR1 = r*R.*gRate/2 - rK*killing*A.*R/2 + mu*R + diffusionRate*D(R1);
    dS = S.*gRate - killing*A.*S + diffusionRate*D(S) - mu*S;
    
    v = [dS ; dR ; dR1 ; dA; dC];
    
end
