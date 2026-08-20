function outState = solveForwardR(inState,T,parameters)

    N = parameters.N;
    options = odeset('NonNegative',ones(5*N,1));

    [~,soln] = ode23(@(t,x)fR(x,parameters),[0 T],inState,options);
    outState = soln(end,:);
    outState = outState(:);

end