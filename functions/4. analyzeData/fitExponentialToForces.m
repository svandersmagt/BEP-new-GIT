function fit = fitExponentialToForces(bead, zmags, configVariable, ind)
%%% Fits the found forces for different magnet heights to a double
%%% exponential function. This is not physical, but gives an empirical
%%% formula.

%%% Input: (bead, zmags, configVariable)
%%% - struct containing the time trace and position traces of the beads
%%% - vector containing the magnetheight for every plateau

%%% Output: fit
%%% - struct containing the result of the fitting
%%
    plotThings = configVariable.plotThings;
    %%% Plot forces;
    if plotThings;
        figure(7);
        plot(zmags(ind),bead(1).forceLong(ind),'b.');
        hold on
        plot(zmags(ind),bead(1).forceShort(ind),'r.');
        title('Force on the magnetic bead vs magnet height');
        xlabel('magnet height (mm)');
        ylabel('force (pN)');
        legend('Long pendulum direction','Short pendulum direction');
    end
    
    %%%Defining the exponential function
    exponential = @(delta, alpha0, zeta0, alpha1, zeta1, magnetHeight)...
        (delta + alpha0.*exp(-(magnetHeight)./zeta0) + alpha1.*exp(-(magnetHeight)./zeta1));
    options = optimoptions('lsqnonlin','MaxFunEvals',10000,'MaxIter',10000);
    
    %%%Fitting Daldrop force long pendulum direction
    display('Daldrop long pendulum direction');
    [parLong] = lsqnonlin(@(par) exponential(par(1),par(2),par(3),par(4),par(5),zmags(ind))...
        - bead(1).forceLong(ind),[0,50,1,50,1],[],[],options);
    fit.deltaLong = parLong(1); fit.alpha0Long = parLong(2); fit.zeta0Long = parLong(3); 
    fit.alpha1Long = parLong(4); fit.zeta1Long = parLong(5);
    
    %%%Fitting Daldrop force short pendulum direction
    display('Daldrop short pendulum direction');
    [parShort] = lsqnonlin(@(par) exponential(par(1),par(2),par(3),par(4),par(5),zmags(ind))...
        - bead(1).forceShort(ind),[0,50,1,50,1],[],[],options);
    fit.deltaShort = parShort(1); fit.alpha0Short = parShort(2); fit.zeta0Short = parShort(3); 
    fit.alpha1Short = parShort(4); fit.zeta1Short = parShort(5);
    
    exponentialLong = exponential(fit.deltaLong, fit.alpha0Long, fit.zeta0Long, fit.alpha1Long,...
        fit.zeta1Long, zmags(ind));
    exponentialShort = exponential(fit.deltaShort, fit.alpha0Short, fit.zeta0Short, fit.alpha1Short,...
        fit.zeta1Short, zmags(ind));
    
    if plotThings;
        plot(zmags(ind),exponentialLong,'b');
        plot(zmags(ind),exponentialShort,'r');
        hold off
    end
end