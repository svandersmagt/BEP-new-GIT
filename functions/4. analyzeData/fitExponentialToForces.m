function fit = fitExponentialToForces(bead, zmags, configVariable)
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
        plot(zmags,bead(1).DaldropForceLong,'b.');
        hold on
        plot(zmags,bead(1).DaldropForceShort,'r.');
        plot(zmags,bead(1).ZhongboForceShort,'g.');
        title('Force on the magnetic bead vs magnet height');
        xlabel('magnet height (mm)');
        ylabel('Force (pN)');
        legend('Daldrops method long pendulum direction','Daldrops method short pendulum direction','Zhongbos method short pendulum direction');
    end
    
    %%%Defining the exponential function
    exponential = @(delta, alpha0, zeta0, alpha1, zeta1, magnetHeight)...
        (delta + alpha0.*exp(-(magnetHeight)./zeta0) + alpha1.*exp(-(magnetHeight)./zeta1));
    options = optimoptions('lsqnonlin','MaxFunEvals',10000,'MaxIter',10000);
    
    %%%Fitting Daldrop force long pendulum direction
    display('Daldrop long pendulum direction');
    [parLong] = lsqnonlin(@(par) exponential(par(1),par(2),par(3),par(4),par(5),zmags)...
        - bead(1).DaldropForceLong,[0,50,1,50,1],[],[],options);
    fit.deltaLong = parLong(1); fit.alpha0Long = parLong(2); fit.zeta0Long = parLong(3); 
    fit.alpha1Long = parLong(4); fit.zeta1Long = parLong(5);
    
    %%%Fitting Daldrop force short pendulum direction
    display('Daldrop short pendulum direction');
    [parShort] = lsqnonlin(@(par) exponential(par(1),par(2),par(3),par(4),par(5),zmags)...
        - bead(1).DaldropForceShort,[0,50,1,50,1],[],[],options);
    fit.deltaShort = parShort(1); fit.alpha0Short = parShort(2); fit.zeta0Short = parShort(3); 
    fit.alpha1Short = parShort(4); fit.zeta1Short = parShort(5);
    
    %%%Fitting Zhongbo force short pendulum direction
    display('Zhongbo short pendulum direction');
    [parZh] = lsqnonlin(@(par) exponential(par(1),par(2),par(3),par(4),par(5),zmags)...
        - bead(1).ZhongboForceShort,[0,50,1,50,1],[],[],options);
    fit.deltaZh = parZh(1); fit.alpha0Zh = parZh(2); fit.zeta0Zh = parZh(3); 
    fit.alpha1Zh = parZh(4); fit.zeta1Zh = parZh(5);
    
    exponentialLong = exponential(fit.deltaLong, fit.alpha0Long, fit.zeta0Long, fit.alpha1Long,...
        fit.zeta1Long, zmags);
    exponentialShort = exponential(fit.deltaShort, fit.alpha0Short, fit.zeta0Short, fit.alpha1Short,...
        fit.zeta1Short, zmags);
    exponentialZh = exponential(fit.deltaZh, fit.alpha0Zh, fit.zeta0Zh, fit.alpha1Zh,...
        fit.zeta1Zh, zmags);
    
    
    if plotThings;
        plot(zmags,exponentialLong,'b');
        plot(zmags,exponentialShort,'r');
        plot(zmags,exponentialZh,'g');
        hold off
    end
end