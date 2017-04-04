function [fit] = fitExponentialToForces(bead, zmags, configVariable)
    
    plotThings = configVariable.plotThings;
    %%% Plot forces;
    if plotThings;
        figure(7);
        plot(zmags,bead(1).DaldropXForce,'b.');
        hold on
        plot(zmags,bead(1).DaldropYForce,'r.');
        plot(zmags,bead(1).ZhongboYForce,'g.');
        title('Force on the magnetic bead vs magnet height');
        xlabel('magnet height (mm)');
        ylabel('Force (pN)');
        legend('Daldrops method x-direction','Daldrops method y-direction','Zhongbos method y-direction');
    end
    
    %%%Defining the exponential function
    exponential = @(delta, alpha0, zeta0, alpha1, zeta1, magnetHeight)...
        (delta + alpha0.*exp(-(magnetHeight)./zeta0) + alpha1.*exp(-(magnetHeight)./zeta1));
    options = optimoptions('lsqnonlin','TolFun',1E-10,'MaxFunEvals',10000);
    
    %%%Fitting Daldrop force x-direction
    [parX] = lsqnonlin(@(par) exponential(par(1),par(2),par(3),par(4),par(5),zmags)...
        - bead(1).DaldropXForce,[0,50,1,50,1],[],[],options);
    fit.deltaX = parX(1); fit.alpha0X = parX(2); fit.zeta0X = parX(3); 
    fit.alpha1X = parX(4); fit.zeta1X = parX(5);
    
    %%%Fitting Daldrop force y-direction
    [parY] = lsqnonlin(@(par) exponential(par(1),par(2),par(3),par(4),par(5),zmags)...
        - bead(1).DaldropYForce,[0,50,1,50,1],[],[],options);
    fit.deltaY = parY(1); fit.alpha0Y = parY(2); fit.zeta0Y = parY(3); 
    fit.alpha1Y = parY(4); fit.zeta1Y = parY(5);
    
    %%%Fitting Daldrop force y-direction
    [parZh] = lsqnonlin(@(par) exponential(par(1),par(2),par(3),par(4),par(5),zmags)...
        - bead(1).ZhongboYForce,[0,50,1,50,1],[],[],options);
    fit.deltaZh = parZh(1); fit.alpha0Zh = parZh(2); fit.zeta0Zh = parZh(3); 
    fit.alpha1Zh = parZh(4); fit.zeta1Zh = parZh(5);
    
    exponentialX = exponential(fit.deltaX, fit.alpha0X, fit.zeta0X, fit.alpha1X,...
        fit.zeta1X, zmags);
    exponentialY = exponential(fit.deltaY, fit.alpha0Y, fit.zeta0Y, fit.alpha1Y,...
        fit.zeta1Y, zmags);
    exponentialZh = exponential(fit.deltaZh, fit.alpha0Zh, fit.zeta0Zh, fit.alpha1Zh,...
        fit.zeta1Zh, zmags);
    
    
    if plotThings;
        plot(zmags,exponentialX,'b');
        plot(zmags,exponentialY,'r');
        plot(zmags,exponentialZh,'g');
        hold off
    end