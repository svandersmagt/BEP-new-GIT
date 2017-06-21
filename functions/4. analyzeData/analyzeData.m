function [bead, forcesExponentialFit] = analyzeData(bead, plat, zmags, configVariable)
%%% STEP 4: ANALYZE THE BEAD DATA FOR EACH PLATEAU
%%% --------------------------------------------------------------------
%%% In this step the bead traces are analyzed. See analyzeTrace.m and
%%% fitExponentialToForces.m

%%% Input: (bead, plat, zmags, configVariable)
%%% - struct containing the time trace and position traces of the beads
%%% - struct containing the first and last index for each plateau
%%% - vector containing the magnetheight for every plateau

%%% Output: [bead, forcesExponentialFit]
%%% - result of the fits is added to the bead struct and returned
%%% - struct containing the result of the double exponential fit to the
%%%   forces vs. magnetheight
    %%
    %%% Some variables from the configuration file
    kT = configVariable.kT; %pN nm
    viscosity = configVariable.viscosity; %viscosity in pN s/nm^2
    beadRadius = configVariable.beadRadius; %bead radius in nm
    sampleFreq = configVariable.sampleFreq;
    plotThings = configVariable.plotThings;
    plotSpectrum = configVariable.plotSpectrum;
    ind = zmags <= configVariable.maxMagnetHeigth;
    plat = plat(ind);
    nPlat = length(plat);

    bead.extensionDNA = zeros(1,nPlat);
    bead.L = zeros(1,nPlat);
    bead.timeDomainForceLong = zeros(1,nPlat);
    bead.timeDomainForceShort = zeros(1,nPlat);
    bead.fitLong = zeros(1,nPlat);
    bead.forceLong = zeros(1,nPlat);
    bead.radiusLong = zeros(1,nPlat);
    bead.fitShort = zeros(1,nPlat);
    bead.forceShort = zeros(1,nPlat);
    bead.radiusShort = zeros(1,nPlat);
    bead.normalizedForce = zeros(1,nPlat);
    bead.errorLongRadius = zeros(1,nPlat);
    bead.errorLongForce = zeros(1,nPlat);
    bead.errorShortRadius = zeros(1,nPlat);
    bead.errorShortForce = zeros(1,nPlat);
    bead.cornerFreqLong = zeros(1,nPlat);
    bead.cornerFreqShort = zeros(1,nPlat);
    bead.signalLong = zeros(1,nPlat);
    bead.signalShort = zeros(1,nPlat);

    %%% Use the script "analyzeTrace" to determine the force for each trace
    for k=1:nPlat

        [extensionDNA, L, timeDomainForceLong, timeDomainForceShort, fitLong, ...
            forceLong, radiusLong, cornerFreqLong, signalLong, errorLong, fitShort,...
            forceShort, radiusShort, cornerFreqShort, signalShort, errorShort]=...
        analyzeTrace(...
            bead.long(plat(k).first:plat(k).last),... 
            bead.short(plat(k).first:plat(k).last),...
            bead.z(plat(k).first:plat(k).last),...
        sampleFreq, beadRadius, kT, viscosity, plotSpectrum);

        bead.extensionDNA(k) = extensionDNA;
        bead.L(k) = L;
        bead.timeDomainForceLong(k) = timeDomainForceLong;
        bead.timeDomainForceShort(k) = timeDomainForceShort;
        bead.fitLong(k) = fitLong;
        bead.forceLong(k) = forceLong;
        bead.radiusLong(k) = radiusLong;
        bead.fitShort(k) = fitShort;
        bead.forceShort(k) = forceShort;
        bead.radiusShort(k) = radiusShort;
        bead.normalizedForce(k) = forceShort/forceLong;
        bead.errorLongRadius(k) = errorLong(2,2) - radiusLong;
        bead.errorLongForce(k) = errorLong(1,2) - forceLong;
        bead.errorShortRadius(k) = errorShort(2,2) - radiusShort;
        bead.errorShortForce(k) = errorShort(1,2) - forceShort;
        bead.cornerFreqLong(k) = cornerFreqLong;
        bead.cornerFreqShort(k) = cornerFreqShort;
        bead.signalLong(k) = signalLong;
        bead.signalShort(k) = signalShort;

        display(['Finished working on bead, plateau number ' num2str(k)]);
    end

    if plotThings
        figure(5)
        semilogx([0.1, bead.forceLong(1)+1],[1, 1],'k--')
        hold on
        semilogx(bead.forceLong(ind),bead.normalizedForce(ind),'ro')
        axis([0.1, bead.forceLong(1)+1, 0.5, 1.5])
        title('Relative forces, long and short pendulum')
        xlabel('force measured in long pendulum direction (pN)')
        ylabel('normalized force')
        legend('1, when both directions give the same force','Short pendulum direction')
        hold off

        figure(6)
        p1 = semilogx(bead.forceLong(ind),bead.radiusLong(ind),'bo');
        hold on
        errorbar(bead.forceLong(ind),bead.radiusLong(ind),bead.errorLongRadius(ind),'bo');
        p2 = semilogx([0.1, bead.forceLong(1)+1],[sum(bead.radiusLong(ind).*bead.errorLongRadius(ind).^(-1))/sum(bead.errorLongRadius(ind).^(-1)), sum(bead.radiusLong(ind).*bead.errorLongRadius(ind).^(-1))/sum(bead.errorLongRadius(ind).^(-1))],'b--');
        p3 = semilogx(bead.forceShort(ind),bead.radiusShort(ind),'ro');
        errorbar(bead.forceShort(ind),bead.radiusShort(ind),bead.errorShortRadius(ind),'ro');
        p4 = semilogx([0.1, bead.forceLong(1)+1],[sum(bead.radiusShort(ind).*bead.errorShortRadius(ind).^(-1))/sum(bead.errorShortRadius(ind).^(-1)), sum(bead.radiusShort(ind).*bead.errorShortRadius(ind).^(-1))/sum(bead.errorShortRadius(ind).^(-1))],'r--');
        p5 = semilogx([0.1, bead.forceLong(1)+1],[beadRadius, beadRadius],'k--');
        axis([0.1, bead.forceLong(1)+1, mean(bead.radiusLong(ind)) - 300, mean(bead.radiusLong(ind)) + 300])
        title('Radius, fitted for different directions')
        xlabel('force (pN)')
        ylabel('radius (nm)')
        legend([p1 p2 p3 p4 p5],{'Fitted radii, long pendulum direction','Weighted average radius, long pendulum direction',...
            'Fitted radii, short pendulum direction', 'Weighted average radius, short pendulum direction',...
            'Radius given by manufacturer'})
        hold off          
            
        display('Fit force vs. magnet height to double exponential');
        [forcesExponentialFit] = fitExponentialToForces(bead, zmags(ind), configVariable);
    end
end