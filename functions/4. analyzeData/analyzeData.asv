function [bead, forcesExponentialFit] = analyzeData(bead, nPlat, plat, zmags, configVariable)
%%% STEP 4: ANALYZE THE BEAD DATA FOR EACH PLATEAU
%%% --------------------------------------------------------------------
%%% In this step the bead traces are analyzed. See analyzeTrace.m and
%%% fitExponentialToForces.m

%%% Input: (bead, nPlat, plat, zmags, configVariable)
%%% - struct containing the time trace and position traces of the beads
%%% - number of plateaus
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
    nBlock = configVariable.nBlock; %number of blocks to use for blocked fourier transform
    sampleFreq = configVariable.sampleFreq;
    plotThings = configVariable.plotThings;

    for i=1;

        bead(i).extensionDNA = zeros(1,nPlat);
        bead(i).L = zeros(1,nPlat);
        bead(i).timeDomainForceLong = zeros(1,nPlat);
        bead(i).timeDomainForceShort = zeros(1,nPlat);
        bead(i).fitLong = zeros(1,nPlat);
        bead(i).forceLong = zeros(1,nPlat);
        bead(i).radiusLong = zeros(1,nPlat);
        bead(i).fitShort = zeros(1,nPlat);
        bead(i).forceShort = zeros(1,nPlat);
        bead(i).radiusShort = zeros(1,nPlat);
        bead(i).normalizedForce = zeros(1,nPlat);

        %%% Use the script "analyzeTrace" to determine the force for each trace
        for k=1:nPlat

            [extensionDNA, L, timeDomainForceLong, timeDomainForceShort, fitLong, ...
                forceLong, radiusLong, fitShort, forceShort, radiusShort]=...
            analyzeTrace(...
                bead(i).long(plat(k).first:plat(k).last),... 
                bead(i).short(plat(k).first:plat(k).last),...
                bead(i).z(plat(k).first:plat(k).last),...
            sampleFreq, beadRadius, kT, viscosity, nBlock, plotThings);

            bead(i).extensionDNA(k) = extensionDNA;
            bead(i).L(k) = L;
            bead(i).timeDomainForceLong(k) = timeDomainForceLong;
            bead(i).timeDomainForceShort(k) = timeDomainForceShort;
            bead(i).fitLong(k) = fitLong;
            bead(i).forceLong(k) = forceLong;
            bead(i).radiusLong(k) = radiusLong;
            bead(i).fitShort(k) = fitShort;
            bead(i).forceShort(k) = forceShort;
            bead(i).radiusShort(k) = radiusShort;
            bead(i).normalizedForce(k) = forceShort/forceLong;

            display(['Finished working on bead # ' num2str(i) ', plateau number ' num2str(k)]);
        end
        
        %to exclude forces too low for accurate fitting
        ind = bead.forceLong > 0.1;

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
            semilogx(bead.forceLong(ind),bead.radiusLong(ind),'ro')
            hold on
            semilogx([0.1, bead.forceLong(1)+1],[mean(bead.radiusLong(ind)), mean(bead.radiusLong(ind))],'r--')
            semilogx(bead.forceShort(ind),bead.radiusShort(ind),'bo')
            semilogx([0.1, bead.forceLong(1)+1],[mean(bead.radiusShort(ind)), mean(bead.radiusShort(ind))],'b--')
            semilogx([0.1, bead.forceLong(1)+1],[beadRadius, beadRadius],'k--')
            axis([0.1, bead.forceLong(1)+1, beadRadius - 300, beadRadius + 300])
            title('Radius, fitted for different directions')
            xlabel('force (pN)')
            ylabel('radius (nm)')
            legend('Fitted radius, long pendulum direction','Average radius, long pendulum direction',...
                'Fitted radii, short pendulum direction', 'Average radius, short pendulum direction',...
                'Radius given by manufacturer')
            hold off          
        end
            
        display('Fit force vs. magnet height to double exponential');
        [forcesExponentialFit] = fitExponentialToForces(bead, zmags, configVariable, ind);
    end
end