ind = zmags <10;

thisBead.magnetHeight = zmags(ind);
thisBead.forceLong = bead.forceLong;
thisBead.forceShort = bead.forceShort;
thisBead.radiusLong = bead.radiusLong;
thisBead.radiusShort = bead.radiusShort;
thisBead.errorLongRadius = bead.errorLongRadius;
thisBead.errorShortRadius = bead.errorShortRadius;
thisBead.errorLongForce = bead.errorLongForce;
thisBead.errorShortForce = bead.errorShortForce;
thisBead.exponentialFitParameters = forcesExponentialFit;

RdRp19 = thisBead;
save('RdRp19.mat','RdRp19');

%%
load('RdRp1.mat');
load('RdRp2.mat');
load('RdRp3.mat');
load('RdRp4.mat');
load('RdRp8.mat');
load('RdRp14.mat');
load('RdRp17.mat');
load('RdRp19.mat');

Results.RdRp1 = RdRp1;
Results.RdRp2 = RdRp2;
Results.RdRp3 = RdRp3;
Results.RdRp4 = RdRp4;
Results.RdRp8 = RdRp8;
Results.RdRp14 = RdRp14;
Results.RdRp17 = RdRp17;
Results.RdRp19 = RdRp19;

save('Results.mat','Results');


