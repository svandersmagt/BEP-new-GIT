% function [] = stuckBeadTest
configVariable.plotThings = true;

data = load('..\BEPmaarnietgit\stuckReferenceBead\RefBeadForceTestX.txt');
motorData = load('..\BEPmaarnietgit\stuckReferenceBead\RefBeadForceTestX_motors.txt');

stuckBead.x = data(:,3)*1000;
stuckBead.y = data(:,4)*1000;
stuckBead.z = data(:,5)*1000;

zmag = motorData(:,3);


[plat, zmags, nPlat] = plateauFinding(zmag,configVariable);
%%
figure(10)

x = stuckBead.x;
x(end+1) = x(1);
nBlock = 1;

[f, PSD, ~] = calcPSDBlock(x,500,nBlock);
f(1) = []; PSD(1) = [];
f(1) = []; PSD(1) = [];
ind = f < 10;
f=f(ind);PSD=PSD(ind);

exponential = @(a,b,freq)(exp(a-b.*log(freq)));
par  = lsqnonlin(@(par) (exponential(par(1),par(2),f)-PSD).*sqrt(1./PSD),[1,1]);
noise = exponential(par(1),par(2),f);
    
loglog(f,PSD)
hold on
loglog(f,noise,'r')

% for i = 2:nPlat;
%     x = stuckBead.x(plat(i).first:plat(i).last);
%     x(end+1) = x(1);
%     nBlock = 1;
%     
%     [f, PSD, ~] = calcPSDBlock(x,500,nBlock);
%     f(1) = []; PSD(1) = [];
%     
%     loglog(f,PSD);
%     hold on;
%     
%     clear x;
% end


hold off

