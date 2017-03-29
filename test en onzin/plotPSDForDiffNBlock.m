%%
clear all; clc; close all;

%%% Datafile and aquisition frequency
traces_file = '..\BEPmaarnietgit\DaldropData.txt';
fs = 2800; %%% Acquisition frequency in Hz

%%% Constants (pN and nm) 
kT = 4.1; %pN nm
eta = 10E-10; %Viscosity in pN s/nm^2
R = 515; %Bead radius in nm

%%% Read in data
data = load(traces_file);    

%%% Parse the bead data
for i=1;
bead(i).time = 1:length(data(:,1));
bead(i).x = data(:,1)*1000;
bead(i).y = data(:,2)*1000;
bead(i).z = data(:,3)*1000; %data is in mum, but this script uses nm

end

%%% Get an estimate of the force and the corner frequency, using the y
%%% direction, because x gives wrong values because it does not take
%%% rotation into account
Fest = kT*mean(bead.z)/std(bead.y-mean(bead.y))^2; %pN
fc = calcFcorner(Fest,mean(bead.z),R,eta); %Hz
fitgood = fc < fs/2;
bead.x(end+1) = bead.x(1);
[f1, PSD1, ~] = calcPSDBlock(bead.x,fs,1);
f1(1) = []; PSD1(1) = [];
PSD1 = PSD1.*1/2;
goodinds1 = f1 > fc(1)/20;
[f5, PSD5, ~] = calcPSDBlock(bead.x,fs,5);
f5(1) = []; PSD5(1) = [];
PSD5 = PSD5.*5/6;
goodinds5 = f5 > fc(1)/20;
[f40, PSD40, ~] = calcPSDBlock(bead.x,fs,40);
f40(1) = []; PSD40(1) = [];
PSD40 = PSD40.*40/41;
goodinds40 = f40 > fc(1)/20;

figure;
loglog(f1(goodinds1),PSD1(goodinds1),'r-');
hold on
loglog(f5(goodinds5),PSD5(goodinds5),'y-');
loglog(f40(goodinds40),PSD40(goodinds40),'b-');

title('Power spectra for different number of blocks');
xlabel('frequency (Hz)');
ylabel('Power Spectrum (nm^2/Hz)')
legend('1 block', '5 blocks','40 blocks');
hold off