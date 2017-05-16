%run eerst een bead, bijv rdrp1, met maxmagnetheight 16.

P=50;
kbT=4.1;
Lo=0.34*20600;
z=100:5000;
F=kbT./P.*((1./(4.*(1-z./Lo).^2)-1/4+z./Lo));
rmsd=sqrt(z.^2+2*kbT*z./F+2.*P.*Lo./(2+(1-z./Lo).^-3));
figure(11)
plot(z,F)
hold on
plot(rmsd,F,'r')
figure(12)
plot(F,z./rmsd)
hold on
plot(bead.forceLong,bead.L./bead.extensionDNA,'r')
axis([-0.01, 0.3, 0, 1.1])