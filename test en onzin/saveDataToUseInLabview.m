for i = 1:nPlat;
    Testdata = zeros(plat(i).last-plat(i).first+1,4);
    Testdata(:,1) = 114*(1:length(Testdata(:,1)));
    Testdata(:,2) = bead(1).long(plat(i).first:plat(i).last)/1000;
    Testdata(:,3) = bead(1).short(plat(i).first:plat(i).last)/1000;
    Testdata(:,4) = bead(1).z(plat(i).first:plat(i).last)/1000;
    name = ['Testdata_plateau_' num2str(i) '_ext_' num2str(bead(1).ext(i)/1000)]
    dlmwrite(name,Testdata,'\t');

end