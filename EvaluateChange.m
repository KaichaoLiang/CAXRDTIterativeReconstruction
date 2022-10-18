Dir = 'ResultsOSEMTV_l11e-2_lXY1e-2_betaSpe3e-2_betaSpa1e-2FullGSFullRand/';
ChangeSet = zeros(230,1);
for i =1:230
    currentName = strcat(Dir,sprintf('ReconStep%d.mat',i*10));
    load(currentName);
    ReconC = ReconImage(41:200,:,:);
    
    nextName = strcat(Dir,sprintf('ReconStep%d.mat',i*10+210));
    load(nextName);
    ReconN = ReconImage(41:200,:,:);
    ChangeSet(i)=norm(ReconC(:)-ReconN(:))/norm(ReconC(:));
end
plot(ChangeSet)
