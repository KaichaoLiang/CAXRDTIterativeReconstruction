% Select data and settings
clear;clc;
DataRoot = 'ResultsOSEMTV_3ViewTV6e-3';
StartImage = 4000;
IntervalIdeal=100;
ThreadRRMSE=0.001;
NumAngle = 3;
Rounds = round(IntervalIdeal/NumAngle);
Interval=NumAngle*Rounds;
ThreadRRMSE=ThreadRRMSE/IntervalIdeal*Interval;

File = strcat(DataRoot,'/',sprintf('ReconStep%d.mat',StartImage+Interval));
RRMSESet = [];
CurrentImg=StartImage;
while(exist(File))
    File1 = strcat(DataRoot,'/',sprintf('ReconStep%d.mat',CurrentImg));
    load(File1);
    Recon1=ReconImage(46:180,:,:);
    File2 = strcat(DataRoot,'/',sprintf('ReconStep%d.mat',CurrentImg+Interval));
    load(File2);
    Recon2=ReconImage(46:180,:,:);
    
    RRMSE=norm(Recon2(:)-Recon1(:))./norm(Recon1(:));
    RRMSESet(end+1)=RRMSE;
    if(RRMSE<ThreadRRMSE)
        %ReconImage=Recon2;
        S=strcat(DataRoot,'/ReconImage.mat');
        save(S,'ReconImage');
        fprintf('The end image is %d\n',CurrentImg+Interval);
        break;
    end
    fprintf('current image%d, RRMSE %f\n',CurrentImg, RRMSE);
    CurrentImg=CurrentImg+1;
    File = strcat(DataRoot,'/',sprintf('ReconStep%d.mat',CurrentImg+Interval));
end