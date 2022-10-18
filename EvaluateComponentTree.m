addpath('GTData');
load('AdiposeXRDGT.mat');
load('GlanXRDGT.mat');
load('CancerXRDGT.mat');
AdiposeXRD=AdiposeXRD(101:2:360);
GlanXRD=GlanXRD(101:2:360);
CancerXRD=CancerXRD(101:2:360);
Result=zeros(54,54);
Similarity=zeros(54,54);
for x=1:54
    for y=1:54
        if(max(ReconImage(51:180,y,x))>2)
            SimiGlanMax=cov(ReconImage(51:180,y,x),GlanXRD);
            SimiGlan=SimiGlanMax(1,2)/sqrt(SimiGlanMax(1,1)*SimiGlanMax(2,2));
            
            SimiCancerMax=cov(ReconImage(51:180,y,x),CancerXRD);
            SimiCancer=SimiCancerMax(1,2)/sqrt(SimiCancerMax(1,1)*SimiCancerMax(2,2));
            
            SimiAdiMax=cov(ReconImage(51:180,y,x),AdiposeXRD);
            SimiAdi=SimiAdiMax(1,2)/sqrt(SimiAdiMax(1,1)*SimiAdiMax(2,2));
            if(SimiAdi<0.93 || (SimiAdi<SimiGlan) ||(SimiAdi<SimiCancer))
                SimiMatrix=[SimiCancer,SimiGlan];
                [m,p]=max(SimiMatrix);
            
                Result(y,x)=p+1;
                Similarity(y,x)=m;
            else
                Result(y,x)=1;
                Similarity=SimiAdi;
            end
            
        end
    end
end
imshow(Result,[])
           
        