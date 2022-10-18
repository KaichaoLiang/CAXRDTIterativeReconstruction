load AdiposeXRDGT;
load CancerXRDGT;
load GlanXRDGT;
Phantom = zeros(200,540,540);
Mask = zeros(540,540); PixelSize=0.1;
AdiposeXRD=AdiposeXRD(2:2:400);
GlanXRD=GlanXRD(2:2:400);
CancerXRD=CancerXRD(2:2:400);

PositionX=([1:540]'-270.5)*PixelSize;
PositionY=([1:540]-270.5)*PixelSize;

%AdiposeBK
CenterAdiX = 0;
CenterAdiY = 0;
SizeR = 25;
Mask(find((PositionX-CenterAdiX).^2+(PositionY-CenterAdiY).^2<=SizeR.^2))=2;


%Carcinoma
CenterCancerX = -15;
CenterCancerY = 0;
SizeR = 4;
Mask(find((PositionX-CenterCancerX).^2+(PositionY-CenterCancerY).^2<=SizeR.^2))=3;
imshow(Mask,[])
%Carcinoma
CenterCancerX = -15*cos(45/180*pi);
CenterCancerY = 15*sin(45/180*pi);
SizeR = 2;
Mask(find((PositionX-CenterCancerX).^2+(PositionY-CenterCancerY).^2<=SizeR.^2))=3;
imshow(Mask,[])
%Carcinoma
CenterCancerX = -15*cos(90/180*pi);
CenterCancerY = 15*sin(90/180*pi);
SizeR = 1;
Mask(find((PositionX-CenterCancerX).^2+(PositionY-CenterCancerY).^2<=SizeR.^2))=3;
imshow(Mask,[])

%Glandular
CenterCancerX = -15*cos(pi);
CenterCancerY = 15*sin(pi);
SizeR = 4;
Mask(find((PositionX-CenterCancerX).^2+(PositionY-CenterCancerY).^2<=SizeR.^2))=4;
imshow(Mask,[])
%Carcinoma
CenterCancerX = -15*cos(225/180*pi);
CenterCancerY = 15*sin(225/180*pi);
SizeR = 2;
Mask(find((PositionX-CenterCancerX).^2+(PositionY-CenterCancerY).^2<=SizeR.^2))=4;
imshow(Mask,[])
%Carcinoma
CenterCancerX = -15*cos(270/180*pi);
CenterCancerY = 15*sin(270/180*pi);
SizeR = 1;
Mask(find((PositionX-CenterCancerX).^2+(PositionY-CenterCancerY).^2<=SizeR.^2))=4;
imshow(Mask,[])

Phantom=reshape(Phantom,200,540*540);
Phantom(:,find(Mask(:)==2))=repmat(AdiposeXRD,1,numel(find(Mask(:)==2)));
Phantom(:,find(Mask(:)==3))=repmat(CancerXRD,1,numel(find(Mask(:)==3)));
Phantom(:,find(Mask(:)==4))=repmat(GlanXRD,1,numel(find(Mask(:)==4)));
Phantom=reshape(Phantom,200,540,540);
imshow(reshape(Phantom(80,:,:),540,540),[]);
Phantom54=zeros(200,54,54);

for x=1:10
    for y=1:10
        Phantom54=Phantom54+Phantom(:,x:10:end,y:10:end)/100;
    end
end
imshow(reshape(Phantom54(56,:,:),54,54),[]);
Phantom=Phantom54;
Phantom(:,27:28,12:13)=repmat(GlanXRD,1,2,2);
Phantom(:,27:28,42:43)=repmat(CancerXRD,1,2,2);
save Phantom Phantom

