%--------------------------------------------------------------------------
% * Kaichao Liang, 2021.1.18
% * Simplified FBP reconstruction function for equal distance fan-beam CT
% * Need: Projection [det,angle] and geoParameter.
%--------------------------------------------------------------------------

function Reconstruction = FanBeamFBP(Projection, geoParameter)
    %---------------------------------------
    % Parameter define 
    LSource2Center = geoParameter.LSource2Center;
    LCenter2Det = geoParameter.LCenter2Det;
    Ndet = geoParameter.Ndet; 
    BinSize = geoParameter.BinSize;
    NPhi = geoParameter.NPhi;
    Nx = geoParameter.Nx;
    Ny = geoParameter.Ny;
    PixelSize = geoParameter.PixelSize;
    
    DetStart=(-Ndet/2+0.5)*BinSize;
    DetEnd=(Ndet/2-0.5)*BinSize;
    DetPos=DetStart:BinSize:DetEnd;
    DPhi = 2*pi/NPhi;
    PhiSet=double(1:NPhi)/NPhi*2*pi;
    
    xx=repmat([1:Nx],Ny,1)-Nx/2-0.5;
    yy=repmat([1:Ny]',1,Nx)-Ny/2-0.5;
    FOVIndex=find(xx.^2+yy.^2<(Nx/2-0.5).^2);
    FOV=zeros(Ny,Nx);
    FOV(FOVIndex)=1;

    PosObX=([1:Nx]-Nx/2-0.5).*PixelSize;
    PosObY=([1:Ny]-Ny/2-0.5).*PixelSize;


    %Step one: weighting
    R1=LSource2Center+LCenter2Det;
    Weights=R1./sqrt(DetPos.^2+R1.^2);
    Weights=repmat(reshape(Weights,Ndet,1),1,NPhi);%��ͬ�Ƕ���Weights��ͬ
    Weights=single(Weights);
    PrjW=Projection.*Weights;

    %Step two: filtering along x direction
    %RL�˲�
    h = zeros(2*Ndet-1,1);
    h(Ndet) = 1/8/BinSize^2;
    h = zeros(2*Ndet-1,1);
    h(Ndet) = 1/8/BinSize^2;
    if mod(Ndet,2)==0
        h(1:2:2*Ndet-1) = -1/2 ./( pi^2 * (-Ndet+1:2:Ndet-1).^2 * BinSize^2 );
    else
        h(2:2:2*Ndet-2) = -1/2 ./( pi^2 * (-Ndet+2:2:Ndet-2).^2 * BinSize^2 );
    end
    h=h.*((LSource2Center+LCenter2Det)/LSource2Center)^2;%��������̽�����ߴ�BinSize��Ҫ���ŵ���ԭ��ֱ����
    h=single(h);
    PrjF=conv2(PrjW,h,'same');

    %Step three: Weighted backprojection
    Reconstruction=zeros(Ny,Nx,'single');
    %PrjF=fliplr(PrjF);
    parfor x=1:Nx
        Position_x=PosObX(x);
        for y=1:Ny
            Position_y=PosObY(y);
            if (FOV(y,x)>0)
                for angle=1:NPhi
                    beta=PhiSet(angle);
                    U=1+(sin(beta)*Position_y+cos(beta)*Position_x) /LSource2Center;
                    CoWeight=DPhi*1./(U.^2)*BinSize*LSource2Center/(LSource2Center+LCenter2Det);%����dtheta��theta���֣�����DetSize_x��̽�������֣�������Ҫ��̽�����ߴ����ŵ���ԭ�����
                    t_a=(LSource2Center+LCenter2Det)*(Position_y*cos(beta)-Position_x *sin(beta))/(LSource2Center+sin(beta)*Position_y+cos(beta)*Position_x);
                
                    index_a=t_a/BinSize+Ndet/2+0.5; %̽����ƽ��˫���Բ�ֵ
                    index_a_min=min(max(floor(index_a),1),Ndet);
                    index_a_max=min(Ndet,index_a_min+1);                              
               
                    Reconstruction(y,x)=Reconstruction(y,x)+single(PrjF(index_a_min,angle)*(1-index_a+index_a_min)*CoWeight);
                    Reconstruction(y,x)=Reconstruction(y,x)+single(PrjF(index_a_max,angle)*(index_a-index_a_min)*CoWeight);
                end
            end
        end
    end
end




