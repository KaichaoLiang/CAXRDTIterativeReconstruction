%--------------------------------------------------------------------------
% * Kaichao Liang, 2022.06.13
% * Calculate the Raylscattering photon production rate based on geometry.
% * Parameter:
%   SourcePos: position from source to rotation center.
%   FanAngle: the fan-beam angle of X-ray source
%   PixelSize: the size of object pixel.
%   DetSize: the size of detector pixel.
%--------------------------------------------------------------------------
function RaylRate = CalculateRalyProductionRate(SourcePos,FanAngle,DetPos,PixelSize,DetSize,dEnergy)
    CosAngleIntergral=0;
    discreteAngleNum=100;
    FanAngle=FanAngle/180*pi;
    deltaAngle=FanAngle/discreteAngleNum;
    
    
    for angleIndex = 0:discreteAngleNum-1
        curAngle=-FanAngle/2+0.5*deltaAngle+angleIndex*deltaAngle;
        CosAngleIntergral=CosAngleIntergral+cos(curAngle)*deltaAngle;
    end
    CentralPixAngle=PixelSize/SourcePos;
    ScatteringSolidAngle=(DetSize/DetPos)^2;
    
    re=2.82e-13;%cm
    NA=6.02e23;
    
    %rate Inp*re^2*NA*Solid/10
    RaylRate = CentralPixAngle/CosAngleIntergral*(re)^2*NA*ScatteringSolidAngle*(PixelSize)*dEnergy/10;
end

