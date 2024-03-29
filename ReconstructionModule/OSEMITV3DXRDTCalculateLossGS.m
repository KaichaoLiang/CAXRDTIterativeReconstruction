%--------------------------------------------------------------------------
% * Kaichao Liang, 2022.10.18
% * Poisson likelyhood with TV relulization optimization for 3D
% coded aperture reconstruction based on online projection and backprojection.
% The first dimension is spectral dimension with L2Norm regulization, the 
% other two dimensions are spatial dimension with L1Norm TV. 
% * Split Bregman method with EM integrated. 
%-----------------Objective  function-------------------
% * min mu: l(g,f) + lambda1/2*(f-mu-bmu)^2 + betaTVSpe/2*(divmuSpe)^2 + 
% betaTVSpa*|dY^2+dX^2| + lambdaSpa/2*(dY-divmuSpaY-bdivmuSpaY)^2 + lambdaSpa/2*(dZ-divmuSpaZ-bdivmuSpaX)^2
%-------------------------------------------------------
% * Reference: "Pascal et.al. 2012, Rudin-Osher-Fatemi Total Variation Denoising 
% using Split Bregman"
%--------------------------------------------------------------------------

function estimation = OSEMITV3DXRDTCalculateLossGS(rawData, sysParameter, TransSys, GeoSys, DiffSys, ResponseMatrix,parameter)
%------------------------------Parameter-----------------------------------
% * rawData: panel data cell, each cell restores the rawData under an angle.
% rawData{ang} size [Energy, DetZ, DetY]
% * SysParameter: key parameter modeling the system.
% * TransSys: transmission spectrum cell, each cell restores the transmission
% data under an angle. TransSys{ang} size [Enery*DetY, PixY*PixX].
% * GeoSys: coded aperture geometry cell, each cell restores the geometry
% transpotation under an angle. GeoSys{ang} size [DetZ*DetY,PixY*PixX].
% * ResponseMatrix: Detector energy response [Energy, Energy]
% * parameter: EMTV reconstruction hyper-parameters.
%--------------------------------------------------------------------------
eps = 1e-4;
AngleSet = sysParameter.AngleSet;
numAngle = numel(AngleSet);

%物体参数
PixX = sysParameter.PixX;
PixY = sysParameter.PixY;

%探测器参数
DetY = sysParameter.DetY;
DetZ = sysParameter.DetZ;

%谱参数
SpeSet = sysParameter.SpeSet;
MTSet= sysParameter.MTSet;

reconSize = [numel(MTSet),PixY,PixX];

%-------------------------
% parameter define
%-------------------------
% TV weighting spectrum dimension
if(isfield(parameter,'betaTVSpe'))
    betaTVSpe = parameter.betaTVSpe;
else
    betaTVSpe = 1;
end

% TV weighting 
if(isfield(parameter,'betaTVSpa'))
    betaTVSpa = parameter.betaTVSpa;
else
    betaTVSpa = 1;
end

% Split Bregman lambda for mu
if(isfield(parameter,'lambda1'))
    lambda1 = parameter.lambda1;
else
    lambda1 = 1;
end

% Split Bregman lambda for divY
if(isfield(parameter,'lambdaSpa'))
    lambdaSpa = parameter.lambdaSpa;
else
    lambdaSpa = 1;
end

% Split Bregman steps
if(isfield(parameter,'totalSPStep'))
    totalSPStep = parameter.totalSPStep;
else
    totalSPStep = 1000;
end

% EM subproblem steps
if(isfield(parameter,'subEMStep'))
    subEMStep = parameter.subEMStep;
else
    subEMStep = 1;
end

% Gauss-Seidel subproblem steps
if(isfield(parameter,'subGSStep'))
    subGSStep = parameter.subGSStep;
else
    subGSStep = 1;
end

if(isfield(parameter,'savePath'))
    savePath = parameter.savePath;
else
    savePath = 'ResultsOSEMTV';
end

if(isfield(parameter,'ResumeStep'))
    load(strcat(savePath,'/',sprintf('ReconStep%d.mat',parameter.ResumeStep)));
    parameter.initImage=ReconImage;
    baseStep=parameter.ResumeStep;
else
    baseStep=0;
end


mkdir(savePath);
%-------------------------
% variable initialization
%-------------------------
% estimation = Hsys'*rawData./(max((Hsys'*ones(size(rawData))).^2,1e-8));
transScale = zeros(reconSize);
imageBK = zeros(reconSize);
rawDataOnes = ones(numel(SpeSet),DetZ,DetY);
transScaleAngSet = cell(numAngle,1);

for angle = 1:numAngle
    estimationAng = CodedXRDTBackward(rawData{angle},PixX,PixY,...
    DetY,DetZ,SpeSet,AngleSet(angle),DiffSys,GeoSys{angle},TransSys{angle},ResponseMatrix);

    transScaleAng = CodedXRDTBackward(rawDataOnes,PixX,PixY,...
    DetY,DetZ,SpeSet,AngleSet(angle),DiffSys,GeoSys{angle},TransSys{angle},ResponseMatrix);
    
    imageBK = imageBK + estimationAng;
    transScale = transScale + transScaleAng;
    transScaleAngSet{angle} = max(transScaleAng,1); %pre-restore;
end
%estimation = imageBK./max(transScale,1e-8);
estimation = ones(reconSize);
if(isfield(parameter,'initImage'))
    estimation = parameter.initImage;
end
    
f = estimation;
bmu = zeros(size(estimation));

% spatial differential Y
dY = zeros(size(f));
dY(:, 1:end-1,:) = estimation(:,2:end,:)-estimation(:,1:end-1,:);
dY(:,end,:) = 0;
bdivmuSpaY = zeros(size(dY));

% spatial differential X
dX = zeros(size(f));
dX(:,:,1:end-1) = estimation(:,:,2:end)-estimation(:,:,1:end-1);
dX(:,:,end) = 0;
bdivmuSpaX = zeros(size(dX));


%-------------------------
% iterative reconstruction
%-------------------------
CostSet = [];
currentAngle = 0;
for totalstep = baseStep+1:totalSPStep
    
    currentAngle = currentAngle+1;
    if(currentAngle>numAngle)
        currentAngle =1;
    end
    fprintf('Current step: %d, current angle%d--------\n',totalstep, currentAngle);
    
    %-----------------------------------
    % 1. Subproblem f, EM algorithm,updata each angle
    % l(g,f) + lambda1/2*(f-mu-bmu)^2
    %-----------------------------------
    for emstep = 1:subEMStep
        angle = currentAngle;
        % E step
        % hidden variable Z sumed by rawData dimension
        %Z = (Hsys'*(rawData./max(Hsys*f(:),1e-4))).*f(:)./(max(Hsys'*ones(size(rawData)),1e-4));
        forward = CodedXRDTForward(f,PixX,PixY, ...
        DetY,DetZ,SpeSet,AngleSet(angle),DiffSys,GeoSys{angle},TransSys{angle},ResponseMatrix);
        forward = max(forward,1);
        scaleRaw = rawData{angle}./forward;
        scaleBK = CodedXRDTBackward(scaleRaw,PixX,PixY,...
        DetY,DetZ,SpeSet,AngleSet(angle),DiffSys,GeoSys{angle},TransSys{angle},ResponseMatrix);
            
        Z = scaleBK.*f./transScaleAngSet{angle};

        % M step
        S = estimation-1/lambda1+bmu;
        f = S/2+sqrt((S/2).^2+1/lambda1*Z);
        
        cost1 = sum(forward(:)-rawData{angle}(:).*log(forward(:)))+sum(lambda1/2*(f(:)-estimation(:)-bmu(:)).^2);
        fprintf('Current step %d, mu subproblem cost: %f\n',totalstep,cost1);
        
    end
   
    %-----------------------------------
    % 2. Subproblem mu, One sweep of Gauss-Seidel iteration
    % lambda1/2*(f-mu-bmu)^2 + betaTVSpe/2*(divmuSpe)^2
    %  + lambdaSpa/2*(dY-divmuSpaY-bdivmuSpaY)^2 + lambdaSpa/2*(dX-divmuSpaX-bdivmuSpaX)^2
    %-----------------------------------
    % Padding
    dYExt=zeros(size(dY,1),size(dY,2)+1,size(dY,3));
    dYExt(:,2:end,:)=dY;
    dYExt(:,1,:)=dY(:,1,:);
    
    bdivmuSpaYExt=zeros(size(bdivmuSpaY,1),size(bdivmuSpaY,2)+1,size(bdivmuSpaY,3));
    bdivmuSpaYExt(:,2:end,:)=bdivmuSpaY;
    bdivmuSpaYExt(:,1,:)=bdivmuSpaY(:,1,:);
    
    dXExt=zeros(size(dX,1),size(dX,2),size(dX,3)+1);
    dXExt(:,:,2:end)=dX;
    dXExt(:,:,1)=dX(:,:,1);
    
    bdivmuSpaXExt=zeros(size(bdivmuSpaX,1),size(bdivmuSpaX,2),size(bdivmuSpaX,3)+1);
    bdivmuSpaXExt(:,:,2:end)=bdivmuSpaX;
    bdivmuSpaXExt(:,:,1)=bdivmuSpaX(:,:,1);
    
    % Gauss-Seidel iterationestimationExt=zeros(size(estimation,1)+2,1);
    for gsstep = 1:subGSStep  
        estimationSpeExt=zeros(size(estimation,1)+2,size(estimation,2),size(estimation,3));
        estimationSpeExt(2:end-1,:,:)=estimation;
        estimationSpeExt(1,:,:)=estimation(1,:,:);
        estimationSpeExt(end,:,:)=estimation(end,:,:);
       
        
        for x=1:PixX
            for y=1:PixY
                if(x==1)
                    beforeIndex = 1;
                else
                    beforeIndex = x-1;
                end
                
                if(x==PixX)
                    nextIndex = PixX;
                else
                    nextIndex = x+1;
                end
                
                if(y==1)
                    leftIndex = 1;
                else
                    leftIndex = y-1;
                end
                
                if(y==PixY)
                    rightIndex = PixY;
                else
                    rightIndex = y+1;
                end
                
                estimation(:,y,x) = (lambda1/(2*lambdaSpa+2*lambdaSpa+2*betaTVSpe+lambda1))*(f(:,y,x)-bmu(:,y,x))...
                +(lambdaSpa/(2*lambdaSpa+2*lambdaSpa+2*betaTVSpe+lambda1))*(estimation(:,leftIndex,x)+ estimation(:,rightIndex,x)-dYExt(:,y+1,x)+dYExt(:,y,x)+bdivmuSpaYExt(:,y+1,x)-bdivmuSpaYExt(:,y,x)) ...
                +(lambdaSpa/(2*lambdaSpa+2*lambdaSpa+2*betaTVSpe+lambda1))*(estimation(:,y,beforeIndex)+ estimation(:,y,nextIndex)-dXExt(:,y,x+1)+dXExt(:,y,x)+bdivmuSpaXExt(:,y,x+1)-bdivmuSpaXExt(:,y,x))...
                + (betaTVSpe/(2*lambdaSpa+2*lambdaSpa+2*betaTVSpe+lambda1))*(estimationSpeExt(1:end-2,y,x)+estimationSpeExt(3:end,y,x));
            end
        end
    end
    
    divmuSpe = zeros(size(estimation));
    divmuSpe(1:end-1,:,:) = diff(estimation);
    divmuSpe(end,:,:) =0;
    
    divmuSpaY = zeros(size(estimation));
    divmuSpaY(:,1:end-1,:) = estimation(:,2:end,:)-estimation(:,1:end-1,:);
    divmuSpaY(:,end,:) =0;
    
    divmuSpaX = zeros(size(estimation));
    divmuSpaX(:,:,1:end-1) = estimation(:,:,2:end)-estimation(:,:,1:end-1);
    divmuSpaX(:,:,end) =0;
    
    cost2 = sum(sum(sum(lambda1/2*(f-estimation-bmu).^2)))+sum(sum(sum(lambdaSpa/2*(dX-divmuSpaX-bdivmuSpaX).^2)))+sum(sum(sum(lambdaSpa/2*(dY-divmuSpaY-bdivmuSpaY).^2)))+sum(sum(sum(betaTVSpe/2*(divmuSpe).^2)));
    fprintf('Current step %d, mu subproblem cost: %f\n',totalstep,cost2);
   
    %-----------------------------------
    % 3. Subproblem dY, dX, shrinkage
    % betaTVSpa*|dY^2+dX^2| + lambdaSpa/2*(dY-divmuSpaY-bdivmuSpaY)^2 + lambdaSpa/2*(dX-divmuSpaX-bdivmuSpaX)^2
    %-----------------------------------
    S = sqrt((divmuSpaX+bdivmuSpaX).^2+(divmuSpaY+bdivmuSpaY).^2);
    dX = max(0,(S-betaTVSpa/lambdaSpa)).*(divmuSpaX+bdivmuSpaX)./(S+eps);
    dY = max(0,(S-betaTVSpa/lambdaSpa)).*(divmuSpaY+bdivmuSpaY)./(S+eps);
    
    cost3 = sum(sum(sum(betaTVSpa*sqrt(dX.^2+dY.^2))))+sum(sum(sum(lambdaSpa/2*(dX-divmuSpaX-bdivmuSpaX).^2)))++sum(sum(sum(lambdaSpa/2*(dY-divmuSpaY-bdivmuSpaY).^2)));
    fprintf('Current step %d, dX subproblem cost: %f\n',totalstep,cost3);
    
    bdivmuSpaX = bdivmuSpaX+divmuSpaX-dX;
    bdivmuSpaY = bdivmuSpaY+divmuSpaY-dY;
    bmu=bmu+estimation-f;
    if((totalstep>1200) ||(mod(totalstep,10)==0))
        SaveName = sprintf("ReconStep%d.mat",totalstep);
        SaveName = strcat(savePath,'/',SaveName);
        ReconImage =estimation;
        save(SaveName,'ReconImage');
    end
    CostSet(end+1)=cost1+cost2+cost3;
    SaveName = sprintf("CostSet.mat");
    SaveName = strcat(savePath,'/',SaveName);
    save(SaveName,'CostSet');
end
end
    