%--------------------------------------------------------------------------
% * Kaichao Liang, 2021.1.14
% * Poisson likelyhood with TV relulization optimization for 3D
% coded aperture reconstruction based on online projection and backprojection.
% The first dimension is spectral dimension with L2Norm regulization, the 
% other two dimensions are spatial dimension with L1Norm TV. 
% * Split Bregman method with EM integrated. 
%-----------------Objective  function-------------------
% * min mu: l(g,f) + lambda1/2*(f-mu-bmu)^2 + betaTVSpe/2*(divmuSpe)^2 + 
% betaTVSpaY*|dY| + lambdaY/2*(dY-divmuSpaY-bdivmuSpaY)^2  + 
% betaTVSpaZ*|dX| + lambdaX/2*(dZ-divmuSpaZ-bdivmuSpaX)^2
%-------------------------------------------------------
% * Reference: "Pascal et.al. 2012, Rudin-Osher-Fatemi Total Variation Denoising 
% using Split Bregman"
%--------------------------------------------------------------------------

function estimation = OSEMTV3DXRDT(rawData, sysParameter, TransSys, GeoSys, DiffSys, parameter)
%------------------------------Parameter-----------------------------------
% * rawData: panel data cell, each cell restores the rawData under an angle.
% rawData{ang} size [Energy, DetZ, DetY]
% * SysParameter: key parameter modeling the system.
% * TransSys: transmission spectrum cell, each cell restores the transmission
% data under an angle. TransSys{ang} size [Enery*DetY, PixY*PixX].
% * GeoSys: coded aperture geometry cell, each cell restores the geometry
% transpotation under an angle. GeoSys{ang} size [DetZ*DetY,PixY*PixX].
% * parameter: EMTV reconstruction hyper-parameters.
%--------------------------------------------------------------------------

AngleSet = sysParameter.AngleSet;
numAngle = numel(AngleSet);
SourcePos = sysParameter.SourcePos;

%物体参数
PixX = sysParameter.PixX;
PixY = sysParameter.PixY;
PixelSize = sysParameter.PixelSize;

%探测器参数
DetPos = sysParameter.DetPos;
DetY = sysParameter.DetY;
DetZ = sysParameter.DetZ;
offsetY = sysParameter.offsetY;
offsetZ = sysParameter.offsetZ;
DetSize = sysParameter.DetSize;
DetMask = sysParameter.DetMask;

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

% TV weighting spatial dimension Y
if(isfield(parameter,'betaTVSpaY'))
    betaTVSpaY = parameter.betaTVSpaY;
else
    betaTVSpaY = 1;
end

% TV weighting spatial dimension X
if(isfield(parameter,'betaTVSpaX'))
    betaTVSpaX = parameter.betaTVSpaX;
else
    betaTVSpaX = 1;
end

% Split Bregman lambda for mu
if(isfield(parameter,'lambda1'))
    lambda1 = parameter.lambda1;
else
    lambda1 = 1;
end

% Split Bregman lambda for divY
if(isfield(parameter,'lambdaY'))
    lambdaY = parameter.lambdaY;
else
    lambdaY = 1;
end

% Split Bregman lambda for divX
if(isfield(parameter,'lambdaX'))
    lambdaX = parameter.lambdaX;
else
    lambdaX = 1;
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
    DetY,DetZ,SpeSet,AngleSet(angle),DiffSys,GeoSys{angle},TransSys{angle});

    transScaleAng = CodedXRDTBackward(rawDataOnes,PixX,PixY,...
    DetY,DetZ,SpeSet,AngleSet(angle),DiffSys,GeoSys{angle},TransSys{angle});
    
    imageBK = imageBK + estimationAng;
    transScale = transScale + transScaleAng;
    transScaleAngSet{angle} = max(transScaleAng,1e-4); %pre-restore;
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
currentAngle = 0;
for totalstep = 1:totalSPStep
    
    currentAngle = currentAngle+1;
    if(currentAngle>=numAngle)
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
        DetY,DetZ,SpeSet,AngleSet(angle),DiffSys,GeoSys{angle},TransSys{angle});
        forward = max(forward,1e-4);
        scaleRaw = rawData{angle}./forward;
        scaleBK = CodedXRDTBackward(scaleRaw,PixX,PixY,...
        DetY,DetZ,SpeSet,AngleSet(angle),DiffSys,GeoSys{angle},TransSys{angle});
            
        Z = scaleBK.*f./transScaleAngSet{angle};

        % M step
        S = estimation-1/lambda1+bmu;
        f = S/2+sqrt((S/2).^2+1/lambda1*Z);
    end
   
    %-----------------------------------
    % 2. Subproblem mu, One sweep of Gauss-Seidel iteration
    % lambda1/2*(f-mu-bmu)^2 + betaTVSpe/2*(divmuSpe)^2
    %  + lambdaY/2*(dY-divmuSpaY-bdivmuSpaY)^2 + lambdaX/2*(dX-divmuSpaX-bdivmuSpaX)^2
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
        
        estimationSpaYExt=zeros(size(estimation,1),size(estimation,2)+2,size(estimation,3));
        estimationSpaYExt(:,2:end-1,:)=estimation; 
        estimationSpaYExt(:,1,:)=estimation(:,1,:);
        estimationSpaYExt(:,end,:)=estimation(:,end,:);
        
        estimationSpaXExt=zeros(size(estimation,1),size(estimation,2),size(estimation,3)+2);
        estimationSpaXExt(:,:,2:end-1)=estimation; 
        estimationSpaXExt(:,:,1)=estimation(:,:,1);
        estimationSpaXExt(:,:,end)=estimation(:,:,end);
        
        estimation = (lambda1/(2*lambdaX+2*lambdaY+2*betaTVSpe+lambda1))*(f-bmu)...
        +(lambdaY/(2*lambdaX+2*lambdaY+2*betaTVSpe+lambda1))*(estimationSpaYExt(:,1:end-2,:)+ estimationSpaYExt(:,3:end,:)-dYExt(:,2:end,:)+dYExt(:,1:end-1,:)+bdivmuSpaYExt(:,2:end,:)-bdivmuSpaYExt(:,1:end-1,:)) ...
        +(lambdaX/(2*lambdaX+2*lambdaY+2*betaTVSpe+lambda1))*(estimationSpaXExt(:,:,1:end-2)+ estimationSpaXExt(:,:,3:end)-dXExt(:,:,2:end)+dXExt(:,:,1:end-1)+bdivmuSpaXExt(:,:,2:end)-bdivmuSpaXExt(:,:,1:end-1))...
        + (betaTVSpe/(2*lambdaX+2*lambdaY+2*betaTVSpe+lambda1))*(estimationSpeExt(1:end-2,:,:)+estimationSpeExt(3:end,:,:));
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
    
    cost = sum(sum(sum(lambda1/2*(f-estimation-bmu).^2)))+sum(sum(sum(lambdaX/2*(dX-divmuSpaX-bdivmuSpaX).^2)))+sum(sum(sum(lambdaY/2*(dY-divmuSpaY-bdivmuSpaY).^2)))+sum(sum(sum(betaTVSpe/2*(divmuSpe).^2)));
    fprintf('Current step %d, mu subproblem cost: %f\n',totalstep,cost);
   
    %-----------------------------------
    % 3. Subproblem dY, shrinkage
    % betaTVSpaY*|dY| + lambdaY/2*(dY-divmuSpaY-bdivmuSpaY)^2
    %-----------------------------------
    dY = (divmuSpaY+bdivmuSpaY)./(max(1e-4,abs(divmuSpaY+bdivmuSpaY))).*max(abs(divmuSpaY+bdivmuSpaY)-betaTVSpaY/lambdaY,0);
    
    cost = sum(sum(sum(betaTVSpaY*abs(dY))))+sum(sum(sum(lambdaY/2*(dY-divmuSpaY-bdivmuSpaY).^2)));
    fprintf('Current step %d, dY subproblem cost: %f\n',totalstep,cost);
    
    %-----------------------------------
    % 4. Subproblem dX, shrinkage
    % betaTVSpaX*|dX| + lambdaX/2*(dX-divmuSpaX-bdivmuSpaX)^2
    %-----------------------------------
    dX = (divmuSpaX+bdivmuSpaX)./(max(1e-4,abs(divmuSpaX+bdivmuSpaX))).*max(abs(divmuSpaX+bdivmuSpaX)-betaTVSpaX/lambdaX,0);
    
    cost = sum(sum(sum(betaTVSpaX*abs(dX))))+sum(sum(sum(lambdaX/2*(dX-divmuSpaX-bdivmuSpaX).^2)));
    fprintf('Current step %d, dX subproblem cost: %f\n',totalstep,cost);
    
    bdivmuSpaX = bdivmuSpaX+divmuSpaX-dX;
    bdivmuSpaY = bdivmuSpaY+divmuSpaY-dY;
    bmu=bmu+estimation-f;
    if(mod(totalstep,10)==0)
        SaveName = sprintf("ReconStep%d.mat",totalstep);
        SaveName = strcat(savePath,'/',SaveName);
        ReconImage =estimation;
        save(SaveName,'ReconImage');
    end
end
end
    