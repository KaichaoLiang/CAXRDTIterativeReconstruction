%--------------------------------------------------------------------------
% * Kaichao Liang, 2021.1.14
% * EMTV reconstruction of multi-angle coded XRDT panel data.
% * Including attenuation effect, first order coherent scattering and
% Compton scattering in Data folder.
% * Need 1. panel data, 2. GeoMatrix, 3. DecayMatrix and 
% 4. source spectrum, 5. Phantom attenuation data
%--------------------------------------------------------------------------
clear;clc;
dbstop if error;
addpath('GTData'); % Ground-truth cross section data.
addpath('DecayMatrix15View');
addpath('CodeGeoMatrix15ViewSimple');
addpath('TransMatrix');
addpath('RaylMatrix');
addpath('ReconstructionModule');
addpath('PanelData15View');

%p = parpool(18);
%=============================
% Parameter define
%=============================
numAngle = 15;
dAngle = pi/numAngle*2;
AngleSet = [1:numAngle]*dAngle;

SourcePos = 400;

%物体参数
PixX = 54;
PixY = 54;
PixelSize =1;

%探测器参数
DetPos = 297.5;
DetY = 64;
DetZ = 64;
offsetY = 0;
offsetZ = 0;
DetSize = 1.6;

SpeSet = 21:85;
MTSet= 0.01:0.02:4;

%=============================
% Energy response matrix
%=============================
% Energy response matrix
load ResponseMatrix

%=============================
% Attenuation estimation, use ground-truth firstly
% the attenuation can be from FBP reconstruction of 
% transmission signal.
%=============================
%Attenuation 
load PhantomDecay;

%=============================
% Generate transmission matrix
%=============================
load Spectrum;
Spectrum=Spectrum/sum(Spectrum)*0.9;
RaylRate = CalculateRalyProductionRate(SourcePos,9,DetPos,PixelSize,DetSize,1);
PhotonPerView=1e8/numAngle*21*600;
Spectrum=Spectrum*RaylRate*PhotonPerView;
mkdir('TransMatrix');
fprintf('calculate transmission model....\n======================================\n');
TransSys = cell(numAngle,1);
for angle = 1:numAngle %Calculate transmission spectrum for each angle, and restore the data
    DecayName = sprintf('DecaySysAngle%d.mat',angle);
    load(DecayName);
    eval(sprintf('DecaySys = DecaySysAngle%d;',angle));
    eval(sprintf('clear DecaySysAngle%d;', angle));
    TransMatrix = zeros(numel(SpeSet),DetY,PixY,PixX);
    
    fprintf('Calculate transmission spectrum, angle:%d\n',angle);
    
    for energy =1:numel(SpeSet)
        PhantomDecayE = PhantomDecay(energy,:,:);
        TransE = exp(-reshape(DecaySys*PhantomDecayE(:),1,DetY,PixY,PixX)).*Spectrum(energy);
        TransMatrix(energy,:,:,:)=TransE;
    end
    clear DecaySys;
    TransMatrix=reshape(TransMatrix,numel(SpeSet)*DetY,PixY*PixX);
    TransSys{angle} = TransMatrix;
    eval(sprintf('TransSysAngle%d = TransMatrix;',angle));
    eval(sprintf('save TransMatrix/TransSysAngle%d TransSysAngle%d;', angle, angle));
    eval(sprintf('clear TransSysAngle%d;', angle));
end
fprintf('calculate transmission model done....\n======================================\n');

%=============================
% Load geosys matrix and panelData, save in cell structure.
%=============================
fprintf('Load geosys matrix ..\n');
GeoSys = cell(numAngle,1);
for angle = 1:numAngle 
    GeoName = sprintf('GeoSysAngle%d.mat',angle);
    load(GeoName);
    eval(sprintf('GeoSys{%d} = GeoSysAngle%d;',angle,angle));
    eval(sprintf('clear GeoSysAngle%d;', angle));
end

fprintf('Load raw data ...\n');
PanelData = cell(numAngle,1);
for angle = 1:numAngle 
    DataName = sprintf('PanelDataAngle%d.mat',angle);
    load(DataName);
    eval(sprintf('PanelData{%d} = PanelDataAngle%d;',angle,angle));
    PanelData{angle}=PanelData{angle};
    eval(sprintf('clear PanelDataAngle%d;', angle));
end
load DiffSys4W;

%=============================
% Reconstruction
%=============================
fprintf('EMTV reconstruction....\n======================================\n');
DetMask=ones(DetZ,DetY);
DetMask(31:34,:)=0;

sysParameter.AngleSet = AngleSet;
sysParameter.SourcePos = SourcePos;

%物体参数
sysParameter.PixX = PixX;
sysParameter.PixY = PixY;
sysParameter.PixelSize = PixelSize;

%探测器参数
sysParameter.DetPos = DetPos;
sysParameter.DetY = DetY;
sysParameter.DetZ = DetZ;
sysParameter.offsetY = offsetY;
sysParameter.offsetZ = offsetZ;
sysParameter.DetSize = DetSize;
sysParameter.DetMask = DetMask;

%谱参数
sysParameter.SpeSet = SpeSet;
sysParameter.MTSet= MTSet;

%重建参数
parameter.subEMStep=1;
parameter.lambda1 = 1e-2;
parameter.betaTVSpe = 1e-2;
parameter.lambdaY = 6e-3;
parameter.lambdaX = parameter.lambdaY;
parameter.betaTVSpaVlambdaY=1e-2;
parameter.betaTVSpaY = parameter.betaTVSpaVlambdaY*parameter.lambdaY;
parameter.betaTVSpaX = parameter.betaTVSpaY;
parameter.totalSPStep = 25000;
parameter.savePath = 'ResultsOSEMTVSimple_15ViewTV6e-3';
ReconImage = OSEMTV3DXRDTCalculateLossGS(PanelData, sysParameter, TransSys, GeoSys, DiffSys, ResponseMatrix,parameter);

save ReconImage ReconImage;
clear;
