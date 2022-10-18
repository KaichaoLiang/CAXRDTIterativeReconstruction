%--------------------------------------------------------------------------
% * Kaichao Liang, 2021.1.18
% * Process the MC simulation data, the scatter data is restored in
% PanelData directory by each angle. The transmission data of different
% angles is merged together
% * Need: SimuData
%--------------------------------------------------------------------------
clear; clc;
mkdir('PanelData31View');
mkdir('TransData');
addpath('PanelData31View');
addpath('SimuData31ViewRect');
addpath('TransData');
addpath('GTData')

%---------------------------
% Parameter
%---------------------------
numAngle = 31;
NdetZ=64; NdetY=64; 
FullSpe = 120; SpeStart = 21; SpeEnd =85;

geoParameter.LSource2Center = 400;
geoParameter.LCenter2Det = 125;
geoParameter.Ndet = 44;
geoParameter.BinSize = 1.6;
geoParameter.NPhi = numAngle;
geoParameter.Nx = 54;
geoParameter.Ny = 54;
geoParameter.PixelSize = 1;

%---------------------------
% Pre-load data
%---------------------------

% Energy response matrix
load ResponseMatrix;
fprintf('Read bin data\n---------------------------\n');
DecayData = zeros(SpeEnd-SpeStart+1,NdetY,numAngle);
for angle = 1:numAngle
    fprintf('Read bin data angle%d\n',angle);
    DataName = sprintf('DataAngle%d.dat',angle-1);
    fid=fopen(DataName,'rb');
    PanelData = double(fread(fid,'int32'));
    fclose(fid);
    PanelData = reshape(PanelData,NdetY,NdetZ,FullSpe);
    PanelData = permute(PanelData,[3,2,1]);
    PanelData = PanelData(SpeStart:SpeEnd,:,:);
    
    %Add energy response
    PanelData=reshape(PanelData,65,NdetZ*NdetY);
    PanelData=ResponseMatrix*PanelData;
    PanelData=reshape(PanelData,65,NdetZ,NdetY);
    PanelData(:,31:34,:)=0;
    
    eval(sprintf('PanelDataAngle%d = PanelData;',angle));
    SaveName = sprintf('PanelData31View/PanelDataAngle%d.mat',angle);
    FileName = sprintf('PanelDataAngle%d',angle);
    save(SaveName,FileName);
    eval(sprintf('clear PanelData%d;',angle));
    
end
clear;clc;