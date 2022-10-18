%--------------------------------------------------------------------------
% * Kaichao Liang, 2022.01.11
% * Generate coded aperture forward projection PanelData = CodedSys*Object
% online.
% * Target domain: [Energy，DetZ，DetY], source domain: [PixY, PixX].
%--------------------------------------------------------------------------
function PanelData = CodedXRDTForward(objectXRD,PixX,PixY, ...
    DetY,DetZ,SpeSet,RotAngle,DiffMatrix,GeoMatrix,TransMatrix, ResponseMatrix)
    %%-----------------------------paramter--------------------------------
    % objectXRD: the Raylscattering distribution, the size [numel(MTSet),PixY, PixX]
    % PixX: the number of object pixels in X direction.
    % PixY: the number of object pixels in Y direction.
    % DetY: the number of detector bins in Y direction.
    % DetZ: the number of detector bins in Z direction.
    % SpeSet: the discrete vector for spectrum.
    % RotAngle: the object rotate angle, [0, 2*pi).
    % DiffMatrix: diffraction model restored in cell [numel(SpeSet)*DetZ*DetY,numel(MTSet)*PixY*PixX].
    % GeoMatrix:Pre-generated Coded aperture geometry factor [DetZ*DetY,PixY*PixX].
    % TransMatrix:Pre-generated transmission spectrum
    % [numel(SpeSet)*DetY,PixY*PixX].assume the same attenuation in DetZ
    % dimension.
    % ResponseMatrix: Detector energy response [Energy, Energy]
    %%---------------------------------------------------------------------
    tic;
    PanelData = zeros(numel(SpeSet),DetZ,DetY);
    for x = 1:PixX
        for y = 1:PixY
            
            %%calculate the corresponding position without rotation
            rotx = (x-PixX/2-0.5)*cos(RotAngle)-(y-PixY/2-0.5)*sin(RotAngle)+PixX/2+0.5;
            roty = (y-PixY/2-0.5)*cos(RotAngle)+(x-PixX/2-0.5)*sin(RotAngle)+PixY/2+0.5;
            rotx = floor(rotx+0.5); roty = floor(roty+0.5);
            if(~(rotx>=1 && rotx<=PixX &&roty>=1&&roty<=PixY))
                continue;
            end
            
            DiffSysPix = DiffMatrix{roty,rotx};
            PanelDataPix = reshape(DiffSysPix*objectXRD(:,y,x),numel(SpeSet),DetZ/2,DetY);
            PanelDataPixF = flip(PanelDataPix,2); 
            PanelDataPix = cat(2, PanelDataPix,PanelDataPixF);
            
            %%Code geometry
            geoFactor = GeoMatrix(:,y+PixY*(x-1));
            geoFactor = repmat(reshape(geoFactor,1,DetZ,DetY),numel(SpeSet),1,1);
            PanelDataPix = PanelDataPix.*geoFactor;
            
            %%Transmission spectrum
            speFactor = TransMatrix(:,y+PixY*(x-1));
            speFactor = repmat(reshape(speFactor,numel(SpeSet),1,DetY),1,DetZ,1);
            PanelDataPix = PanelDataPix.*speFactor;
            
            PanelData=PanelData+PanelDataPix;
        end
    end
    PanelData = reshape(PanelData,numel(SpeSet),DetZ*DetY);
    PanelData = ResponseMatrix*PanelData;
    PanelData = reshape(PanelData,numel(SpeSet),DetZ,DetY);
    toc;
end

    

    
