%--------------------------------------------------------------------------
% * Kaichao Liang, 2022.01.11
% * Generate coded projection backward transpotation  online.
% * Object = CodedSys'*PanelData
% * Source domain: [PixY, PixX],target domain: [Energy，DetZ，DetY]
%--------------------------------------------------------------------------
function objectXRD = CodedXRDTBackward(panelData,PixX,PixY,...
    DetY,DetZ,SpeSet,RotAngle,DiffMatrix,GeoMatrix,TransMatrix, ResponseMatrix)
    %%-----------------------------paramter--------------------------------
    % panelData: the Raylscattering raw data, the size [numel(SpeSet),DetZ, DetY]
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
    ResponseMatrix=ResponseMatrix';
    panelData = reshape(panelData,numel(SpeSet),DetZ*DetY);
    panelData = ResponseMatrix*panelData;
    panelData = reshape(panelData,numel(SpeSet),DetZ,DetY);
    
    for x=1:PixX
        for y=1:PixY
            PanelDataPix=panelData;
            
            %%Code geometry
            geoFactor = GeoMatrix(:,y+PixY*(x-1));
            geoFactor = repmat(reshape(geoFactor,1,DetZ,DetY),numel(SpeSet),1,1);
            PanelDataPix = PanelDataPix.*geoFactor;
            
            %%Transmission spectrum
            speFactor = TransMatrix(:,y+PixY*(x-1));
            speFactor = repmat(reshape(speFactor,numel(SpeSet),1,DetY),1,DetZ,1);
            PanelDataPix = PanelDataPix.*speFactor;
            
            rotx = (x-PixX/2-0.5)*cos(RotAngle)-(y-PixY/2-0.5)*sin(RotAngle)+PixX/2+0.5;
            roty = (y-PixY/2-0.5)*cos(RotAngle)+(x-PixX/2-0.5)*sin(RotAngle)+PixY/2+0.5;
            rotx = floor(rotx+0.5); roty = floor(roty+0.5);
            if(~(rotx>=1 && rotx<=PixX &&roty>=1&&roty<=PixY))
                continue;
            end
            DiffSysPix = DiffMatrix{roty,rotx};
            
            PanelDataPixUp = PanelDataPix(:,1:DetZ/2,:);
            PanelDataPixDown = flip(PanelDataPix(:,DetZ/2+1:end,:),2);
            
            XRDPattern = DiffSysPix'*PanelDataPixUp(:);
            XRDPattern = XRDPattern+DiffSysPix'*PanelDataPixDown(:);
            
            objectXRD(:,y,x) = XRDPattern;
 
        end
    end
    toc;
end

    

    
