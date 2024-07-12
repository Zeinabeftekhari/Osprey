function out=io_loadspec_mat(filename);
csi = load(filename, 'csi').csi;
ReadInInfo = load(filename , 'ReadInInfo').ReadInInfo;
Par = load(filename , 'Par').Par;
fids = csi(30,30,15,:); %mask information later on will add here,
reshape_size = size(fids);
new_size = [reshape_size(4) 1 1];
fids = reshape(fids, new_size); 
dims.t = 1;
specs = fftshift(fft(fids,[],dims.t),dims.t);
sz=size(fids);
%calculate ppm, t ,spectralwidth, dwelltime, txfrq, B0
txfrq = ReadInInfo.Par.LarmorFreq;
B0 = txfrq/42.577;
dwelltime = ReadInInfo.Par.Dwelltime;
spectralwidth=1/dwelltime;
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=f/(B0*42.577);
% Siemens data assumes the center frequency to be 4.7 ppm:
centerFreq = 4.7;
ppm=ppm + centerFreq;
t=[0:dwelltime:(sz(1)-1)*dwelltime];
date = '';
dims.coils = 2; %this can be 0.
dims.averages = 3; %this can be 0 if crash.
dims.subSpecs = 0;
dims.extra = 0;
averages = ReadInInfo.Par.nAve; %number of voxel in the region 
rawAverages = ReadInInfo.Par.nAve; %number of voxel in the region 
subspecs = 1;
rawSubspecs = 1; 
seq = 'MRSI_Vienna'; 
if abs(7-B0) < 0.5  
   TE = 1.3;
else
    TE = 0.8; %3T TR write it here
end
if abs(7-B0) < 0.5  
   TR = 460;
else
    TR = 950; %3T TR write it here
end
leftshift = 0; %i am not sure, find it later.
%these are stored in raw.geometry, I donot know how to store them in raw.geometry
geometry.size.VoI_RoFOV     = ReadInInfo.Par.ReadVOI; % Voxel size in readout direction [mm]
geometry.size.VoI_PeFOV     = ReadInInfo.Par.PhaseVOI; % Voxel size in phase encoding direction [mm]
geometry.size.VoIThickness  = ReadInInfo.Par.SliceVOI; % Voxel size in slice selection direction [mm]
geometry.pos.PosCor         = Par.CSI.PosVOI_Cor; % Coronal coordinate of voxel [mm]
geometry.pos.PosSag         = Par.CSI.PosVOI_Sag; % Sagittal coordinate of voxel [mm]
geometry.pos.PosTra         = Par.CSI.PosVOI_Tra; % Transversal coordinate of voxel [mm]
geometry.pos.TablePosSag    = 0; % Sagittal table position [mm]. donot know
geometry.pos.TablePosCor    = 0; % Coronal table position [mm]. donot know
geometry.pos.TablePosTra    = 0; % Transversal table position [mm]. donot know
geometry.rot.VoI_InPlaneRot = ReadInInfo.Par.InPlaneRotation; % Voxel rotation in plane
geometry.rot.NormCor        = ReadInInfo.Par.SliceNormalVector_x; % Coronal component of normal vector of voxel. donot know
geometry.rot.NormSag        = ReadInInfo.Par.SliceNormalVector_y; % Sagittal component of normal vector of voxel. donot know
geometry.rot.NormTra        = ReadInInfo.Par.SliceNormalVector_z; % Transversal component of normal vector of voxel. donot know
% nucleus = '1H';
% software = 'vd syngo MR E12'; %not know about 3T
% Manufacturer = 'Siemens'; 

%FILLING IN DATA STRUCTURE
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.ppm=ppm;
out.t=t;
out.spectralwidth=spectralwidth;
out.dwelltime=dwelltime;
out.txfrq=txfrq;
out.date=date;
out.dims=dims;
out.Bo=B0;
out.averages=averages;
out.rawAverages=rawAverages;
out.subspecs=subspecs;
out.rawSubspecs=rawSubspecs;
out.seq=seq;
out.te=TE;
out.tr=TR;
out.pointsToLeftshift=leftshift;
out.centerFreq = centerFreq;
out.geometry = geometry;

%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=0;
out.flags.addedrcvrs=0;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.coilsCombined=1;
if out.dims.subSpecs==0
    out.flags.isFourSteps=0;
else
    out.flags.isFourSteps=(out.sz(out.dims.subSpecs)==4);
end
