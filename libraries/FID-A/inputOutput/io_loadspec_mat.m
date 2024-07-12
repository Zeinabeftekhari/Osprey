function out=io_loadspec_mat(filename);
csi = load(filename, 'csi').csi;
ReadInInfo = load(filename , 'ReadInInfo').ReadInInfo;
Par = load(filename , 'Par').Par;

fids = load_reshape_fids(csi);

dims.t = 1;
dims.coils = 0; %this can be 0.
dims.averages = 0; %this can be 0 if crash.
dims.subSpecs = 0;
dims.extras = 0;

specs = fftshift(fft(fids,[],dims.t),dims.t);
sz=size(fids);
B0 = txfrq/42.577;
ppm = load_ppm(ReadInInfo);
t=[0:dwelltime:(sz(1)-1)*dwelltime];
averages = ReadInInfo.Par.nAve; %number of voxel in the region
rawAverages = ReadInInfo.Par.nAve; %number of voxel in the region
subspecs = 1;
rawSubspecs = 1;

seq = 'MRSI_Vienna';

date = '';

if abs(7-B0) < 0.5
    TE = 1.3;
else
    TE = 0.8;
end
if abs(7-B0) < 0.5
    TR = 460;
else
    TR = 950;
end

leftshift = 0; %i am not sure, find it later.

geometry = load_geometry(ReadInInfo,Par);

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
end

function ppm = load_ppm(ReadInInfo)
txfrq = ReadInInfo.Par.LarmorFreq/1e6;
dwelltime = ReadInInfo.Par.Dwelltime;
B0 = txfrq/42.577;
spectralwidth = 1e9 / dwelltime;
frequency_step = spectralwidth / sz(1);
frequency_start = -spectralwidth / 2 + frequency_step / 2;
frequency_end = spectralwidth / 2 - frequency_step / 2;
f = frequency_start:frequency_step:frequency_end;
ppm = f / txfrq;
% Siemens data assumes the center frequency to be 4.7 ppm: % 4.65 from Vienna script
centerFreq = 4.65;
ppm=ppm + centerFreq;
end

function geometry = load_geometry(ReadInInfo,Par)
geometry.size.VoI_RoFOV     = ReadInInfo.Par.ReadVOI; % Voxel size in readout direction [mm]
geometry.size.VoI_PeFOV     = ReadInInfo.Par.PhaseVOI; % Voxel size in phase encoding direction [mm]
geometry.size.VoIThickness  = ReadInInfo.Par.SliceVOI; % Voxel size in slice selection direction [mm]
geometry.pos.PosCor         = Par.CSI.PosVOI_Cor; % Coronal coordinate of voxel [mm]
geometry.pos.PosSag         = Par.CSI.PosVOI_Sag; % Sagittal coordinate of voxel [mm]
geometry.pos.PosTra         = Par.CSI.PosVOI_Tra; % Transversal coordinate of voxel [mm]
geometry.pos.TablePosSag    = 0; % Sagittal table position [mm].
geometry.pos.TablePosCor    = 0; % Coronal table position [mm].
geometry.pos.TablePosTra    = 0; % Transversal table position [mm].
geometry.rot.VoI_InPlaneRot = ReadInInfo.Par.InPlaneRotation; % Voxel rotation in plane
geometry.rot.NormCor        = ReadInInfo.Par.SliceNormalVector_x; % Coronal component of normal vector of voxel.
geometry.rot.NormSag        = ReadInInfo.Par.SliceNormalVector_y; % Sagittal component of normal vector of voxel.
geometry.rot.NormTra        = ReadInInfo.Par.SliceNormalVector_z; % Transversal component of normal vector of voxel.
end

function fids = load_reshape_fids(csi)
fids = csi(30,30,15,:); %mask information later on will add here,
reshape_size = size(fids);
new_size = [reshape_size(4) 1 1];
fids = reshape(fids, new_size);
end