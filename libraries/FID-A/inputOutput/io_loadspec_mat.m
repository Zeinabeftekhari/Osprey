function out = io_loadspec_mat(filename)
[csi, ReadInInfo, Par] = load_mat(filename);
mask = load_mask(filename);


dims = set_dims();

[txfrq, B0, dwelltime, spectralwidth, centerFreq, TE, TR] = load_parameters(ReadInInfo);

[fids, number_selected_voxels] = load_reshape_fids(csi, mask);
%[fids] = linear_baseline_fitting(fids,dims);
sz = size(fids);
t = 0:dwelltime:((sz(1)-1)*dwelltime); % time points for fids
specs = fftshift(fft(fids,[],dims.t),dims.t);
ppm = calculate_ppm(txfrq,sz,spectralwidth,centerFreq); % ppm points for specs
averages = number_selected_voxels; %number of voxel in the region
rawAverages = number_selected_voxels; %number of voxel in the region
subspecs = 1;
rawSubspecs = 1;

seq = 'MRSI_Vienna';
date = '';

leftshift = 0;

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
out.centerFreq=centerFreq;
out.geometry=geometry;

% Sequence flags
out.flags.isUnEdited = 1;
out.flags.isMEGA = 0;
out.flags.isHERMES = 0;
out.flags.isHERCULES = 0;
out.flags.isPRIAM = 0;
out.flags.isMRSI = 0;

%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=0;
out.flags.addedrcvrs=1;
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

function [csi,ReadInInfo,Par] = load_mat(filename)
csi = load(filename, 'csi').csi;
ReadInInfo = load(filename , 'ReadInInfo').ReadInInfo;
Par = load(filename , 'Par').Par;
end

function [txfrq, B0, dwelltime, spectralwidth, centerFreq, TE, TR] = load_parameters(ReadInInfo)
txfrq = ReadInInfo.Par.LarmorFreq/1e6;
B0 = txfrq/42.577;
dwelltime = ReadInInfo.Par.Dwelltime;
centerFreq = 4.65; % Siemens data assumes the center frequency to be 4.7 ppm: % 4.65 from Vienna script
spectralwidth = (1e9 / dwelltime)/2;
[TE,TR] = get_TE_TR(B0);
end

function [TE,TR] = get_TE_TR(B0)
if abs(7-B0) < 0.5
    TE = 1.3;
    TR = 460;
else
    TE = 0.8;
    TR = 950;
end
end

function ppm = calculate_ppm(txfrq,sz,spectralwidth,centerFreq)
frequency_step = spectralwidth / sz(1);
frequency_start = -spectralwidth / 2 + frequency_step / 2;
frequency_end = spectralwidth / 2 - frequency_step / 2;
f = frequency_start:frequency_step:frequency_end;
ppm = f / txfrq;
ppm = ppm + centerFreq;
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

function [fids, number_selected_voxels] = load_reshape_fids(csi, mask)
mask = logical(mask);
%mask(:,:,:) = 0;  
%mask(33:34,25:26,12:13) = 1;
number_selected_voxels = sum(mask,'all');
number_timepoints = size(csi,4);
mask_4D = repmat(mask,1,1,1,number_timepoints);
selected_voxels = csi(mask_4D);
reshaped_voxels = reshape(selected_voxels,number_selected_voxels,number_timepoints);
fids = transpose(reshaped_voxels);
fids = double(fids);
fids = fids ./ fids(3,:) .* abs(fids(3,:)); %this is the phase correction, using third point in fids. 
%fids = flip(fids,1); %it might need to be flip for Osprey.
end

function [fids] = linear_baseline_fitting(fids,dims)
specs = fftshift(fft(fids,[],dims.t),dims.t);
startmean = mean(specs(1:100,:)); 
endmean = mean(specs(end-99:end,:));
baseline = zeros(size(specs));
for voxel = 1:size(specs,dims.averages)
   baseline(:,voxel) = linspace(startmean(:,voxel),endmean(:,voxel),size(specs,dims.t));
end
specs = specs-baseline;
fids = ifft(fftshift(specs,dims.t),[],dims.t);
end

function dims = set_dims()
dims.t = 1;
dims.coils = 0; %this can be 0.
dims.averages = 2; %this can be 0 if crash.
dims.subSpecs = 0;
dims.extras = 0;
end

function mask = load_mask(CSI_filename)
dirname = fileparts(CSI_filename);
new_mask_name = 'SVS_mask.nii.gz';
mask_filename = fullfile(dirname, new_mask_name);
nii = nii_tool('load', mask_filename);
mask = nii.img;
end