function out = io_loadspec_mat(filename, SVS_mask_filename, B0_map_filename)
    %% Read Data
    [csi, ReadInInfo, Par] = load_mat(filename);
    % mask = load_mask_osprey(filename); %reading SVS mask for Osprey
    mask = logical(load_nifti(SVS_mask_filename));
    B0_map = load_nifti(B0_map_filename);
    [fids, n_selected_voxels] = load_reshape_fids(csi, mask);

    %% Load and set parameters
    [txfrq, B0, dwelltime, spectralwidth, centerFreq, TE, TR] = load_parameters(ReadInInfo);
    geometry = load_geometry(ReadInInfo, Par);
    dims = set_dims();
    sz = size(fids);
    n_fid_points = sz(1);
    t = 0:dwelltime:((n_fid_points - 1) * dwelltime); % time points for fids
    ppm = calculate_ppm(txfrq, n_fid_points, spectralwidth, centerFreq); % ppm points for specs
    averages = n_selected_voxels; %number of voxel in the region
    rawAverages = n_selected_voxels; %number of voxel in the region
    subspecs = 1;
    rawSubspecs = 1;
    seq = 'MRSI_Vienna';
    date = '';
    leftshift = 0;

    %% Corrections
    %[fids] = linear_baseline_fitting(fids,dims);
    fids = frequency_correction(fids, B0_map, dims, mask, spectralwidth, txfrq);

    %% Calculate spectrum from fids
    specs = fftshift(fft(fids, [], dims.t), dims.t);

    %% Filling in data structure for Osprey
    out.fids = fids;
    out.specs = specs;
    out.sz = sz;
    out.ppm = ppm;
    out.t = t;
    out.spectralwidth = spectralwidth;
    out.dwelltime = dwelltime;
    out.txfrq = txfrq;
    out.date = date;
    out.dims = dims;
    out.Bo = B0;
    out.averages = averages;
    out.rawAverages = rawAverages;
    out.subspecs = subspecs;
    out.rawSubspecs = rawSubspecs;
    out.seq = seq;
    out.te = TE;
    out.tr = TR;
    out.pointsToLeftshift = leftshift;
    out.centerFreq = centerFreq;
    out.geometry = geometry;

    % Sequence flags
    out.flags.isUnEdited = 1;
    out.flags.isMEGA = 0;
    out.flags.isHERMES = 0;
    out.flags.isHERCULES = 0;
    out.flags.isPRIAM = 0;
    out.flags.isMRSI = 0;

    %FILLING IN THE FLAGS
    out.flags.writtentostruct = 1;
    out.flags.gotparams = 1;
    out.flags.leftshifted = 0;
    out.flags.filtered = 0;
    out.flags.zeropadded = 0;
    out.flags.freqcorrected = 0;
    out.flags.phasecorrected = 0;
    out.flags.averaged = 0;
    out.flags.addedrcvrs = 1;
    out.flags.subtracted = 0;
    out.flags.writtentotext = 0;
    out.flags.downsampled = 0;
    out.flags.coilsCombined = 1;

    if out.dims.subSpecs == 0
        out.flags.isFourSteps = 0;
    else
        out.flags.isFourSteps = (out.sz(out.dims.subSpecs) == 4);
    end

end

function [csi, ReadInInfo, Par] = load_mat(filename)
    csi = load(filename, 'csi').csi;
    ReadInInfo = load(filename, 'ReadInInfo').ReadInInfo;
    Par = load(filename, 'Par').Par;
end

function [txfrq, B0, dwelltime, spectralwidth, centerFreq, TE, TR] = load_parameters(ReadInInfo)
    txfrq = ReadInInfo.Par.LarmorFreq / 1e6;
    B0 = txfrq / 42.577;

    if abs(7 - B0) < 0.5
        dwelltime = 2 * ReadInInfo.Par.Dwelltime;
        TE = 1.3;
        TR = 460;
    else
        dwelltime = ReadInInfo.Par.Dwelltime;
        TE = 0.8;
        TR = 950;
    end

    centerFreq = 4.65; % Siemens data assumes the center frequency to be 4.7 ppm: % 4.65 from Vienna script
    spectralwidth = (1e9 / dwelltime);
end

function ppm = calculate_ppm(txfrq, n_fid_points, spectralwidth, centerFreq)
    frequency_step = spectralwidth / n_fid_points;
    frequency_start = -spectralwidth / 2 + frequency_step / 2;
    frequency_end = spectralwidth / 2 - frequency_step / 2;
    f = frequency_start:frequency_step:frequency_end;
    ppm = f / txfrq;
    ppm = ppm + centerFreq;
end

function geometry = load_geometry(ReadInInfo, Par)
    geometry.size.VoI_RoFOV = ReadInInfo.Par.ReadVOI; % Voxel size in readout direction [mm]
    geometry.size.VoI_PeFOV = ReadInInfo.Par.PhaseVOI; % Voxel size in phase encoding direction [mm]
    geometry.size.VoIThickness = ReadInInfo.Par.SliceVOI; % Voxel size in slice selection direction [mm]
    geometry.pos.PosCor = Par.CSI.PosVOI_Cor; % Coronal coordinate of voxel [mm]
    geometry.pos.PosSag = Par.CSI.PosVOI_Sag; % Sagittal coordinate of voxel [mm]
    geometry.pos.PosTra = Par.CSI.PosVOI_Tra; % Transversal coordinate of voxel [mm]
    geometry.pos.TablePosSag = 0; % Sagittal table position [mm].
    geometry.pos.TablePosCor = 0; % Coronal table position [mm].
    geometry.pos.TablePosTra = 0; % Transversal table position [mm].
    geometry.rot.VoI_InPlaneRot = ReadInInfo.Par.InPlaneRotation; % Voxel rotation in plane
    geometry.rot.NormCor = ReadInInfo.Par.SliceNormalVector_x; % Coronal component of normal vector of voxel.
    geometry.rot.NormSag = ReadInInfo.Par.SliceNormalVector_y; % Sagittal component of normal vector of voxel.
    geometry.rot.NormTra = ReadInInfo.Par.SliceNormalVector_z; % Transversal component of normal vector of voxel.
end

function [fids, n_selected_voxels] = load_reshape_fids(csi, mask)
    %mask(:,:,:) = 0;
    %mask(33:34,25:26,12:13) = 1;
    n_selected_voxels = sum(mask, 'all');
    n_time_points = size(csi, 4);
    mask_4D = repmat(mask, 1, 1, 1, n_time_points);
    selected_voxels = csi(mask_4D);
    reshaped_voxels = reshape(selected_voxels, n_selected_voxels, n_time_points);
    fids = transpose(reshaped_voxels);
    fids = double(fids);
    fids = fids ./ fids(3, :) .* abs(fids(3, :)); %this is the phase correction, using third point in fids.
    %fids = flip(fids,1); %it might need to be flip for Osprey.
end

function dims = set_dims()
    dims.t = 1;
    dims.coils = 0; %this can be 0.
    dims.averages = 2; %this can be 0 if crash.
    dims.subSpecs = 0;
    dims.extras = 0;
end

function mask = load_mask_osprey(CSI_filename)
    dirname = fileparts(CSI_filename);
    new_mask_name = 'SVS_mask.nii.gz';
    mask_filename = fullfile(dirname, new_mask_name);
    nii = nii_tool('load', mask_filename);
    mask = nii.img;
end

function img = load_nifti(filename)
    B0_map = nii_tool('load', filename);
    img = B0_map.img;
end
