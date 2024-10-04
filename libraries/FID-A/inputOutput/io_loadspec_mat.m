function out = io_loadspec_mat(filename, SVS_mask_filename, B0_map_filename)
    %% Read Data
    [csi, ReadInInfo, Par] = load_mat(filename);
    % mask = load_mask_osprey(filename); %reading SVS mask for Osprey
    mask = logical(load_nifti(SVS_mask_filename));
    B0_map = load_nifti(B0_map_filename);
    fids = load_reshape_fids(csi, mask);

    %% Load and set parameters
    params = load_parameters(ReadInInfo, Par, fids);

    %% Corrections
    %[fids] = linear_baseline_fitting(fids,dims);
    frequency_correction_on_fids = true;
    fids = frequency_correction(fids, B0_map, mask, params, frequency_correction_on_fids); %Frequency correction for each voxel before averaging

    %% Calculate spectrum from fids
    specs = fftshift(fft(fids, [], params.dims.t), params.dims.t);

    %% Filling in output structure
    out = params;
    out.fids = fids;
    out.specs = specs;
    out.flags = get_osprey_flags();
end

function fids = load_reshape_fids(csi, mask)
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

function [csi, ReadInInfo, Par] = load_mat(filename)
    csi = load(filename, 'csi').csi;
    ReadInInfo = load(filename, 'ReadInInfo').ReadInInfo;
    Par = load(filename, 'Par').Par;
end

function img = load_nifti(filename)
    B0_map = nii_tool('load', filename);
    img = B0_map.img;
end

function mask = load_mask_osprey(CSI_filename)
    dirname = fileparts(CSI_filename);
    new_mask_name = 'SVS_mask.nii.gz';
    mask_filename = fullfile(dirname, new_mask_name);
    mask = load_nifti(mask_filename);
end

function ppm = calculate_ppm(txfrq, n_fid_points, spectralwidth, centerFreq)
    frequency_step = spectralwidth / n_fid_points;
    frequency_start = -spectralwidth / 2 + frequency_step / 2;
    frequency_end = spectralwidth / 2 - frequency_step / 2;
    f = frequency_start:frequency_step:frequency_end;
    ppm = f / txfrq;
    ppm = ppm + centerFreq;
    ppm = flip(ppm);
end

function params = load_parameters(ReadInInfo, Par, fids)
    n_fid_points = size(fids, 1);
    n_selected_voxels = size(fids, 2);
    params.txfrq = ReadInInfo.Par.LarmorFreq / 1e6;
    params.Bo = params.txfrq / 42.577;

    if abs(7 - params.Bo) < 0.5
        params.dwelltime = 2 * ReadInInfo.Par.Dwelltime;
        params.TE = 1.3;
        params.TR = 460;
    else
        params.dwelltime = ReadInInfo.Par.Dwelltime;
        params.TE = 0.8;
        params.TR = 950;
    end

    params.centerFreq = 4.65; % Siemens data assumes the center frequency to be 4.7 ppm: % 4.65 from Vienna script
    params.spectralwidth = (1e9 / params.dwelltime);
    params.geometry = load_geometry(ReadInInfo, Par);
    params.dims = set_dims();
    params.sz = size(fids);
    params.t = 0:params.dwelltime:((n_fid_points - 1) * params.dwelltime); % time points for fids
    params.ppm = calculate_ppm(params.txfrq, n_fid_points, params.spectralwidth, params.centerFreq); % ppm points for specs
    params.averages = n_selected_voxels; %number of voxel in the region
    params.rawAverages = n_selected_voxels; %number of voxel in the region
    params.subspecs = 1;
    params.rawSubspecs = 1;
    params.seq = 'MRSI_Vienna';
    params.date = '';
    params.pointsToLeftshift = 0;
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

function flags = get_osprey_flags()
    % Sequence flags
    flags.isUnEdited = 1;
    flags.isMEGA = 0;
    flags.isHERMES = 0;
    flags.isHERCULES = 0;
    flags.isPRIAM = 0;
    flags.isMRSI = 0;

    %FILLING IN THE FLAGS
    flags.writtentostruct = 1;
    flags.gotparams = 1;
    flags.leftshifted = 0;
    flags.filtered = 0;
    flags.zeropadded = 0;
    flags.freqcorrected = 0;
    flags.phasecorrected = 0;
    flags.averaged = 0;
    flags.addedrcvrs = 1;
    flags.subtracted = 0;
    flags.writtentotext = 0;
    flags.downsampled = 0;
    flags.coilsCombined = 1;
    flags.isFourSteps = 0;
end

function dims = set_dims()
    dims.t = 1;
    dims.coils = 0; %this can be 0.
    dims.averages = 2; %this can be 0 if crash.
    dims.subSpecs = 0;
    dims.extras = 0;
end
