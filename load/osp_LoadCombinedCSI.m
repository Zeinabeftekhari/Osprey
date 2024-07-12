function [MRSCont] = osp_LoadCombinedCSI(MRSCont)
%% [MRSCont] = osp_LoadCombinedCSI(MRSCont)
%   This function reads raw (un-combined) data from MRSI Vienna files.
%
%   USAGE:
%       [MRSCont] = osp_LoadCombinedCSI(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Korbinian Eckstein and Zeinab Eftekhari
%       z.eftekhari@uq.edu.au
%
%   HISTORY:
%       2024-06-14: First version of the code.
% in this script I changed the io_loadspec_twix to io_loadspect_mat and changed the te to 1.3 (mrsi te)
% Close any remaining open figures
close all;

%loading combinedcsi.mat file,
%filename_combinedCSI = MRSCont.file{1,1};
%CSI = load(filename_combinedCSI, 'csi');
raw  = io_loadspec_mat(MRSCont.files{1,1});
%raw  = op_leftshift(raw,raw.pointsToLeftshift); % it changed all ppm
%values to 4.something
MRSCont.raw{1,1} = raw;
MRSCont.flags.coilsCombined     = 1;
% Add NIfTI-MRS information
%raw   = osp_add_nii_mrs_field(raw,MRSCont.ver.Osp);