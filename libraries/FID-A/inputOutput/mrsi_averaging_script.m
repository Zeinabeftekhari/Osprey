filename = 'D:\Osprey_new\mrsi_data\CombinedCSI.mat'; %set to your combined.mat path
out_dir = 'D:\Osprey_new\mrsi_data\output'; %set your output path
basis_path = 'D:\Osprey_new\mrsi_data\fid_1.300000ms.basis'; %it is for 7T
LCM_ProgramPath = 'D:\Ospery\osprey-develop\libraries\LCModel\win\win10\LCModel_win_win10\LCModel.exe';
LCM_ControlFile = 'D:\Osprey_new\mrsi_data\LCModel_Control_Template.m';

%% Load Data, apply mask, phase correct, remove baseline
out = io_loadspec_mat(filename);

%% Voxel averaging
Data.csi = mean(out.fids, 2);

%% Reassign Info for LCM
Paths.out_dir = out_dir;
Paths.basis_file = {basis_path};
Paths.LCM_ProgramPath = LCM_ProgramPath;
MetaInfo.DimNames = {'x','y','z'};
MetaInfo.pat_name = 'anonymous'; % patient_name
MetaInfo.LarmorFreq = double(out.txfrq*1e6);
MetaInfo.dwelltime = 2 * double(out.dwelltime);

%% Write LCM Files
control_file = Write_LCM_files(Data, Paths, MetaInfo, LCM_ControlFile); % data is averaged fids

%% Run LCM
command = sprintf('"%s" < "%s"', LCM_ProgramPath, control_file);
system(command)
