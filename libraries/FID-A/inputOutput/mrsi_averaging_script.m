filename = 'D:\Osprey_new\mrsi_data\CombinedCSI.mat'; %set to your combined.mat path
out_dir = 'D:\Osprey_new\mrsi_data\output'; %set your output path
basis_path = 'D:\Osprey_new\mrsi_data\fid_1.300000ms.basis'; %it is for 7T
LCM_ProgramPath = 'D:\Ospery\osprey-develop\libraries\LCModel\win\win10\LCModel_win_win10\LCModel.exe';
batchdir = 'D:\Osprey_new\mrsi_data\output';
out = io_loadspec_mat(filename);
Data.csi = mean(out.fids);
Paths.out_dir = out_dir;
Paths.basis_file = {basis_path};
Paths.LCM_ProgramPath = LCM_ProgramPath;
Paths.batchdir = batchdir;
MetaInfo.DimNames = {'x','y','z'};
MetaInfo.pat_name = 'anonymous'; %patient_name
MetaInfo.LarmorFreq = double(out.txfrq*1e6);
MetaInfo.dwelltime = 2 * double(out.dwelltime);

Write_LCM_files(Data,Paths,MetaInfo) %data is averaged fids 
command = sprintf()
