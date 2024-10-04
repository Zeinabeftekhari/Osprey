function mrsi_averaging(MRSI_filename, SVS_mask, B0_map, basis_path, LCM_ControlFile, LCM_ProgramPath, out_dir, pat_name)
    %% Load Data, apply mask, phase correct, remove baseline
    out = io_loadspec_mat(char(MRSI_filename), char(SVS_mask), char(B0_map));
                
    %% Voxel averaging
    Data.csi = mean(out.fids, 2);
    %Data.csi = out.fids (:,1); % for one voxel without averaging

    %% Reassign Info for LCM
    Paths.out_dir = char(out_dir);
    Paths.basis_file = {char(basis_path)};
    Paths.LCM_ProgramPath = char(LCM_ProgramPath);
    MetaInfo.DimNames = {'x','y','z'};
    MetaInfo.pat_name = char(pat_name); % lcmodel_output
    MetaInfo.LarmorFreq = double(out.txfrq*1e6);
    MetaInfo.dwelltime = double(out.dwelltime);

    %% Write LCM Files
    control_file = Write_LCM_files(Data, Paths, MetaInfo, LCM_ControlFile); % data is averaged fids
    
    %% Run LCM
    command = sprintf('"%s" < "%s"', LCM_ProgramPath, control_file);
    system(command)
end
