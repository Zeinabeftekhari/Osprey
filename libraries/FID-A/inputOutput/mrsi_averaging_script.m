main_path = '/Volumes/MRSHB2021-Q4401/';

for subject = ["sub-002","sub-006","sub-007","sub-008","sub-009"]
    for session = ["ses01","ses02"]
        for location = ["ll","ul"]
            for field = ["3T","7T"]
                MRSI_filename = strcat (main_path,'Data-3t-7t/',field,'/',subject,'/ProcessedData/MRSI_',session,'/CombinedCSI.mat') %set to your combined.mat path
                SVS_mask = strcat (main_path, 'Data-3t-7t/analysis/MRSI/SVS_MRSI/SVS_MRSI/output/',subject,'-',session,'_',field,'/mri/',subject,'_',session,'_',location,'_mask_',field,'_register.nii.gz') %set to your mask path
                B0_map =  strcat (main_path,'Data-3t-7t/',field,'/',subject,'/ProcessedData/MRSI_',session,'/maps/Extra/shift_map.nii')
                out_dir = strcat (main_path, 'Data-3t-7t/analysis/MRSI/SVS_MRSI/MRSI_averaged') %set your output path

                if  strcmp('3T',field)
                    basis_path = '/Users/uqzeftek/Data/fid_0.800000ms.basis'; %it is for 3T
                    LCM_ControlFile = '/Users/uqzeftek/Data/LCModel_Control_Template_3T.m'

                else 
                    basis_path = '/Users/uqzeftek/Data/fid_1.300000ms.basis'; %it is for 7T
                    LCM_ControlFile = '/Users/uqzeftek/Data/LCModel_Control_Template.m';
                end
                LCM_ProgramPath = '/Users/uqzeftek/Github/LCModel/source/lcmodel';
                    
                %% Load Data, apply mask, phase correct, remove baseline
                out = io_loadspec_mat(MRSI_filename, SVS_mask, B0_map);
                
                %% Voxel averaging
                Data.csi = mean(out.fids, 2);
                %Data.csi = out.fids (:,1); % for one voxel without averaging

 
                %% Reassign Info for LCM
                Paths.out_dir = out_dir;
                Paths.basis_file = {basis_path};
                Paths.LCM_ProgramPath = LCM_ProgramPath;
                MetaInfo.DimNames = {'x','y','z'};
                MetaInfo.pat_name = char(strcat (subject,'_',session,'_',location,'_',field)); % lcmodel_output
                MetaInfo.LarmorFreq = double(out.txfrq*1e6);
                if  strcmp('3T',field)
                    MetaInfo.dwelltime = double(out.dwelltime);
                else
                    MetaInfo.dwelltime = double(out.dwelltime)*2;
                end              
                %% Write LCM Files
                control_file = Write_LCM_files(Data, Paths, MetaInfo, LCM_ControlFile); % data is averaged fids
                
                %% Run LCM
                command = sprintf('"%s" < "%s"', LCM_ProgramPath, control_file);
                system(command)

          end
        end
    end
end

