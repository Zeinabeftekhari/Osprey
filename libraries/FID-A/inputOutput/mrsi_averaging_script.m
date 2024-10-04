main_path = "/Volumes/MRSHB2021-Q4401/Data-3t-7t";

for subject = ["sub-002" "sub-006" "sub-007" "sub-008" "sub-009"]
    for session = ["ses01" "ses02"]
        for location = ["ll" "ul"]
            for field = ["3T" "7T"]
                MRSI_filename = fullfile(main_path, field, subject, "ProcessedData", "MRSI_" + session, "CombinedCSI.mat"); % set to your combined.mat path
                SVS_mask = fullfile(main_path, "analysis/MRSI/SVS_MRSI/SVS_MRSI/output/" + subject + "-" + session + "_" + field, "mri", subject + "_" + session + "_" + location + "_mask_" + field + "_register.nii.gz"); % set to your mask path
                B0_map =  fullfile(main_path, field, subject, "ProcessedData", "MRSI_" + session, "maps/Extra/shift_map.nii");
                out_dir = fullfile(main_path, "analysis/MRSI/SVS_MRSI/MRSI_averaged"); % set your output path

                if  strcmp("3T", field)
                    basis_path = "/Users/uqzeftek/Data/fid_0.800000ms.basis"; % it is for 3T
                    LCM_ControlFile = "/Users/uqzeftek/Data/LCModel_Control_Template_3T.m";
                else
                    basis_path = "/Users/uqzeftek/Data/fid_1.300000ms.basis"; % it is for 7T
                    LCM_ControlFile = "/Users/uqzeftek/Data/LCModel_Control_Template.m";
                end
                LCM_ProgramPath = "/Users/uqzeftek/Github/LCModel/source/lcmodel";
                pat_name = subject + '_' + session + '_' + location + '_' + field;

                mrsi_averaging(MRSI_filename, SVS_mask, B0_map, basis_path, LCM_ControlFile, LCM_ProgramPath, out_dir, pat_name);
            end
        end
    end
end
