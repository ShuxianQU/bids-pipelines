%%% this is a script for using bids apps to process BOLD fMRI using offline
%%% recon such as those with amri_epi, where no dicoms are available. 
%%%
%%% created by Xiaoping on 10/21/2024

addpath(genpath('~/opt/utils'));
addpath(genpath('~/opt/bids'));
addpath(genpath('~/Documents/moco-recon'));
% add paths of '3dMPRAGEize' and 'AFNI' to the $PATH

% 049+             054-             037-             050-              039+
% session_20250130 session_20250206 session_20250213 session_20250220 session_20250227
% % sub-037's tissue segmentations and registrations show odd
% which are corrected by re-running smriprep part, but sub-050's remain wrong.
% % sub-054 and sub-037 still imperfection up to now.
subID= '039'; 
sessionName= 'session_20250227';

parentDir_containers= '~/opt/my_images';
parentDir_bids= '/home/range6-raid17/wulab/105t/';
path2bids= fullfile(parentDir_bids,'bids');
projectName= 'amri_epi';
path2bidsNewProj= fullfile(path2bids,filesep,projectName);

%
path2dicoms= ['/home/range6-raid17/wulab/105t/',sessionName]; 
% path2dicoms= ['/home/range6-raid4/xpwu-data/Siemens105T/',sessionName]; 
path2configJson= fullfile(path2bidsNewProj,'code','config.json');

%
path2raw= fullfile(path2dicoms,'raw'); % '/home/range6-raid4/xpwu-data/Siemens105T/session_20241024/raw';
path2offlineRecon= fullfile(path2raw,'result'); % '/home/range6-raid4/xpwu-data/Siemens105T/session_20241024/raw/result'; 
mid_offlineRecon= {'134' '138'}; %110 114
PELabel= {'AP' 'PA'};
mid_refScan= '132';
reconMode= 'nomoco'; %'moco': motion and B0 correction; 'nomoco': no rigid motion correction; 'noco': no motion or B0 correction
taskMode= 'EO'; %EO
switch reconMode
    case 'test'
        modifier= 'result_MoCo';
        taskID= 'restTest';
    case 'moco'
        modifier= 'result_MoCo';
        taskID= 'rest';
    case 'nomoco'
        modifier= 'result_gB0Co';
        taskID= 'restNomoco';
    case 'noco'
        modifier= 'result_noB0MoCo';
        taskID= 'restNoco';
    otherwise

end
taskID= [taskID,taskMode];


bSkipSmriprep= false;
bSkipFmriprep= false;
bSkipXcp= false;

% If use mp2rage, mprageize it for t1w.
bSkipT1wPrep= false; % true: t1w is ready, skipping mprageizing; false: to prep t1w by mprageizing mp2rage.
seid_MP2RAGEinv2='009'; % MR-SE<xx>
seid_MP2RAGEuni='010';

subDir= fullfile(path2bidsNewProj,['sub-',subID]);
%% step 1: dcm2bids
%% step 1.0:  'apptainer exec -e --containall -B /path/to/bids:/bids dcm2bids.sif dcm2bids_scaffold -o bids/new_project
% -- to generate project directory structure and configuration file templates 
% (NOTE: this should just be run once per project)
mycmd= ['apptainer exec -e --containall -B ', path2bids,':/bids ', ...
    parentDir_containers,filesep,'dcm2bids.sif dcm2bids_scaffold -o /bids/', projectName];
system(mycmd);
disp('-> dcm2bids_scaffold is done ...')

%% step 1.1:  'apptainer exec -e --containall -B /path/to/dicoms:/dicoms:ro -B /path/to/bids/new_project:/bids dcm2bids.sif dcm2bids_helper -o /bids -d /dicoms
% -- to convert DICOM to NIFTI and save them in '../bids/new_project/tmp_dcm2bids/helper/'
% (run this when processing a new subject)
mycmd= ['apptainer exec -e --containall -B ', path2bidsNewProj,':/bids -B ', ...
    path2dicoms,':/dicoms:ro ', ...
    parentDir_containers,filesep,'dcm2bids.sif dcm2bids_helper -o /bids -d /dicoms ', ...
    '--force_dcm2bids'];
system(mycmd);
disp('-> dcm2bids_helper is done ...')

%% create a config.json and copy it to path/to/
%%% for offline recon, perhaps to use a fake bold to have a placeholder,
%%% which was tested to not work since no func dir was created when picking
%%% a gre dicom
% -- create configuration json file in '../bids/new_project/code/'
% try placeholdder for func and fmap but failed
clear config
config.descriptions = { 
    struct( ...
        "datatype", "anat", ...
        "suffix", "T1w", ...
        "criteria", ...
            struct( ...
                "SeriesDescription", "mp2rage*UNI*" ...
            ), ... 
        "sidecar_changes", ...
            struct( ...
                "ProtocolName", "T1" ...
            ) ... 
    ) ...
    % struct( ...
    %     "id", "task_rest", ...
    %     "datatype", "func", ...
    %     "suffix", "bold", ...
    %     "custom_entities", "task-rest", ...      
    %     "criteria", ...
    %         struct( ...
    %             "ProtocolName", "func_task-*", ...
    %             "ImageType", ["ORIG*", "PRIMARY", "M", "MB", "ND", "MOSAIC"] ...
    %         ) ...
    % ), ...
    % struct( ...
    %     "datatype", "fmap", ...
    %     "suffix", "fmap", ...
    %     "criteria", ...
    %         struct( ...
    %             "ProtocolName", "*field_mapping*" ...
    %         ), ...
    %     "sidecar_changes", ... 
    %         struct( ...
    %             "IntendedFor", "task_rest" ...
    %         ) ...
    % ) ...
};

bids.util.jsonencode(path2configJson, config);

% update README file if it's empty for empty files are not allowed in BIDS
path2Readme = [path2bidsNewProj,filesep,'README'];
if isfile(path2Readme)
    fileInfoReadme = dir(path2Readme);
    if fileInfoReadme.bytes == 0
        fid = fopen(path2Readme,'w');
        fprintf(fid, ['This project aims to study how motion robust ' ...
                      '3D EPI with amri_epi would help improve BOLD ' ...
                      'fMRI of the human brain at 10.5T.\n']);
        fclose(fid);
        disp('README file has been updated.');
    else
        % disp('')
    end
end

disp('-> config.json for dcm2bids is generated ...')

%% step 1.2:  'apptainer run -e --containall -B /path/to/dicoms:/dicoms:ro
% -B /path/to/config.json:/config.json:ro
% -B /path/to/bids/new_project:/bids
% dcm2bids.sif -o /bids -d /dicoms -c /config.json -p 001
% --auto_extract_entities --bids_validate
% (run this when processing a new subject)
mycmd= ['apptainer run -e --containall -B ', path2bidsNewProj,':/bids -B ', ...
    path2dicoms,':/dicoms:ro -B ', path2configJson, ':/config.json:ro ', ...
    parentDir_containers,filesep,'dcm2bids.sif --auto_extract_entities --bids_validate ',...
    '-o /bids -d /dicoms -c /config.json --force_dcm2bids -p ',subID];
system(mycmd);
disp('-> dcm2bids is done ...')

%% mprageize mp2rage for t1w.
% (run this when processing a new subject)
if bSkipT1wPrep
    disp('=> Skipping T1w prep...')
else
    disp('=> Starting T1w prep ...')
    helperDir= [path2bidsNewProj,filesep,'tmp_dcm2bids',filesep,'helper'];
    filename_inv2= [seid_MP2RAGEinv2,'*.nii*'];
    inv2= dir([helperDir,filesep,filename_inv2]);
    filename_uni= [seid_MP2RAGEuni,'*.nii*'];
    uni= dir([helperDir,filesep,filename_uni]);

    file_inv2= [inv2.folder,filesep,inv2.name];
    file_uni= [uni.folder,filesep,uni.name];
    %mycmd=['3dMPRAGEize -i ', file_inv2, ' -u ', file_uni];
    mycmd=['3dMPRAGEize -i ', file_inv2, ' -u ', file_uni, ' -r 1']; % need to used the rebiased version otherwise artifacts in T1w will result from fmriprep
    system(mycmd);
    disp('-> 3dMPRAGEize is done ...')

    %% copy the mprageized uni file
    %filename_uniCleaned= [seid_MP2RAGEuni,'*unbiased_clean.nii*'];
    filename_uniCleaned= [seid_MP2RAGEuni,'*rebiased_clean.nii*'];
    uniCleaned= dir(filename_uniCleaned);

    anatDir= [subDir,filesep,'anat'];
    filename_t1w= [anatDir,filesep,'sub-',subID,'_T1w.nii*'];
    t1w= dir(filename_t1w);
    if ~isempty(t1w)
        delete([t1w.folder,filesep,t1w.name])
    end

    %movefile(uniCleaned.name, [anatDir,filesep,'sub-',subID,'_T1w.nii']);
    copyfile(uniCleaned.name, [anatDir,filesep,'sub-',subID,'_T1w.nii']);

    disp('-> 3dMPRAGEized MP2RAGE is copied ...')
end

%% step 2: smriprep/fmriprep
%% step 2.0: smriprep
% (run this when processing a new subject) 
if bSkipSmriprep
    disp('=> Skipping smriprep...')
else
    disp('=> Starting smriprep...')
    mycmd= ['singularity run --cleanenv ', parentDir_containers,filesep,'smriprep_0.17.0.sif ', path2bidsNewProj, ' ',...
        path2bidsNewProj,filesep,'derivatives participant --participant-label ',subID, ...
        ' --fs-license-file ', parentDir_containers,filesep, 'fs-license.txt -w ',path2bids,'/workdir_smriprep',...
        ' --output-spaces MNI152NLin2009cAsym --fs-no-reconall --stop-on-first-crash'];
    %%% Test the command without '--fs-no-reconall' 
    %%%    - error, because of missing 'mritotal' of freesurfer in both 'smriprep-0.16.1.simg' and 'smriprep_0.17.0.sif'
    %%% Adding '--cifti-output' was found to result in misaligned tissue segmentation and unusable surface recon, etc.
    system(mycmd);
    disp('-> smriprep is done ...')
end
%%%
%% cp offline recon to path2bidsNewProj/sub-<subID>/func/
% (run this when processing a new subject / a new condition) 
funcDir= [subDir,filesep,'func'];
if ~exist(funcDir,'dir')
    mkdir(subDir,'func')
end

for i = 1:length(mid_offlineRecon)   
    filename_offlineRecon= [path2offlineRecon,filesep,'im',mid_offlineRecon{i},'*',modifier,'*.nii.gz'];
    offlineRecon= dir(filename_offlineRecon);

    file_src= [offlineRecon.folder,filesep,offlineRecon.name];
    file_des= [funcDir,filesep,'sub-',subID,'_task-',taskID,'_dir-',PELabel{i},'_bold.nii.gz'];
    if exist(file_des,'file')
        disp(['-> The file ',file_des,' exists ...'])
        disp('-> Skipping copying offline recon ...')
    else
        copyfile(file_src, file_des);
        disp(['-> Offline recon with ',PELabel{i},' direction is copied ...'])
    end

    %%% remember to create a .json file incluidng repetition time, PhaseEncodingDirection and EffectiveEchoSpacing.
    name_spec.modality = 'func';
    name_spec.suffix = 'bold';
    name_spec.ext = '.json';
    name_spec.entities = struct('sub', subID, ...
                                'task', taskID, ...
                                'dir',PELabel{i});
    bidsFile_func = bids.File(name_spec, 'use_schema', true);
    jsonName_func = fullfile(funcDir, bidsFile_func.filename);
    
    json_func.TaskName = taskID;
    json_func.RepetitionTime = 2.34;
    if PELabel{i}=='AP'
        json_func.PhaseEncodingDirection = "j";
    elseif PELabel{i}=='PA'
        json_func.PhaseEncodingDirection = "j-";
    end
    json_func.EffectiveEchoSpacing = 0.00018;
    json_func.EchoTime = 0.021;
    
    bids.util.mkdir(fileparts(jsonName_func));
    bids.util.jsonencode(jsonName_func,json_func);
            
    disp(['-> json file for the offline recon with ',PELabel{i},' direction is generated ...'])

end

% disp('=> NOTE: make sure the associated .json file is ready before you run fmriprep...')

%% Recon B0 maps and save to path2bidsNewProj/sub-<subID>/fmap/ and create associated .json file. 
% (run this just when processing a new subject -- 
% field maps are calcultaed from the reference scan by default) 
fmapDir= [subDir,filesep,'fmap'];
if ~exist(fmapDir,'dir')
    mkdir(subDir,'fmap')
end
    
filename_refScan = [path2raw,filesep,'meas*',mid_refScan,'*.dat'];
refScan= dir(filename_refScan);

path2refScan = [refScan.folder,filesep,refScan.name];

file_src_b0map = [path2raw,filesep,'mid',mid_offlineRecon{1},'_b0fieldmap.nii.gz'];
file_src_mag = [path2raw,filesep,'mid',mid_offlineRecon{1},'_magnitude.nii.gz'];

file_des_b0map = [fmapDir,filesep,'sub-',subID,'_fieldmap.nii.gz'];
file_des_mag = [fmapDir,filesep,'sub-',subID,'_magnitude.nii.gz'];

methodofFM='useOE4B0map';
if exist(file_des_b0map,'file')
    disp(['-> The file ',file_des_b0map,' exists...'])
    disp('-> Skipping generating B0 field maps...')
else
    if ismember(methodofFM,{'useAE4B0map', 'useOE4B0map', 'useREB40map'})      
        B0_fieldmap(path2refScan,mid_offlineRecon{1},methodofFM);    

    elseif ismember(methodofFM,'useSDE4B0map')
        numFM=14;
        B0_fieldmap(path2refScan,mid_offlineRecon{1}, ...
            methodofFM,[path2bidsNewProj,'/tmp_dcm2bids/helper'],numFM);

    else
        error('Wrong calculation method of field maps');
    end
    copyfile(file_src_b0map, file_des_b0map);
    copyfile(file_src_mag, file_des_mag);
    disp('-> dB0 map is generated and copied...')
end

%% create .json file for field map
% (repeat this section whenever modifying func to ensure fmap.json is updated)
fmapDir= [subDir,filesep,'fmap'];
name_spec.modality = 'fmap';
name_spec.suffix = 'fieldmap';
name_spec.ext = '.json';
name_spec.entities = struct('sub',subID);

bidsFile_b0map = bids.File(name_spec, 'use_schema', true);
jsonName_b0map = fullfile(fmapDir, bidsFile_b0map.filename);

% schema = bids.Schema;
% def = schema.get_definition(taskID);

filelist = dir(fullfile(funcDir,'*nii.gz'));
clear filenameFunc
for i=1:length(filelist)
    filenameFunc{i} = fullfile('func',filelist(i).name);
end
json_b0map.Units = 'Hz';
json_b0map.B0FieldIdentifier = 'b0map_fmap0';
json_b0map.IntendedFor = filenameFunc;

bids.util.mkdir(fileparts(jsonName_b0map));
bids.util.jsonencode(jsonName_b0map,json_b0map);

disp('-> json file for fmap is generated...')

%% step 2.1: fmriprep
% 'singularity run --cleanenv ~/my_containers/fmriprep.simg
% testData4bids_apps/ testData4bids_apps/output participant --participant-label 01 --fs-license-file ~/my_containers/fs-license.txt -w ~/work --clean-workdir
% -d smriprep=testData4bids_apps/output/smriprep --skip_bids_validation'
% disp('=> Starting smriprep...')
% mycmd= ['singularity run --cleanenv ', parentDir_containers,filesep,'fmriprep.simg ', path2bidsNewProj, ' ',...
%     path2bidsNewProj,filesep,'derivatives participant --participant-label ',subID, ...
%     ' --fs-license-file ', parentDir_containers,filesep, 'fs-license.txt -w ',path2bids,'/workdir -t ',taskID,...
%     ' --clean-workdir --anat-only --skip-bids-validation --fs-no-reconall'];
% system(mycmd);
% disp('-> smriprep is done...')

if bSkipFmriprep
    disp('=> Skipping fmriprep...')
else  
    disp('=> Starting fmriprep...')
    mycmd= ['singularity run --cleanenv ', parentDir_containers,filesep,'fmriprep-24.1.1.simg ', path2bidsNewProj, ' ',...
        path2bidsNewProj,filesep,'derivatives participant --participant-label ',subID, ...
        ' --fs-license-file ', parentDir_containers,filesep, 'fs-license.txt -w ',path2bids,'/workdir_fmriprep -t ',taskID,...
        ' --ignore slicetiming --cifti-output --skip-bids-validation --stop-on-first-crash ',...
        ' -d smriprep=',path2bidsNewProj,filesep,'derivatives',filesep,'smriprep'];
    %   ' --fs-subjects-dir ',path2bidsNewProj,filesep,'derivatives',filesep,'sourcedata',filesep,'freesurfer']; 
    % 
    % add '--cifti-output' so that xcp_d can be used handily.
    % add '--use-syn-sdc' for fieldmap-less susceptibility distortion
    % correction. However this one should not be intended for oblique images. 
    % use of '--force-syn' in addtion to field map sdc appear to result in distortion
    % in preprocessed fMRI timeseries. 
    system(mycmd);
    disp('-> fmriprep is done...')
end

%% step 3: postprocessing with xcp_d
% 'apptainer run --cleanenv ~/my_containers/xcp_d.simg
% amri_epi/derivatives/ amri_epi/derivatives/xcp_d participant --participant-label 052 --mode linc'
%% create the .json file of bids filter for xcp_d
jsonName = fullfile(path2bidsNewProj,'derivatives','xcp_d-bidsFilterFile.json');
json_postfilter.lh_pial_surf = struct('desc','preproc');
json_postfilter.lh_wm_surf = struct('desc','preproc');
json_postfilter.rh_pial_surf = struct('desc','preproc');
json_postfilter.rh_wm_surf = struct('desc','preproc');
bids.util.mkdir(fileparts(jsonName)); 
bids.util.jsonencode(jsonName,json_postfilter);

%%
if bSkipXcp
    disp('=> Skipping xcp-d...')
else    
    disp('=> Starting xcp-d...')
     mycmd= ['apptainer run --cleanenv ', parentDir_containers,filesep,'xcp_d-0.10.1.simg ', path2bidsNewProj,filesep,'derivatives ',...
        path2bidsNewProj,filesep,'derivatives/xcp_d participant --participant-label ',subID,' --mode linc -t ',taskID,...
        ' -w ',path2bids,'/workdir_xcp-d ',...
        ' --warp-surfaces-native2std --bids-filter-file ',path2bidsNewProj,filesep, 'derivatives/xcp_d-bidsFilterFile.json'];
    % mycmd= ['apptainer run --cleanenv ', parentDir_containers,filesep,'xcp_d-0.10.1.simg ', path2bidsNewProj,filesep,'derivatives ',...
    %     path2bidsNewProj,filesep,'derivatives/xcp_d participant --participant-label ',subID,' --mode abcd -t ',taskID,...
    %     ' -w ',path2bids,'/workdir_xcp-d --combine-runs y --linc-qc --motion-filter-type lp --band-stop-min 6',...
    %     ' --warp-surfaces-native2std --bids-filter-file ',path2bidsNewProj,filesep, 'derivatives/xcp_d-bidsFilterFile.json'];
%%% % use '--atlases' to specify which atlases are used, or use '--skip-parcellation' to skip this step entirely
    % use '--nuisance-regressors none' to skip the denoising step completely
    % --abcc-qc: to create the DCAN executive summary or QC 
    % --create-matrices all
    % --motion-filter-type lp --band-stop-min 6 ---- if scrubbing is
    % applied, this should be enabled, especially for infants and toddler
    % cohorts
    system(mycmd);
    disp('-> xcp-d is done...')
end

