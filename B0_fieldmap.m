function B0_fieldmap(path2refScan,mid_offlineRecon,methodofFM,varargin)
% %
% Calculate B0 maps (mag & pha) from reference scan in MoCo or Siemens
% double-echo GRE and store files in specified path and name following 
% standard BIDS.
% 
% The function incorporates
% 1) functions from Jiaen Liu's matlab package 'MoCo':
% https://github.umn.edu/wulab/moco-recon
% 2) a Julia package 'MRIFieldmaps.jl' from Jeff Fessler:
% https://github.com/MagneticResonanceImaging/MRIFieldmaps.jl
% 3) functions from Xiaoping Wu's matlab package:
% https://github.umn.edu/wulab/training-database-creation
%
% -- Shuxian, 1/6/2025
    
    %
    originalDir = pwd;
    cleanup = onCleanup(@() cd(originalDir));
    
    [filepath,filename,~] = fileparts(path2refScan);
    cd(filepath)
    
    if ismember(methodofFM,{'useAE4B0map', 'useOE4B0map', 'useREB40map'})
        expr = 'MID00*(\d+)_';
        extrNum = regexp(filename,expr,'tokens');
        midRef = str2num(extrNum{1}{1});
        midRecon = str2num(mid_offlineRecon);
    
        %
        [imgRef,paraRef,imgRefUncomb] = recon_amri_epi(midRef);
        te_arr = 1e-3*paraRef.te_contr;

        imgRef_e1 = imgRef(:,:,:,1,1,1);
        [~,thrd]=mask1d(imgRef_e1(:),0.1,0.08);
        mask=imgRef_e1>thrd;
        dimsM = ndims(mask);
        if dimsM>3
            subs = [repmat({':'},1,3), repmat({1},1,dimsM-3)];
            mask3D = mask(subs{:});
        elseif dimsM==3
            mask3D = mask;
        else
            mask3D = mask;
            warning('******The mask is not a 3D martix!******')
        end

        % calc dB0 maps
        switch methodofFM
            case 'useAE4B0map' % all echoes            
                b0_data_e3 = squeeze(sum(imgRefUncomb(:,:,:,:,:,1:end).*conj(imgRefUncomb(:,:,:,:,:,1)),5));
                te_e3 = paraRef.te_contr(1:end)*1e-3;
                % b0fieldmap_op = b0_map(b0_data_e3,te_e3).*mask3D;
                b0fieldmap_op = b0_map(b0_data_e3,te_e3);
    
            case 'useOE4B0map' % odd echoes
                b0_data_e2 = squeeze(sum(imgRefUncomb(:,:,:,:,:,1:2:end).*conj(imgRefUncomb(:,:,:,:,:,1)),5));
                te_e2 = paraRef.te_contr(1:2:end)*1e-3;
                % b0fieldmap_op = b0_map(b0_data_e2,te_e2).*mask3D;
                b0fieldmap_op = b0_map(b0_data_e2,te_e2);
    
            case 'useREB40map' % regularized estimation
                % % Initial field map        
                % nPhaRef = size(imgRefUncomb); % [nx,ny,nz,nmeas,nch,necho]
                % nMeas =nPhaRef(4);            
                % initFM = zeros(nPhaRef([1,2,3,4]));
                % % Combine phase images
                % B0_tmp= sum(imgRefUncomb(:,:,:,:,:,3).*conj(imgRefUncomb(:,:,:,:,:,1)),5);
                % if nMeas == 1
                %    initFM(:,:,:) = unwrapper_3d_mask( ...
                %                     angle(squeeze(B0_tmp(:,:,:))),mask3D) ...
                %                     /2/pi/(te_arr(3)-te_arr(1));
                % else
                %     error('More measurements than 1');
                % end
                
                % b1S = load(['result/mid',num2str(midRecon),'.b1.mat']);
                % b1=b1S.b1;
                imgRefUncomb = squeeze(imgRefUncomb);
                save('refUncomb.mat','imgRefUncomb','paraRef','mask3D');
        
                path2RegularizedFM = '/home/naxos2-raid27/qu000109/opt/JeffFessler/MRIFieldmaps/b0map2.jl';
                cmd = ['julia  ',path2RegularizedFM,' ',filepath,filesep,'refUncomb.mat',' ',filepath];
                % system(cmd); 
                system(['unset LD_LIBRARY_PATH; ',cmd]);

                str = load([filepath,filesep,'b0fieldmap.mat']);
                b0fieldmap_op = str.b0fieldmap(:,:,:,end);
            
        end  
        
        % Save nifti
        
        path2b0field = ['mid',num2str(midRecon),...
                            '_b0fieldmap.nii.gz'];
        siem_to_nifti(path2b0field,b0fieldmap_op,paraRef,1,1);
    
        imgRef_op = squeeze(imgRef(:,:,:,1,1,1));
        path2mag = ['mid',num2str(midRecon),...
                            '_magnitude.nii.gz'];
        siem_to_nifti(path2mag,imgRef_op,paraRef,1,1);

        return;

    elseif ismember(methodofFM,'useSDE4B0map')
        if length(varargin)<2
            error('Missing path or file number to source nifti file of field maps!');
        else 
            path2FM = varargin{1};
            numFM = num2str(varargin{2});
            numFM_mag = num2str(varargin{2}-1);              
        end

        filenameFM_pha = fullfile(path2FM,['0',numFM,'*.nii.gz']);
        filenameFM_mag = fullfile(path2FM,['0',numFM_mag,'*e1.nii.gz']);
        filenameFM_json = fullfile(path2FM,['0',numFM,'*.json']);
        FM_pha = dir(filenameFM_pha);
        FM_mag = dir(filenameFM_mag);
        FM_json = dir(filenameFM_json);
        path2FM_phase = fullfile(FM_pha.folder,FM_pha.name);
        path2FM_mag = fullfile(FM_mag.folder,FM_mag.name);
        path2FM_json = fullfile(FM_json.folder,FM_json.name);

        % phaDiff = niftiread(path2FM_phase);
        % mag = niftiread(path2FM_mag);
        % phaNinfo = niftiinfo(path2FM_phase); 
        % phaJinfo = bids.util.jsondecode(path2FM_json); 
        nii = load_untouch_nii(path2FM_phase);
        % mag = load_untouch_nii(path2FM_mag);
        phaJinfo = bids.util.jsondecode(path2FM_json);
        
        %
        file_des_mag = fullfile(filepath,['mid',mid_offlineRecon,'_magnitude.nii.gz']);
        copyfile(path2FM_mag, file_des_mag);


        % [~,magthrd]=mask1d(mag(:),0.1,0.03);
        % magmask=mag>magthrd;
        dte = (phaJinfo.EchoTime2 - phaJinfo.EchoTime1)*1e3;
        [~, freqmap] = calculate_b0map(nii.img, dte);        

        % save nifti
        % path2b0field = fullfile(filepath,['mid',mid_offlineRecon,...
        %                     '_b0fieldmap']);
        % phaNinfo.Datatype = class(b0fieldmap_op);
        % niftiwrite(b0fieldmap_op,path2b0field,phaNinfo);
        % gzip([path2b0field,'.nii']);
        % delete([path2b0field,'.nii']);
        nii2save = nii;
        nii2save.img = single(freqmap);
        
        nii2save.hdr.dime.datatype= 16;
        nii2save.hdr.dime.bitpix= 32;
        nii2save.hdr.dime.scl_slope= 1;
        nii2save.hdr.dime.scl_inter= 0;
        
        nii2save.hdr.dime.glmax= 4000;
        nii2save.hdr.dime.glmin= -4000;
        
        nii2save.img= single(freqmap);
        path2b0field = fullfile(filepath,['mid',mid_offlineRecon,...
                             '_b0fieldmap.nii.gz']);
        save_untouch_nii(nii2save,path2b0field)
    

        % path2mag = fullfile(filepath,['mid',mid_offlineRecon,...
        %                     '_magnitude']);
        % phaNinfo.Datatype = class(mag);
        % niftiwrite(mag,path2mag,phaNinfo);
        % gzip([path2mag,'.nii']);
        % delete([path2mag,'.nii']);

        return;

    end
        
   
    
    
    

