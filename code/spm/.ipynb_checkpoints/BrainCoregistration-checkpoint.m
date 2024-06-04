% CL Septembre 2022, landelle.caroline@gmail.com // caroline.landelle@mcgill.ca

% Toolbox required: Matlab, SPM12


%% 
function Norm=Norm_BMPD(filename_anat, filename_funcmean,inputDir,filenames_func4D)

%______________________________________________________________________
%% Initialization 
%______________________________________________________________________
SPM_Dir='/cerebro/cerebro1/dataset/bmpd/derivatives/thibault_test/code/toolbox/spm12/'
addpath(SPM_Dir); % Add SPM12 to the path
if ~isempty(filenames_func4D);
f = spm_select('ExtFPList',fullfile(inputDir),filenames_func4D,Inf); % listes d'images 3D 
elseif filenames_func4D == "";
f='';
end;
f
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = filename_anat
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {filename_funcmean}
matlabbatch{1}.spm.spatial.coreg.estwrite.other = cellstr(f);
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

Norm='Coregistration Done'
clear matlabbatch
end   

