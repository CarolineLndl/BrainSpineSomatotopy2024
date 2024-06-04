# -*- coding: utf-8 -*-
import os, glob, shutil
import json
import pandas as pd
import numpy as np
import nibabel as nib

import matlab.engine

import brsc_utils

#plotting:
import matplotlib
import matplotlib.pyplot as plt

from nilearn import image
   
class Preprocess_BrSc:
    '''
    The Preprocess class is used to compute Brain and spinal cord preprocessings simultaneously
    Slice timing & image cropping are availables
    Attributes
    ----------
    config : dict
    '''
    
    def __init__(self, config, verbose=True):
        self.config = config # load config info
        self.participant_IDs= self.config["participants_IDs"] # list of the participants to analyze
        self.main_dir=self.config["main_dir"] # main drectory of the project
        self.raw_dir=self.config["main_dir"]+self.config["raw_dir"] # directory of the raw data
        self.preprocess_dir=self.config["main_dir"]+config["preprocess_dir"]["main_dir"] # directory for the preprocessings
        
        if verbose==True:    
            print("Lines 0-10 of config files should be manually modified first")
            print("All the raw data should be store in BIDS format")
            print(" ")
             
        # Create participant directories (if not already existed)
        for ID in self.participant_IDs:
            ID_preproc_dir=self.preprocess_dir + "sub-" + ID
            if not os.path.exists(ID_preproc_dir):
                os.mkdir(ID_preproc_dir)
                
                # create 1 folder per session if there are multiple sessions (exemple multiple days of acquisition)
                for ses_name in self.config["design_exp"]['ses_names']:
                    ses_dir="/ses-" + ses_name if int(self.config["design_exp"]["ses_nb"])>1 else ""
                    if ses_dir != "":
                        os.mkdir(ID_preproc_dir +  ses_dir)
                    os.mkdir(ID_preproc_dir +  ses_dir + "/anat/")
                    os.mkdir(ID_preproc_dir +  ses_dir + "/anat/brain/")
                    os.mkdir(ID_preproc_dir + ses_dir + "/func/")
                         
                print("New folders in preprocess dir have been created")
                
            # if there are multiple runs create a folder for each runs in func folder:
            for ses_name in self.config["design_exp"]['ses_names']:
                ses_dir="ses-" + ses_name if int(self.config["design_exp"]["ses_nb"])>1 else ""
                for task_name in self.config["design_exp"]['task_names']:
                    task_dir='task-'+task_name if int(self.config["design_exp"]["task_nb"])>1 else ""
                    if not os.path.exists(ID_preproc_dir +"/" +ses_dir +"/func/" + task_dir):
                        os.mkdir(ID_preproc_dir +"/" +ses_dir +"/func/" + task_dir)
            
            # Check if there is specificity for anat or func filename:
            if ID in config['files_specificities']["anat"]:
                print("sub-" + ID + " have a anat filename specitity: " + config['files_specificities']["anat"][ID])
            
            if ID in config['files_specificities']["func"]:
                print("sub-" + ID + " have a anat filename specitity: " + config['files_specificities']["func"][ID])
         
            
    
    def stc(self,ID=None,i_img=None,json_f=None,ses_name='',task_name='',tag='',t_custom=False,down=None,odd=None, redo=False,verbose=True):
        '''
        Slice timing correction is applied in order to minimize the effect of slice ordering in the acquisition of the images.  
        https://poc.vl-e.nl/distribution/manual/fsl-3.2/slicetimer/index.html
        
        Attributes:
        ----------
        ID: name of the participant
        i_img: input filename of functional images (str, default:None, the file will be defined by default in the function)
        json_f: json file with relevant info (str, default:None, the file will be defined by default in the function)
        ses_name: if the ses have any specific name (ses- in BIDS format)
        task_name: if the run have any specific name (task- in BIDS format)
        tag: if there is any specific tag related to the filename. eg. run-02
        tr: Specify TR of data - /!\ (str, default:3)
        down:  reverse slice indexing (If slices were acquired from the top of the SC to the bottom)  
        odd:  use it for interleaved acquisition (str, default:)
        tcustom: set True to use this option > filename of single-column custom interleave timing file. The units should be in TRs, with 0.5 corresponding to no shift. Therefore a sensible range of values will be between 0 and 1. option parameters provide in the json file "SliceTiming" can be used. The array is in seconds, one number per slice is provide and the maximum value is always be less that the TR.
This array should be transform in TRs units (value/TR)  before using it with slicetimer it was not the case in this script but it should not change the results for resting state correlation analyses </font>


        Outputs: 
        ----------
        Slice time corrected image *stc.nii
        slice timings file : *stc.txt  
        nb: One value (ie for each slice) on each line of a text file. The units are in TRs, with 0.5 corresponding to no shift. 
        '''
        
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")

        # Select the default input filename if it is not provided
        if i_img==None:
            #print(self.raw_dir + "sub-" + ID+ "/" + ses_name + "/func/sub-" + ID + "*" + task_name +"*" +'*' + tag + "*.nii*")
            i_img=glob.glob(self.raw_dir + "sub-" + ID+ "/" + ses_name + "/func/sub-" + ID + "*" + task_name +"*" +'*' + tag + "*.nii*")[0]
        if json_f==None:
            json_f=glob.glob(i_img.split(".")[0] + "*.json")[0]
        #print("input image is: " + i_img)
        
        # Create output directories if no existed
        ID_dir=self.preprocess_dir + "/sub-" + ID  
        stc_dir=ID_dir + '/'+ses_name+"/func/"+task_name+ '/' +self.config["preprocess_dir"]["func_stc"]

        if not os.path.exists(stc_dir):
            os.mkdir(stc_dir)
            os.mkdir(stc_dir + "/brain")
            os.mkdir(stc_dir + "/spinalcord")
        
        # Define output name:
        o_img=stc_dir + os.path.basename(i_img.split(".")[0] + "_stc.nii.gz")

        if not os.path.exists(o_img) or redo==True:
            print(">>>>> slice timing correction is running for sub-" + ID)
            
            # read the json file to extract some info
            with open(json_f) as g:
                params = json.load(g)
            tr=params["RepetitionTime"] # extract the time repetition value
            if t_custom ==True:
                o_txt=stc_dir + os.path.basename(i_img.split(".")[0] + "_stc.txt") # output with slicetiming info

                stc_info=params["SliceTiming"] # provide info about interleave slice order

                with open(o_txt, 'w') as f:
                    for item in stc_info:
                        item_tr=item/tr # Transform slicetiming in secs in st in TRs units
                        f.write('{}\n'.format(item_tr))
                f.close()
                del f, item

            # run slice timing correction:
                string="slicetimer -i "+ i_img + " -o " + o_img +" -r " + str(tr) + " --odd --tcustom=" + o_txt
                os.system(string)
                print("done")
                
            else:
                raise Warning("No other option that interleaved acquisition for sct was implemented yet, you should cutomise the code here to add an option")
        
        if os.path.exists(stc_dir) and redo==False:
            print(">>>>> slice timing correction was already completed for sub-" + ID)
        
        return o_img

    def crop_img(self,ID=None,i_img=None,o_folder=None,tsv_f=None,structure=None,ses_name='',task_name='',tag='',img_type="func",redo=False,verbose=True):
        
        '''
        This function will help to crop brain and spinale cord from an image included both in one single FOV
        to improve: the option to use automatic detection of the brain and the spinal cord
        
        Attributes:
        ----------
        ID: name of the participant
        i_img: input filename of functional images (str, default:None, an error will be raise)
        o_img: output folder name filename (str, default:None, the input folder will be used)
        
        tsv_f: tsv file with relevant info (str, default:None, the file will be defined by default in the function)
            should include at least the following columns: participant_id	anat_crop_brain	anat_crop_spine	func_crop_brain	func_crop_spine
            
        ses_name: if the session have any specific name (ses- in BIDS format)
        task_name: if the run have any specific name (task- in BIDS format)
        img_type: the type of input image should be specify "func" or "anat"

        Outputs: 
        ----------
        Image cropped
        '''
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if i_img==None:
            raise Warning("Please provide filename of the input file")
      
        # Select the default output directory (input directory) 
        if o_folder==None:
            o_folder=os.path.dirname(i_img)
            
        #read tsv file
        if tsv_f==None:
            tsv_f=glob.glob(self.preprocess_dir + "*participants.tsv")[0]
        df = pd.read_csv(tsv_f, sep='\t') # read the file
        
        # this lines should be remove for next version
        if img_type=="func" and structure==None:
            o_img_brain=o_folder + "/brain/" + os.path.basename(i_img.split(".")[0] + "_brain.nii.gz")
            o_img_sc=o_folder + "/spinalcord/" + os.path.basename(i_img.split(".")[0] + "_sc.nii.gz")
        else :
            o_img_brain=o_folder  + os.path.basename(i_img.split(".")[0] + "_brain.nii.gz")
            o_img_sc=o_folder  + os.path.basename(i_img.split(".")[0] + "_sc.nii.gz")
        
        # crop brain:
        if not os.path.exists(o_img_sc) or redo==True:
            if structure==None or structure=="brain":
                print(">>>>> brain and spinal cord cropping is running for sub-" + ID + " " + img_type)
                z_max=int(df[img_type+"_zmax_brain"][df["participant_id"]==ID].values[0]) if img_type+"_zmax_brain" in df.columns else nib.load(i_img).shape[2]
                z_min=int(df[img_type+"_zmin_brain"][df["participant_id"]==ID].values[0]) if img_type+"_zmin_brain" in df.columns else 0

                string_br='sct_crop_image -i ' + i_img+ ' -o '+ o_img_brain+' -zmin ' + str(z_min) + ' -zmax ' + str(z_max)
                os.system(string_br) # run the string as a command line
            
            if structure==None or structure=="spinalcord":
                # crop spinalcord:
                z_max=int(df[img_type+"_zmax_sc"][df["participant_id"]==ID].values[0]) if img_type+"_zmax_sc" in df.columns else nib.load(i_img).shape[2]
                z_min=int(df[img_type+"_zmin_sc"][df["participant_id"]==ID].values[0]) if img_type+"_zmin_sc" in df.columns else 0
                string_sc='sct_crop_image -i ' + i_img+ ' -o '+ o_img_sc+' -zmin '+str(z_min)+' -zmax ' + str(z_max)
                os.system(string_sc) # run the string as a command line

        return o_img_brain, o_img_sc
    
    def smooth_img(self,i_img=None,o_folder=None,fwhm=[0,0,0],ses_name='',task_name='',tag='_s.nii.gz',redo=False,verbose=True):            
        '''
        This function will smooth the 3D input image using nilearn see:
        https://nilearn.github.io/stable/modules/generated/nilearn.image.smooth_img.html#nilearn.image.smooth_img
          
        Attributes:
        ----------
        i_img <filename>, mendatory, default=None : input filename of functional images (str, default:None, an error will be raise)
        o_img <filename>, optional, default=None : output folder name filename (str, default:None, the input folder will be used)
        fwhm <array>, optional, default=[0,0,0]: it must have 3 elements, giving the FWHM along each axis. If any of the elements is 0 or None, smoothing is not performed along that axis.
                
        ses_name: if the run have any specific name (ses- in BIDS format)
        task_name: if the run have any specific name (task- in BIDS format)
        img_type: the type of input image should be specify "func" or "anat"

        Outputs: 
        ----------
        Image cropped
        '''
        
        if i_img==None:
            raise Warning("Please provide filename of the input file")
        
        # Default ouput folder will be i_img folder:
        if o_folder==None:
            o_folder=os.path.dirname(i_img)
            
        o_img=o_folder + "/"+ os.path.basename(i_img).split('.')[0] + tag
        if not os.path.exists(o_img) or redo==True:
            smoothed_image=image.smooth_img(i_img, fwhm)
            smoothed_image.to_filename(o_img)
            string='fslmaths '+o_img+' -Tmean '+o_img.split('.')[0] + '_mean.nii.gz'
            os.system(string)
            # save smoothing parameters                    
            with open(o_img.split('.')[0] + '.json', 'w') as f:
                json.dump(fwhm, f) # save info

            print("Smoothing done: " + os.path.basename(o_img))   
        
        return o_img
########################################################################            
########################################################################
class Preprocess_Br:
    '''
    The Preprocess_Br class is used to compute brain preprocessings 
    Motion correction, Segmentation, vertebral labelling and coregistration steps are available
    
    Attributes
    ----------
    config : dict
    '''
    def __init__(self, config, verbose=True):
        self.config = config # load config info
        self.participant_IDs= self.config["participants_IDs"] # list of the participants to analyze
        self.main_dir=self.config["main_dir"] # main drectory of the project
        self.raw_dir=self.config["main_dir"]+self.config["raw_dir"] # directory of the raw data
        self.preprocess_dir=self.config["main_dir"]+config["preprocess_dir"]["main_dir"] # directory for the preprocessings
        
    
            
    def moco(self,ID=None,i_img=None,o_folder=None,params=None, ses_name='',task_name='',tag='',plot_show=True,redo=False,verbose=True):
        
        '''
        This function will correct the brain images for motion using the fsl
        Calculate framewise displacement (FD)
        Plot the motion parameters
        https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/MCFLIRT
         
        Attributes:
        ----------
        ID: name of the participant
        i_img: input filename of functional images 4D (str, default:None, an error will be raise)
        o_folder: output folder (str, default:None, the input folder will be used)
        params: see sct toolbox for more details, should not be modified across participants
        ses_name: if the run have any specific name (ses- in BIDS format)
        task_name: if the run have any specific name (task- in BIDS format)
        
        Outputs: 
        ----------
        The outputs of the motion correction process are:
            - the motion-corrected fMRI volumes: *_sc_moco.nii & one temporary file for each volume
            - the time average of the corrected fMRI volumes: - Time average of corrected vols: *_sc_moco_mean.nii
            - a time-series with 1 voxel in the XY plane, for the X and Y motion direction (two separate files).
            - a TSV file with the slice-wise average of the motion correction for XY (one file), that can be used for Quality Control.

- Motion corrected volumes: 

- Time-series with 1 voxel in the XY plane: *_sc_moco_params_X.nii & params_Y.nii
- Slice-wise average of MOCO for XY : moco_params.tsv & moco_params.txt
        '''
            
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if i_img==None:
            raise Warning("Please provide filename of the input file")
       
        #1. create ouput folder
        if o_folder is None : # gave the default folder name if not provided
            o_folder=self.preprocess_dir + "sub-" + ID+ "/" + ses_name + "/func/"+task_name+ "/" +self.config["preprocess_dir"]["func_moco"]

        if not os.path.exists(o_folder):
            os.mkdir(o_folder)
        if not os.path.exists(o_folder+ "/brain/"):
            os.mkdir(o_folder + "/brain/")
            
        mask_img= self.preprocess_dir + "sub-" + ID+ "/" + ses_name + "/func/"+task_name+ "/" +self.config["preprocess_dir"]["func_brain_seg"]+ os.path.basename(i_img).split(".")[0] + "_masked.nii.gz"
        moco_file= o_folder+"/brain/"+ os.path.basename(i_img).split(".")[0] + "_moco.nii.gz"
        moco_mean_file= o_folder+"/brain/"+ os.path.basename(i_img).split(".")[0] + "_moco_mean.nii.gz"

        if not os.path.exists(moco_file) or redo==True:
            print(">>>>> Moco is running for the brain image of the sub-" + ID + " ")
       
            # 1. Create a binary mask before apply moco
            string_mask="bet "+i_img+ " "+mask_img+" -F" 
            os.system(string_mask)
            
            # 2. compute motion correction
            string_moco="mcflirt -in "+ mask_img+" -out "+moco_file+" -mats -plots"
            os.system(string_moco)
            
            #3. Calculate the mean image
            string_mean="fslmaths " +moco_file+ " -Tmean " + moco_mean_file
            os.system(string_mean)
            
        if plot_show==True:
            #4. Calculate Framewise displacement
            output_fd=o_folder + '/brain/' + os.path.basename(i_img).split('.')[0] + "_FD_brain"
            if not os.path.exists(output_fd + ".txt") or redo==True:
                print('Compute framewise displacement for sub-' +ID)
                string_fd="fsl_motion_outliers -i "+mask_img+" -o "+o_folder+"/brain/"+" -s "+ output_fd + ".txt -p "+ output_fd + ".png --fd"
                os.system(string_fd)

            FD=pd.read_csv(output_fd + '.txt', delimiter=',',header=None)
            FD_mean = np.mean(FD[0])
            print('Mean FD for sub-' +ID + ": " + str(round(FD_mean,3)) + ' mm')
            
            
            #5. Plot motion parameters
            fig, axs = plt.subplots(2,1, figsize=(18, 6), facecolor='w', edgecolor='k')
            fig.tight_layout()
            fig.subplots_adjust(hspace = .5, wspace=.001)
            moco_params=pd.read_csv(moco_file + '.par',delimiter='  ',header=None,engine='python')
    
            axs[0].plot(moco_params.iloc[:,0:3]) 
            axs[0].set_title(ID + " " + ses_name)
            axs[1].plot(moco_params.iloc[:,3:6])
            
            
            axs[0].set_ylabel("Rotation (rad)")
            axs[1].set_ylabel("Translation (mm)")
            axs[1].set_xlabel("Volumes")
            
            if not os.path.exists(o_folder + "/brain/moco_params.png") or redo==True:
                plt.savefig(o_folder + "/brain/moco_params.png")
                np.savetxt(o_folder + '/brain/FD_mean.txt', [FD_mean])
            plt.show(block=False)
        print(" ")
        return moco_file, moco_mean_file
    
    
    def segmentation(self,ID=None,i_img=None,o_folder=None,ses_name='',task_name='',tag='',img_type="anat",redo=False,verbose=True):
        
        '''
        This function will segment the brain using CAT12 toolbox:

        https://neuro-jena.github.io/cat/
        
        Attributes:
        ----------
        ID: name of the participant
        i_img: input filename of functional images 3D (str, default:None, an error will be raise)
        
        o_folder: output folder (str, default:None, the input folder will be used)
        ses_name: if the session have any specific name (ses- in BIDS format)
        task_name: if the run have any specific name (task- in BIDS format)
        img_type: the type of input image should be specify "func" or "anat"
        
        Outputs: 
        ----------
        
        '''
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if i_img==None:
            raise Warning("Please provide filename of the input file")
        
           
        #1. create ouput folder
        if o_folder is None : # gave the default folder name if not provided
            if img_type=="func":
                o_folder=self.preprocess_dir + "sub-" + ID+ "/" + ses_name + "/func/"+task_name+ "/" +self.config["preprocess_dir"]["func_brain_seg"]
            elif img_type=="anat":
                o_folder=self.preprocess_dir + "sub-" + ID+ "/anat/"+self.config["preprocess_dir"]["anat_brain_seg"]

        if not os.path.exists(o_folder):
            os.makedirs(o_folder)
        
        #2. unzip files
        if os.path.basename(i_img).split(".")[-1] == "gz" :
            unzip_i_img=brsc_utils.unzip_file(i_img,o_folder=o_folder,ext=".nii",zip_file=False, redo=False,verbose=False)
        else:
            unzip_i_img=i_img
            
        ##3. run segmentation with cat12
        if not os.path.exists(o_folder + "/mri/") or redo==True:
            print('Brain segmentation is running for sub-' +ID)
            os.chdir(self.config["tools_dir"]["main_codes"] + '/code/spm/')# need to change the directory to find the CAT12 function
            eng = matlab.engine.start_matlab()
            eng.CAT12_BrainSeg(ID,unzip_i_img, self.config["tools_dir"]["main_codes"] + '/code/spm/Cat12_log/',self.config["tools_dir"]["spm_dir"])
            
            #Mask the brain image using segmented tissues
            mask=o_folder + "/mri/p0*" + os.path.basename(unzip_i_img)
            string_mask='fslmaths '+unzip_i_img +" -mas " + mask + " >0 " +unzip_i_img.split(".")[0] + "_masked.nii.gz"
            os.system(string_mask)
            
            
        elif verbose==True:
            print('Brain segmentation alredy exists for sub-' +ID + " here: " + o_folder)
        
        return o_folder + "/mri/"
        
    def coreg_anat2func(self,ID=None,func_img=None,anat_img=None, other_img=None,wm_file=None,tag="_",
                        ses_name="",task_name="",o_folder=None,redo=False,verbose=True):
        '''
        to do: find a way to do this step with SPM (or the next one with fsl)
        
        Attributes:
        ----------
        ID: name of the participant
        source_img <filename>: functional image filename 3D  (str, default:None, an error will be raise)
        anat_img <filename>: anat image filename 3D  (str, default:None, an error will be raise)
        other_img <filename>: white other file  (str, default:None, optional)
        o_folder <folder dir>: output folder (str, default:None, the input folder will be used)
        
        '''
        
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if anat_img==None or func_img==None or wm_file==None:
            raise Warning("Please provide filename of the source_img and ref_img and wm_file")
        
        if o_folder==None:
            o_folder=self.preprocess_dir + "sub-" + ID+ "/" + ses_name + "/func/"+task_name+ "/" +self.config["preprocess_dir"]["func_coreg"] + "/brain/"
        
        if not os.path.exists(self.preprocess_dir + "sub-" + ID+ "/" + ses_name + "/func/"+task_name+ "/" +self.config["preprocess_dir"]["func_coreg"]):
            os.mkdir(self.preprocess_dir + "sub-" + ID+ "/" + ses_name + "/func/"+task_name+ "/" +self.config["preprocess_dir"]["func_coreg"])
        if not os.path.exists(o_folder):
            os.mkdir(o_folder)# create output folder if not exists
    

        o_filename=o_folder + os.path.basename(func_img).split(".")[0] + tag + ".nii.gz"
        warp_field=o_filename.split(".")[0] +".mat"
        print(warp_field)
        if not os.path.exists(o_filename) or redo==True:
            string_flirt="flirt -in "+func_img+" -ref "+anat_img+" -o "+o_filename+" -omat "+warp_field+" -dof 6 -cost bbr -wmseg " + wm_file 
            os.system(string_flirt)
            
        
            if other_img:
                #Create inverse warping field:
                inv_warp=warp_field.split('.mat')[0] + '_inv.mat' # inverse transformation matrice (from t1w to func)
                string_2="convert_xfm -omat "+inv_warp+" -inverse "+warp_field
                os.system(string_2)
    
    
                # b. Apply transformation other file in ref to source
                anat_folder=self.preprocess_dir + "sub-" + ID+ "/anat/"+self.config["preprocess_dir"]["anat_brain_seg"] + "/mri/"
                o_other_img=anat_folder + os.path.basename(other_img).split(".")[0] + "_inFunc.nii.gz"
                string_3="flirt -in "+other_img+" -ref " + func_img  + " -o " + o_other_img +" -omat " +inv_warp+" -dof 6 -interp spline"
                os.system(string_3)
                
                # c. Transform the output image in a binary image 
                string_4="fslmaths "+o_other_img+" -thr 0.5 -bin " +o_other_img
                os.system(string_4)
                
        
    def coreg_func2anat(self,ID=None,ref_img=None,source_img=None,other_files=None,tag="",o_folder=None,redo=False,verbose=True):
        '''
        ID: name of the participant
        source_img <filename>: functional image filename 3D  (str, default:None, an error will be raise)
        ref_img <filename>: anat image filename 3D  (str, default:None, an error will be raise)
        other_files <filename>: functional image filename 4D  (str, default:None, optional)
        o_folder <folder dir>: output folder (str, default:None, the input folder will be used)
        dest_img <filename>: destination filename 
        warping_fields <filename>: Transformation(s), which can be warping fields (nifti image) or affine transformation matrix (text file). Separate with space.

        '''
        
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if source_img==None or ref_img==None:
            raise Warning("Please provide filename of the source_img and ref_img")
    
    
        # Unizip anat and mean func data
        if os.path.basename(ref_img).split(".")[-1] == "gz" :
            unzip_ref_img=brsc_utils.unzip_file(ref_img,o_folder="/"+os.path.dirname(ref_img)+"/",ext="_SPM.nii",redo=True,verbose=False)
        else:
            unzip_ref_img=ref_img.split()

        if os.path.basename(source_img).split(".")[-1] == "gz" :
            unzip_source_img=brsc_utils.unzip_file(source_img,o_folder=os.path.dirname(source_img) +"/",ext=tag+".nii",redo=True,verbose=False)
        else:
            unzip_source_img=source_img.split('.')[0]  + tag + ".nii"
            shutil.copy(source_img, unzip_source_img)
        
            
       # Unzip 4d data to apply warpping field
        if other_files:
            if os.path.basename(other_files).split(".")[-1] == "gz":
                unzip_other_files=brsc_utils.unzip_file(other_files,o_folder=os.path.dirname(other_files),ext=tag+".nii",redo=True,verbose=False)
            else:
                unzip_other_files=other_files.split('.')[0] + tag + ".nii"
                shutil.copy(other_files, unzip_other_files)
            unzip_other_filenames=os.path.basename(unzip_other_files)
            unzip_other_files_dir=os.path.dirname(unzip_other_files)
            
        elif other_files==None:
            unzip_other_files=other_files
            unzip_other_filenames=""
            unzip_other_files_dir=""
            
        # output info:
        if o_folder==None:
            o_folder=self.preprocess_dir + "sub-" + ID+ "/" + ses_name + "/func/"+task_name+ "/" +self.config["preprocess_dir"]["func_coreg"] + "/brain/"
        if not os.path.exists(o_folder):
            os.mkdir(o_folder)# create output folder if not exists
        o_filename= o_folder+ "/"+ os.path.basename(unzip_source_img).split('.')[0] +".nii" # output filename
        o_filename_r= o_folder+ "/"+ os.path.basename(unzip_source_img).split('.')[0] +"_r.nii" # output filename
        
        # run coregistration into anat:
        if not os.path.exists(o_filename) or redo==True:
            os.chdir(self.config["tools_dir"]["main_codes"] +'/code/spm/')# need to change the directory to find the BrainCoregistration function
            eng = matlab.engine.start_matlab()
            eng.BrainCoregistration(unzip_ref_img,unzip_source_img,unzip_other_files_dir,unzip_other_filenames)
           
            print('Coregistration into anat space done')
            os.rename(os.path.dirname(unzip_source_img) + "/"+ os.path.basename(unzip_source_img),o_filename) # rename the output
            
            os.rename(os.path.dirname(unzip_source_img) + "/r"+ os.path.basename(unzip_source_img),o_filename_r) # rename the output
            
            if other_files == None:
                os.remove(o_folder+ "/"+ os.path.basename(unzip_source_img).split('.')[0] +tag+".nii") # a default file will be created but we don't want it
                
            if other_files != None:
                o_other_filename= o_folder+ "/"+ os.path.basename(unzip_other_files).split('.')[0] +tag+".nii" # output filename
                o_other_filename_r= o_folder+ "/"+ os.path.basename(unzip_other_files).split('.')[0] +tag+"_r.nii" # output filename
       
                os.rename(os.path.dirname(unzip_other_files) + "/"+ os.path.basename(unzip_other_files),o_other_filename) # rename the output
                os.rename(os.path.dirname(unzip_other_files) + "/r"+ os.path.basename(unzip_other_files),o_other_filename_r) # rename the output
            
            
        else:
            print('Brain functional coregistration into anat space was already done, set redo=True to run again the coregistration')
        
        return o_filename
        
    
                                                   
        if threshold:                                                    
            # b. Transform the output image in a binary image
            string2="fslmaths "+o_img+" -thr "+threshold+" -bin " + o_img
            #os.system(string2)
                
        else:
            if verbose:
                print("Tranformation was already applied put redo=True to redo that step")
            
        return o_img

    
    def normalisation(self,ID=None,warp_file=None,coreg2anat_file=None,o_file=None,brain_mask=None,redo=False):
        '''
            ID: name of the participant
            warp_file <filename>: warping field from anat to MNI space
            coreg2anat_file <filename>: functional image coregister into anat space
            o_file <filename>: output filename
            brain_mask <filename>: apply a mask

        '''
        if not os.path.exists(o_file) or redo==True:
            os.chdir(self.config["tools_dir"]["main_codes"] +'/code/spm/') # need to change the directory to find the function
            eng = matlab.engine.start_matlab()
            print(eng.BrainNormalisation(warp_file,coreg2anat_file)) #for quality check
            o_files_def=os.path.split(coreg2anat_file)[0]+'/w'+coreg2anat_file.split('/')[-1] # default output
    
            os.rename(o_files_def,o_file) #move and rename the file
            string="fslmaths "+o_file+ " -nan " +o_file; os.system(string)# remove nan values
            os.remove(o_file) # remove .nii file and keep only .nii.gz

            if brain_mask is not None:
                string='fslmaths '+o_file+'.gz -mas '+brain_mask+' '+o_file
                os.system(string)
                #os.remove(func_norm_file)
            
            print("output: " + os.path.basename(o_file))

            return print('Normalisation into MNI space done')
        else:
            return print('Normalisation into MNI space was already done, set redo=True to run again the coregistration')
        
########################################################################            
########################################################################
class Preprocess_Sc:
    '''
    The Preprocess class is used to compute spinal cord preprocessings 
    Motion correction, Segmentation and coregistration steps are available
    
    Attributes
    ----------
    config : dict
    '''
    def __init__(self, config, verbose=True):
        self.config = config # load config info
        self.participant_IDs= self.config["participants_IDs"] # list of the participants to analyze
        self.main_dir=self.config["main_dir"] # main drectory of the project
        self.raw_dir=self.config["main_dir"]+self.config["raw_dir"] # directory of the raw data
        self.preprocess_dir=self.config["main_dir"]+config["preprocess_dir"]["main_dir"] # directory for the preprocessings
        
    def moco_mask(self,ID=None,i_img=None,o_folder=None, radius_size=15,task_name='',ses_name='',tag='',manual=False,redo=False,verbose=True):
        
        '''
        This function will create mask arround a centerline
        https://spinalcordtoolbox.com/user_section/command-line.html#sct-get-centerline
        https://spinalcordtoolbox.com/user_section/command-line.html#sct-create-mask
        
        to do: rename name of the output folder by moco_mask
        Attributes:
        ----------
        ID: name of the participant
        i_img: input filename of functional images (str, default:None, an error will be raise)
        o_img: output folder name filename (str, default:None, the input folder will be used)
        radius_size: value of the diameter of the surrounding mask in voxels (default is 15)
        sesname: if the session have any specific name (ses- in BIDS format)
        task_name: if the run have any specific name (task- in BIDS format)
        
        Outputs: 
        ----------
        Image cropped
        '''
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if i_img==None:
            raise Warning("Please provide filename of the input file")
  
        #1. create ouput folder
        if o_folder is None : # gave the default folder name if not provided
            o_folder=self.preprocess_dir + "sub-" + ID+  "/" + ses_name + "/func/"+task_name+ "/" + self.config["preprocess_dir"]["func_mask"]

        if not os.path.exists(o_folder):
            os.mkdir(o_folder)
            os.mkdir(o_folder+ "/spinalcord/")
            
        #2 Create the centerline, if manual centerline is needed change: -method optic by -method viewer
        centerline_f=o_folder + "/spinalcord/" + os.path.basename(i_img).split(".")[0] + "_centerline"
        mask_f=o_folder + "/spinalcord/" + os.path.basename(i_img).split(".")[0] + "_mask.nii.gz"
        
        if not os.path.exists(mask_f) or redo==True:
            
            if manual==True:
                method="viewer"
            else:
                method="optic"

            string_cntr="sct_get_centerline -i "+ i_img +" -o "+centerline_f +" -c t1 -method "+ method+" -centerline-algo bspline"
            os.system(string_cntr)
        
            #3. Create the mask arround the centerline
            string_mask="sct_create_mask -i "+i_img +" -p centerline,"+centerline_f+".nii.gz -size "+ str(radius_size)+ " -o " + mask_f
            os.system(string_mask)
        
            if verbose ==True:
                print("Mask created, check the outputs files in fsleyes by copy and paste:")
                print("fsleyes " + mask_f)
                print("use manual=False if the centerline need manual corrections")

        return mask_f, centerline_f+'.nii.gz'
           
       
            
    def moco(self,ID=None,i_img=None,mask_img=None,o_folder=None,params=None,
             ses_name='',task_name='',tag='',abs_moco_param=False,plot_show=True,redo=False,verbose=True):
        
        '''
        This function will correct the spinal cord images for motion using the sct toolbox
        Calculate framewise displacement (FD)
        Plot the motion parameters
        https://spinalcordtoolbox.com/user_section/command-line.html#sct-fmri-moco
        
        to do: rename name of the output folder by moco_mask
        recalculate volume wise motion parameters : erros in the SCT
        
        Attributes:
        ----------
        ID: name of the participant
        i_img: input filename of functional images 4D (str, default:None, an error will be raise)
        mask_img: Binary mask to limit voxels considered by the registration metric (str, default:None, an error will be raise)
        o_folder: output folder (str, default:None, the input folder will be used)
        params: see sct toolbox for more details, should not be modified across participants
        ses_name: if the session have any specific name (ses- in BIDS format)
        task_name: if the run have any specific name (task- in BIDS format)
        abs_moco_param <Bool> optional, defautl=False: the default moco params were wrongly calculated by the sct toolbox, we propose to calculate the absolute mean value instead
        Outputs: 
        ----------
        The outputs of the motion correction process are:
            - the motion-corrected fMRI volumes: *_sc_moco.nii & one temporary file for each volume
            - the time average of the corrected fMRI volumes: - Time average of corrected vols: *_sc_moco_mean.nii
            - a time-series with 1 voxel in the XY plane, for the X and Y motion direction (two separate files).
            - a TSV file with the slice-wise average of the motion correction for XY (one file), that can be used for Quality Control.

- Motion corrected volumes: 

- Time-series with 1 voxel in the XY plane: *_sc_moco_params_X.nii & params_Y.nii
- Slice-wise average of MOCO for XY : moco_params.tsv & moco_params.txt
        '''
        
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if i_img==None:
            raise Warning("Please provide filename of the input file")
        
        if mask_img==None:
            raise Warning("Please provide filename of the mask file")
  
  
        #1. create ouput folder
        if o_folder is None : # gave the default folder name if not provided
            o_folder=self.preprocess_dir + "sub-" + ID+ "/" + ses_name + "/func/"+task_name+ "/" +self.config["preprocess_dir"]["func_moco"]

        if not os.path.exists(o_folder):
            os.mkdir(o_folder)
        if not os.path.exists(o_folder+ "/spinalcord/"):
            os.mkdir(o_folder + "/spinalcord/")
            
        # Define moco paramaters (should be the same for the entire dataset)
        if params==None:
            params = 'poly=0,smooth=1,metric=MeanSquares,gradStep=1,sampling=0.2'#,numTarget=' + str(num_target)
        
        
        moco_file= o_folder+"/spinalcord/"+ os.path.basename(i_img).split(".")[0] + "_moco.nii.gz"
        moco_mean_file= o_folder+"/spinalcord/"+ os.path.basename(i_img).split(".")[0] + "_moco_mean.nii.gz"
        print(moco_file)
        if not os.path.exists(moco_file) or redo==True:
            #print(mask_img)
            print(">>>>> Moco is running for the spinal cord image of the sub-" + ID + " ")
            string="sct_fmri_moco -i "+i_img+" -m "+mask_img+" -param "+params+" -ofolder "+o_folder + "/spinalcord/"+" -x spline -g 1 -r 1"
            print(string)
            os.system(string)
        
        
        # # Calculating the slice-wise absolute average moco estimate 
        params_txt=os.path.dirname(moco_file) + '/moco_params_abs.txt'
        if abs_moco_param==True and not os.path.exists(params_txt):
            motion_abs={}
            for dim in ["x","y"]:
                img=nib.load(glob.glob(os.path.dirname(moco_file) + "/moco_params_"+dim+".nii.gz")[0])
                data_array = img.get_fdata()
                motion_abs[dim]=np.mean(np.abs(np.squeeze(data_array)),axis=0)
            
            moco_param_2D=np.column_stack((motion_abs["x"], motion_abs["y"])) 
            pd.DataFrame(moco_param_2D).to_csv(params_txt,index=False, header=None)

        # Read moco parameters (absolute or not)
        params_tsv=o_folder + "/spinalcord/" + 'moco_params.tsv'
        data=pd.read_csv(params_tsv, delimiter='\t')
        params_txt=params_tsv.split('.')[0] + '.txt'
        data.to_csv(params_txt,index=False, header=None)
        
        params_txt=params_tsv.split('.')[0] + '_abs.txt' if abs_moco_param else params_tsv.split('.')[0] + '.txt'
        params_data=pd.read_csv(params_txt, delimiter=',', header=None)
        # Plot moco parameters
        fig, axs = plt.subplots(1,1, figsize=(16, 5), facecolor='w', edgecolor='k')
        fig.tight_layout()
        fig.subplots_adjust(hspace = .5, wspace=.001)
        axs.plot(params_data[0]) # add axs[] if more than one run
        axs.plot(params_data[1])
        axs.set_title(ID + " " + task_name)
        
        axs.set_ylabel("Translation (mm)")
        axs.set_xlabel("Volumes")
        
        if not os.path.exists(params_txt.split(".")[0] + ".png") or redo==True:
            plt.savefig(params_txt.split(".")[0] + ".png")
        
        if plot_show==True:
            plt.show(block=False)
            # Calculate Framewise displacement (abs difference of displacement between each volumes)
            print('Framewise displacement, sub-' +ID)
            diff_X = np.abs(np.diff(params_data[0]))
            diff_Y = np.abs(np.diff(params_data[1]))
            meandiff=[np.mean(diff_X),np.mean(diff_Y)]
            print('Diff_X: ' + str(round(meandiff[0],3)) + ' mm')
            print('Diff_Y: ' + str(round(meandiff[1],3)) + ' mm')
            if not os.path.exists(o_folder + '/spinalcord/FD_mean.txt') or redo==True:
                np.savetxt(o_folder + '/spinalcord/FD_mean.txt', [meandiff])
                
                

        return moco_file, moco_mean_file
    
    def segmentation(self,ID=None,i_img=None,ctr_img=None,o_folder=None,task_name='',ses_name='',tag='',img_type="anat",redo=False,verbose=True):
        
        '''
        This function will segment the spinal cord:
            - GM from the anatomical image
            - GM + Wm from the func image
        
        Visual inspection and manual corrections are sometime needed
        https://spinalcordtoolbox.com/user_section/command-line.html#sct-propseg
        
        to do:
        - supress the directory /func/spinalcord/cord/
        
        Attributes:
        ----------
        ID: name of the participant
        i_img: input filename of functional images 3D (str, default:None, an error will be raise)
        ctr_img: Centerline image, need for img_type="func" (str, default:None)
        o_folder: output folder (str, default:None, the input folder will be used)
        ses_name: if the session have any specific name (ses- in BIDS format)
        task_name: if the run have any specific name (task- in BIDS format)
        img_type: the type of input image should be specify "func" or "anat"
        
        Outputs: 
        ----------
        
        '''
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if i_img==None:
            raise Warning("Please provide filename of the input file")
        
        if ctr_img==None and img_type=="func":
            raise Warning("Please provide centerline filename")
        
        #1. create ouput folder
        if o_folder is None : # gave the default folder name if not provided
            if img_type=="func":
                o_folder=self.preprocess_dir + "sub-" + ID+ "/" + ses_name + "/func/"+task_name+ "/" +self.config["preprocess_dir"]["func_seg"] + "/cord"
            elif img_type=="anat":
                o_folder=self.preprocess_dir + "sub-" + ID+ "/" + ses_name + "/anat/"+self.config["preprocess_dir"]["anat_sc_seg"]

        if not os.path.exists(o_folder):
            os.mkdir(o_folder)
        if img_type=="func" and not os.path.exists(o_folder+ "/spinalcord/"):
            os.mkdir(o_folder + "/spinalcord/")
            os.mkdir(o_folder + "/spinalcord/cord/")
        
        #2. run segmentation

        if img_type=="anat":
            o_img= o_folder +os.path.basename(i_img).split('.')[0] + "_seg.nii.gz"
            string= "sct_propseg -i " +i_img + " -c t1 -ofolder " +o_folder
           
        
        elif img_type=="func":
            o_img= o_folder + "/spinalcord/" +os.path.basename(i_img).split('.')[0] + "_seg.nii.gz"
            string="sct_deepseg_sc -i " +i_img +" -c t2s -centerline file -file_centerline "+ctr_img +" -o " +o_img
             
        if not os.path.exists(o_img) or redo==True:
            print(o_img)
            print(">>>>> Segmentation is running for the "+img_type + "image of the sub-" + ID + " ")
            os.system(string)
            print("Check the output and correct manually if needed" )
            print("fsleyes " + o_img)
        
        elif os.path.exists(o_img) and verbose==True:
            print(">>>>> Segmentation file already exists for the "+img_type + "image of the sub-" + ID + " ")
            print("fsleyes " + o_img)
        
        print(" ")    
        return o_img
    
    
            
    def label_vertebrae(self,ID=None,i_img=None,nb_labels=15,o_folder=None,ses_name='',task_name='',tag='',run_local=False,redo=False,verbose=True):
        
        '''
        Labelisation on the image
        
        https://spinalcordtoolbox.com/user_section/command-line.html#sct-label-utils
        
        Attributes:
        ----------
        ID: name of the participant
        i_img: input filename of functional images 3D (str, default:None, an error will be raise)
        nb_labels: total number of labels needed (int, default= 15)
        o_folder: output folder (str, default:None, the input folder will be used)
        ses_name: if the session have any specific name (ses- in BIDS format)
        task_name: if the run have any specific name (task- in BIDS format)
        run_local: to run the command on local laptop
        
        
        Outputs: 
        ----------
        
        '''
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if i_img==None:
            raise Warning("Please provide filename of the input file")
        
              
        #1. create ouput folder
        if o_folder is None : # gave the default folder name if not provided
            o_folder=self.preprocess_dir + "sub-" + ID+ "/"+ses_name+"/anat/"+self.config["preprocess_dir"]["anat_sc_seg"]

        if not os.path.exists(o_folder):
            os.mkdir(o_folder)

        #2. Create the labels
        o_img= o_folder+ os.path.basename(i_img).split(".")[0] +'_label-ivd_mask.nii.gz'
        if not os.path.exists(o_img) or redo==True:
            nb=np.array(range(1,nb_labels+1)) # array with label numbers
            string="sct_label_utils -i " +i_img + " -o "+o_img+" -create-viewer " + ', '.join(map(str, nb)).replace(" ","") #1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"
            print("Place labels at the posterior tip of each inter-vertebral disc for sub-" + ID)
            if run_local==True:
                print("Copy and path the command line locally, with corrected path: ")
                print(string)
            else:
                os.system(string)
        else:
            if verbose==True:
                print(">>>>> Label file already exists for sub-" + ID)
        print(" ")

        return o_img
    
    def coreg_anat2PAM50(self,ID=None,i_img=None,o_folder=None,seg_img=None,labels_img=None,img_type="t2",param=None,ses_name='',task_name='',tag='T2w',redo=False,verbose=True):
        
        '''
        Register anat to template and warp the template into the anatomical space
        
        https://spinalcordtoolbox.com/user_section/command-line.html#sct-register-to-template
        https://spinalcordtoolbox.com/user_section/command-line.html#sct-warp-template
        
        Attributes:
        ----------
        ID: name of the participant
        i_img <filename>: input filename of functional images 3D (str, default:None, an error will be raise)
        
        o_folder <folder dir>: output folder (str, default:None, the input folder will be used)
        ses_name <str>: if the session have any specific name (ses- in BIDS format)
        task_name <str>: if the run have any specific name (task- in BIDS format)
        img_type <filename>: Contrast to use for registration. (default: t2)
        param <list>: Parameters for registration default:
            "step=1,type=seg,algo=centermassrot:step=2,type=im,algo=syn,iter=5,slicewise=1,metric=CC,smooth=0"
        tag <str>: specify tag for the outputs. default T2w
        Outputs: 
        ----------
        
        '''
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if i_img==None or seg_img==None or labels_img==None:
            raise Warning("Please provide filename of the input file (i_img), segmentation file (seg_img) and inter disk labels(labels_img)")
        
        if param==None:
            param= "step=1,type=seg,algo=centermassrot:step=2,type=im,algo=syn,iter=5,slicewise=1,metric=CC,smooth=0"
        
        #1. create ouput folder
        if o_folder is None : # gave the default folder name if not provided
            o_folder=self.preprocess_dir + "sub-" + ID+ "/"+ses_name+"/anat/"+self.config["preprocess_dir"]["anat_sc_coreg"]

        if not os.path.exists(o_folder):
            os.mkdir(o_folder)

        #2. Coregistration between anat and template
        warp_from_anat2PAM50=o_folder+"/sub-"+ID+"_from-"+tag+"_to-PAM50_mode-image_xfm.nii.gz" # warping field form anat to PAM50
        warp_from_PAM502anat=o_folder+"/sub-"+ID+ "_from-PAM50_to-"+tag+"_mode-image_xfm.nii.gz" # warping field form PAM50 to anat
        if not os.path.exists(warp_from_anat2PAM50) or redo==True:
            string_coreg="sct_register_to_template -i "+i_img +" -s " + seg_img+" -ldisc " + labels_img +" -c " + img_type +" -param " + param +" -ofolder " +o_folder
            print(">>>>> Registration step is running for sub-" + ID)
            os.system(string_coreg)
            
            # rename the warping files
            os.rename(os.path.join(o_folder, "warp_anat2template.nii.gz"), warp_from_anat2PAM50)
            os.rename(os.path.join(o_folder, "warp_template2anat.nii.gz"), warp_from_PAM502anat)

            # 3. create the template in anat
            string_template="sct_warp_template -d " + i_img +" -w  " + warp_from_PAM502anat +" -s 1 -ofolder " +o_folder +"/template_in_" + tag + " -a 0 -s 0"

            os.system(string_template)

       
            if verbose==True:
                print("Done here: " + o_folder)
            

        else:
            if verbose==True:
                print(">>>>> Registration file anat2template already exists for sub-" + ID)
        
        
        print(" ")   
        return warp_from_PAM502anat, warp_from_anat2PAM50 
        
    def coreg_func2PAM50(self,ID=None,i_img=None,o_folder=None,func_seg=None,PAM50_cord=None,PAM50_t2=None,initwarp=None,initwarpinv=None,param=None,ses_name='',task_name='',redo=False,verbose=True):
        
        '''
        Register mean functional image to template
        
        https://spinalcordtoolbox.com/user_section/command-line.html#sct-register-multimodal
        
        Attributes:
        ----------
        ID: name of the participant
        i_img <filename>: functional mean image filename 3D (str, default:None, an error will be raise)
        func_seg <filename>: segmentation file of the functional image 3D
        PAM50_cord <filename>: PAM50 cord image cropped in restrictive shape to reduce computational time
        init_wrap <filename>: initial warping field (usually from PAM502anat)
        initwarpinv <filename>: initial inverse warping field (usually from anat2PAM50)
        
        o_folder <folder dir>: output folder (str, default:None, the input folder will be used)
        ses_name <str>: if the session have any specific name (ses- in BIDS format)
        task_name <str>: if the run have any specific name (task- in BIDS format)
        
        param <list>: Parameters for registration default:
            "step=1,type=imseg,algo=centermassrot,rot_method=hog,metric=MeanSquares,smooth=2:step=2,type=im,algo=syn,metric=CC,iter=3,slicewise=1"

        Outputs: 
        ----------
        
        '''
        
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if i_img==None or func_seg==None or initwarp==None or initwarpinv==None:
            raise Warning("Please provide filename of the input file (i_img), segmentation files (func_seg, PAM50_cord) and warping fields(init_wrap, initwarpinv)")
        
        if PAM50_cord==None:
            PAM50_cord=self.config["tools_dir"]["main_codes"] + "/template/"+self.config["PAM50_cord"]
            
        if PAM50_t2==None:
            PAM50_t2=self.config["tools_dir"]["main_codes"] + "/template/"+ self.config["PAM50_t2"]
            
        if param==None:
            param= "step=1,type=imseg,algo=centermassrot,rot_method=hog,metric=MeanSquares,smooth=2:step=2,type=im,algo=syn,metric=CC,iter=3,slicewise=1"
        
        #1. create ouput folder
        if o_folder is None : # gave the default folder name if not provided
            o_coreg=self.preprocess_dir + "sub-" + ID+"/"+ses_name+ "/func/"+ task_name+self.config["preprocess_dir"]["func_coreg"]
            o_folder= o_coreg+ "/spinalcord/"

        if not os.path.exists(o_coreg):
            os.mkdir(o_coreg)
        if not os.path.exists(o_folder):
            os.mkdir(o_folder)

            
        #2. Coregistration between anat and template
        o_img= o_folder+  os.path.basename(i_img).split('.')[0]+'_coreg_in_PAM50.nii.gz'
        o_warp_img=o_folder+ "/sub-"+ID+'_from-PAM50_to_func_mode-image_xfm.nii.gz'
        o_warpinv_img=o_folder+ "/sub-"+ID+'_from-func_to_PAM50_mode-image_xfm.nii.gz'
        
        if not os.path.exists(o_img) or redo==True:

            string_coreg="sct_register_multimodal -i "+PAM50_t2 +" -iseg " + PAM50_cord+" -d " + i_img+" -dseg " + func_seg  +" -param " + param +" -initwarp "+initwarp+ " -initwarpinv " + initwarpinv +" -owarp " + o_warp_img +" -owarpinv "+ o_warpinv_img +" -ofolder "+ o_folder + " -x spline"
            

            print(">>>>> Registration step is running for sub-" + ID)
            os.system(string_coreg)
            
            os.rename(o_folder+  os.path.basename(i_img).split('.')[0]+ "_reg.nii.gz",o_img) 
            #os.rename(o_folder+  'warp_' +os.path.basename(i_img).split('.')[0] + "2" + os.path.basename(PAM50_t2),o_warp_img)
                  
            if verbose==True:
                print("Done here: " + o_folder)
            

        else:
            if verbose==True:
                print(">>>>> Registration between func image and PAM50 already exists for sub-" + ID)
        
        

        
        print(" ")   
        return o_folder
        
    def apply_warp(self,ID=None,i_img=None,o_folder=None,dest_img=None,warping_field=None,ses_name='',task_name='',tag='_w',threshold=None,redo=False,verbose=True):
        '''
        Apply warping field to spinalcord input image in a destination image using sct_apply_transfo (sct toolbox)
        https://spinalcordtoolbox.com/user_section/command-line.html#sct-apply-transfo

        Attributes
        ----------
        ID: name of the participant
        i_img <filename>: input image filename 3D or 3D (str, default:None, an error will be raise)
        o_folder <folder dir>: output folder (str, default:None, the input folder will be used)
        dest_img <filename>: destination filename 
        warping_fields <filename>: Transformation(s), which can be warping fields (nifti image) or affine transformation matrix (text file). Separate with space.
                        
                 
        ses_name <str>: if the run have any specific name (ses- in BIDS format)
        task_name <str>: if the run have any specific name (task- in BIDS format)
        tag <str>: specify a tag for the output
        redo <Bolean> optional, to binarize the output file (default: False)
        redo <Bolean> optional, to rerun the analysis put True (default: False)
                
        return
        ----------
        o_img <filename>
                    
        '''
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if i_img==None or warping_field==None:
            raise Warning("Please provide filename of the input file (i_img), transformation files (warping_fields)")
        if o_folder==None:
            o_folder=os.path.dirname(warping_field) + "/"

        o_img= o_folder +  os.path.basename(i_img).split('.')[0] + tag + ".nii.gz"

        if verbose:
            print(" ")
            print(">>>>> Apply transformation to sub-" + ID)
        

        if not os.path.exists(o_img) or redo==True:
            # a. Apply transformation
            string1="sct_apply_transfo -i "+i_img+" -d "+dest_img+" -o "+o_img+" -w "+warping_field+" -x spline"
            os.system(string1)
            
                                                                                      
            if threshold:                                                    
                # b. Transform the output image in a binary image
                string2="fslmaths "+o_img+" -thr "+str(threshold)+" -bin " + o_img
                os.system(string2)
                
        else:
            if verbose:
                print("Tranformation was already applied put redo=True to redo that step")
            
        return o_img


    def csf_masks(self,ID=None,i_img=None,o_folder=None,task_name='',ses_name='',tag='',redo=False,verbose=True):
        
        '''
        Warp PAM50 csf into func space 
        
        Attributes:
        ----------
        ID: name of the participant
        i_img: input filename of input mask images 3D, should be a cord mask (str, default:None, an error will be raise)
        o_folder: output folder (str, default:None, the input folder will be used)
        ses_name: if the session have any specific name (ses- in BIDS format)
        task_name: if the run have any specific name (task- in BIDS format)
        img_type: the type of input image should be specify "func" or "anat"
        
        Outputs: 
        ----------
        
        '''
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if i_img==None:
            raise Warning("Please provide filename of the input file")
                
        #1. create ouput folder
        if o_folder is None : # gave the default folder name if not provided
            o_folder=os.path.dirname(i_img) + "/"
         
   
        if not os.path.exists(o_folder):
            os.mkdir(o_folder)
               # B. Use CSF+seg from T2 template image and check the output
        run_proc('fslmaths {} -add {} {}'.format(template_out + '/template/PAM50_cord.nii.gz',template_out + '/template/PAM50_csf.nii.gz', template_out + '/template/PAM50_csfwmgm.nii.gz'))

        #2. run dilatation of the mask
        o_dil1_img= o_folder +os.path.basename(i_img).split('.')[0] + "_dilated1.nii.gz"
        o_dil2_img= o_folder +os.path.basename(i_img).split('.')[0] + "_dilated2.nii.gz"
        o_csf_img= o_folder +os.path.basename(i_img).split('.')[0] + "_"+tag+".nii.gz"
        string1="fslmaths " +i_img +" -kernel 2D -dilM "+o_dil1_img
        string2="fslmaths " +o_dil1_img +" -dilM "+o_dil2_img
        string3="fslmaths " +o_dil2_img +" -sub "+i_img + " " +o_csf_img
      
        if not os.path.exists(o_csf_img) or redo==True:
            print(">>>>> Dilatation is running for the  sub-" + ID + " ")
            os.system(string1);os.system(string2)
            os.system(string3);
            print("Check the output and correct manually if needed" )
            print("fsleyes " + o_csf_img)
            os.remove(o_dil1_img);
        elif os.path.exists(o_csf_img) and verbose==True:
            print(">>>>> Surrounded file already exists for sub-" + ID + " ")
            print("fsleyes " + o_csf_img)
        
        print(" ")    
        return o_csf_img
