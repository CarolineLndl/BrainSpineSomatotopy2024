# Main imports ------------------------------------------------------------
import sys,os, gzip
import numpy as np
import nibabel as nib
import json

from nilearn import image

def tmean_img(ID=None,i_img=None,o_img=None,redo=False,verbose=False):
        
        '''
        This function will help to calculate mean images across volumes (tmean).
        use fslmaths
        
        Attributes:
        ----------
        ID: name of the participant
        i_img: input filename of functional images (str, default:None, an error will be raise), 4D image
        o_img: output folder name filename (str, default:None, the input filename will be used as a base)
        
        Outputs: 
        ----------
        Mean image inputfile_tmean.nii.gz
        '''
        if ID==None:
            raise Warning("Please provide the ID of the participant, ex: _.stc(ID='A001')")
        
        if i_img==None:
            raise Warning("Please provide filename of the input file")
      
        # Select the default output directory (input directory) 
        if o_img==None:
            o_img=i_img.split(".")[0] + "_tmean.nii.gz"

        # calculate the tmean:
        if not os.path.exists(o_img) or redo==True:
            string='fslmaths ' + i_img+ ' -Tmean '+ o_img
            os.system(string) # run the string as a command line
            
        if verbose ==True:
            print("Done : check the outputs files in fsleyes by copy and past:")
            print("fsleyes " + o_img)
            
        return o_img
    
def unzip_file(i_file,o_folder=None,ext=".nii",zip_file=False, redo=False,verbose=False):
        '''
        unzip the file to match with SPM
        Attributes
        ----------
        i_file <filename>: input file
        o_img: output folder name filename (str, default:None, the input filename will be used as a base)
        ext <str>: extension after unzip default: ".nii", put ".nii.gz" to zip a file
        zip_file <Bolean>: zip the file instead of unzip a file (default: False)
        redo <Bolean>: to rerun the analysis put True (default: False)
        
        return
        ----------
        o_file: <filename>: file name of unziped or zipped files 
        '''
        if o_folder is not None:
            output_file=o_folder + os.path.basename(i_file).split('.')[0] + ext
            
        else:
            output_file=i_file.split('.')[0] + ext
            
        # Zip file
        if zip_file:
            if not os.path.exists(i_file.split('.')[0] + ext) or redo:
                string= 'gzip ' + i_file
                os.environ(string)
                if o_folder:
                    os.rename(i_file.split('.')[0] + ext, output_file)
            else:
                if verbose:
                    print("Zip was already done please put redo=True to redo that step")
                else:
                    pass
                 
        else:
            if not os.path.exists(i_file.split('.')[0] + ext) or redo:
            
                input = gzip.GzipFile(i_file, 'rb') # load the  .nii.gz
                s = input.read(); input.close()
                unzip = open(i_file.split('.')[0] + ext, 'wb') # save the .nii
                unzip.write(s); unzip.close()
                os.rename(i_file.split('.')[0] + ext, output_file)
                
                if verbose:
                    print('Unzip done for: ' + os.path.basename(i_file))
                else:
                    pass
                
            else :
                if verbose:
                    print("Unzip was already done please put redo=True to redo that step")
                else:
                    pass
        
            
        return output_file
    
def standardize(i_img=None,o_folder=None,json_files=None,mask_img=None,tag="",redo=False,verbose=False):

        '''
        unzip the file to match with SPM
        Attributes
        ----------
        i_img <filename>, mendatory, default: None: input filename
        o_folder <dirname> optional, default None : output directory (e.g: output_file='/mydir/')
        json_file <str>: 
        mask_img <filename> optional, default None, If provided, signal is only standardized from voxels inside the mask. 
        redo <Bolean>: to rerun the analysis put True (default: False)
        
        '''
        
        if i_img==None:
            raise ValueError("Please provide the input filename, ex: _.cleam_images(i_img='/mydir/sub-1_filename.nii.gz')")
     
        timeseries=nib.load(i_img).get_fdata() # extract Time series dats
        signals= timeseries.reshape(-1, timeseries.shape[-1]).T # reshape timeseries (nb_volumes, nb_voxels)

        if signals.shape[0] == 1:
            warnings.warn('Standardization of 3D signal has been requested but '
                              'would lead to zero values. Skipping.')
        else:
            signals= timeseries.reshape(-1, timeseries.shape[-1]).T # reshape timeseries (nb_volumes, nb_voxels)
            std = signals.std(axis=0)
            std[std < np.finfo(np.float64).eps] = 1.  # avoid numerical problems
            signals /= std

        # save into filename
        o_filename=i_img.split('.')[0] + tag + ".nii.gz"
        json_file=o_filename.split('.')[0] + ".json"
        if not os.path.exists(o_filename) or redo==True:
            o_img=image.new_img_like(i_img, signals.T.reshape(timeseries.shape),copy_header=True)
            o_img.to_filename(o_filename) #save image

            if mask_img:
                
                string="fslmaths "+o_filename+" -mas " +mask_img +" "+ o_filename
                os.system(string)

            
            infos={"standardize":True,"mask":mask_img}
            with open(json_file, 'w') as f:
                json.dump(infos, f) # save info
                    
