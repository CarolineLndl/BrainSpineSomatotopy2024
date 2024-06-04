# Main imports ------------------------------------------------------------
import glob, os
import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np
from nilearn.image import math_img,smooth_img
from nilearn.input_data import NiftiMasker
from brsc_preprocess import Preprocess_Sc, Preprocess_Br

def tSNR(config=None,ID=None,i_img=None,mask=None,structure='brain',o_tag="_moco",inTemplate=False,redo=False):
    '''
        This function calculate the tSNR withi the brain or spinal cord
        
        Attributes:
        ----------
        config: load config file
        ID: participant ID
        i_img: 4d func image
        inTemplate: put True to coregister the tSNR map into template space
        redo: put True to re-run the analysis on existing file (default=False)
    
    '''

    if i_img==None:
        i_img= glob.glob(config["main_dir"] + config["moco_files"]["dir"].format(ID,structure) + config["moco_files"]["moco_f"])[0]
    if mask==None:
        mask=glob.glob(config["main_dir"]+ config["seg_files"]["dir"].format(ID,structure) +"/*")[0]

    native_tSNR= config["tSNR"]["dir"] + "sub-"+ ID +"/" + os.path.basename(i_img).split(".")[0] + "_tSNR.nii"
    # compute tSNR *******************************************************************************
    if not os.path.exists(native_tSNR) or redo==True:
        tsnr_func= math_img('img.mean(axis=3) / img.std(axis=3)', img=i_img)
        tsnr_func_smooth = smooth_img(tsnr_func, fwhm=[3,3,6])
        tsnr_func_smooth.to_filename(native_tSNR)

    # extract value inside the mask
    o_txt=config["tSNR"]["dir"] + "sub-"+ ID +"/sub-" + ID + "_" + structure +o_tag +"_tSNR_mean.txt"
    if not os.path.exists(o_txt) or redo==True:
        print("redo")
        print(o_txt)
        masker_stc = NiftiMasker(mask_img=mask,smoothing_fwhm=None,standardize=False,detrend=False) # select the mask
        tSNR_masked=masker_stc.fit_transform(native_tSNR) # mask the image
        mean_tSNR_masked=np.mean(tSNR_masked) # calculate the mean value
        with open(o_txt, 'w') as f:
            f.write(str(mean_tSNR_masked))  # save in a file

        
    # Coregistratin into MNi or PAM50 space **************************************************
    if inTemplate == True:
        #>>>>>>>>>>>> Spinal cord
        if structure=="spinalcord":
            preprocess_Sc=Preprocess_Sc(config)
            outreg_f=native_tSNR.split('.')[0] +"_inTemplate.nii.gz" #tSNR in template space
            if not os.path.exists(outreg_f):
                warp_img= glob.glob(config["main_dir"] + config["warp_files"]["dir_sc"].format(ID,structure) + config["warp_files"]["warp_func2PAM50"])[0] # warping field
                preprocess_Sc.apply_warp(ID=ID,
                                         i_img=native_tSNR +".gz",
                                         o_folder=os.path.dirname(native_tSNR)+ "/",
                                         dest_img=config["tools_dir"]["main_codes"] + "/template/"+config["PAM50_t2"],
                                         warping_field=warp_img,
                                         ses_name='',task_name='',tag='_inTemplate',
                                         threshold=None,redo=True,verbose=True) # apply transformation
            
        #>>>>>>>>>>>> Brain, in 2 steps
        if structure=="brain":
            preprocess_Br=Preprocess_Br(config)
           
            # >> Coreg into anat space
            coreg2anat_f=native_tSNR.split('.')[0] +"_inAnat.nii"
            if not os.path.exists(coreg2anat_f) or redo==True:
                # coregistration to anat space
                anat_f= glob.glob(config["main_dir"] + config["brain_anat"]["file"].format(ID))[0] 
                preprocess_Br.coreg_func2anat(ID=ID,ref_img=anat_f,
                                              source_img=native_tSNR,
                                              other_files=None,
                                              tag="_inAnat",
                                              o_folder=os.path.dirname(native_tSNR))
                  
            # >> Coregitration into MNI space
            if not os.path.exists(coreg2anat_f.split('.')[0] +"_inTemplate.nii.gz"):
                warp_f=glob.glob(config["main_dir"] + config["warp_files"]["dir_brain"].format(ID,structure) + config["warp_files"]["warp_anat2MNI"])[0]
                
                outreg_f=coreg2anat_f.split('.')[0] + "_inTemplate.nii" # default output
                preprocess_Br.normalisation(ID=None,
                                            warp_file=warp_f,
                                            coreg2anat_file=coreg2anat_f,
                                            o_file=outreg_f,
                                            brain_mask=None,
                                            redo=False)
            outreg_f=coreg2anat_f.split('.')[0] + "_inTemplate.nii.gz" #tSNR in template space
            
    return o_txt, native_tSNR, outreg_f if inTemplate == True else None

def tSNR_group(config=None,i_img=None,structure='brain',o_tag="_moco",redo=False):
    
    if i_img==None:
        raise Warning("Provide a list of filenames !")

    # Create 4D image
    o_4d_img=config["tSNR"]["dir"].split("tSNR/")[0] + "/group/" + "4d_n" + str(len(config['participants_IDs'])) +  "_"+structure+"_tSNR.nii.gz" # output filename
    all_files=(' ').join(i_img) # join strings
    
    if not os.path.exists(o_4d_img) or redo==True:
        string="fslmerge -t " + o_4d_img + " " + all_files
        os.system(string)

    # Calculate mean image
    o_mean_img=config["tSNR"]["dir"].split("tSNR/")[0] + "/group/" + "mean_n" + str(len(config['participants_IDs'])) +  "_"+structure+"_tSNR.nii.gz" # output filename
    if not os.path.exists(o_mean_img) or redo==True:
        string="fslmaths " + o_4d_img + " -Tmean " + o_mean_img
        os.system(string)

    return o_mean_img

             
            
          
    
    
def plot_metrics(config,df=None,y=None,index='ID',columns=['structure'],y_title="y_axis",save_plot=False):
        
        '''
        This function will help to plot different metrics
        
        Attributes:
        ----------
        df: dataframe with metrics informations
        y: values of the metrics to plot
        index: index column for each individuals
        columns: columns to plot as separate variables
            
        '''
                            
        intra_plot=False
        jitter=0
        df_counts_wide = df.pivot(index=index, columns=columns, values=y)
        df_counts_wide.reset_index(drop=True, inplace=True)
        df_counts_wide=df_counts_wide[["spinalcord","brain"]]
        columns = df_counts_wide.columns
        colors=['#20b5bf','#ebb80b']#'#efb537',
        form=['o','o']
        
        df_x_jitter = pd.DataFrame(np.random.normal(loc=0, scale=jitter, size=df_counts_wide.shape), columns=columns)
        df_x_jitter += np.arange(len(columns))
                
        fig, ax = plt.subplots(figsize=(4, 6))
                
        columns2list = []
        # Plot data points
        for col_nb in range(0,len(columns)):
            col=columns[col_nb]
            _ = ax.plot(df_x_jitter[col], df_counts_wide[col],form[col_nb],c=colors[col_nb],  zorder=1, ms=11, mew=1,alpha=0.5)#form[col_nb],,c=colors[col_nb])
            columns2list.append(col)
                    
            #Box plot
            box=plt.boxplot(df_counts_wide[col],positions=[col_nb],patch_artist=True, widths=(0.8),showfliers=False)
        
                  
            #set the box colors
            for patch, color in zip(box['boxes'], colors):
                patch.set_alpha(0.6)
                patch.set_facecolor(colors[col_nb])
                patch.set_edgecolor(colors[col_nb])
                patch.set_linewidth(2)
                    
                # set the outline colors
                for median in box['medians']:
                    median.set_color(colors[col_nb])
                    median.set_linewidth(2)
                        
                for caps in box['caps']:
                    caps.set_color(colors[col_nb])
                    caps.set_linewidth(2)
                        
                for whisker in box['whiskers']:
                    whisker.set_color(colors[col_nb])
                    whisker.set_linewidth(2)
             
            # Add lines between points
        if intra_plot==True:
            for i_row in range(0,len(df_counts_wide.index)):
                subject=df_counts_wide.index[i_row]
                for i_col in range(0, len(columns), 2):
                    _ = ax.plot(df_x_jitter.loc[i_row, columns[i_col:i_col+2]], df_counts_wide.loc[subject, columns[i_col:i_col+2]], color = 'grey', linewidth = 0.5, linestyle = '--', zorder=-1)
                
            # Add lines
            for i_row in range(0,len(df_counts_wide.index)):
                subject=df_counts_wide.index[i_row]
                for i_col in range(1, len(columns)+1, 2):
                    _ = ax.plot(df_x_jitter.loc[i_row, columns[i_col:i_col+2]], df_counts_wide.loc[subject, columns[i_col:i_col+2]], color = 'grey', linewidth = 0.5, linestyle = '--', zorder=-1)
        
             
        # Add axis ticks & labels
        _ = ax.set_xticks(range(len(columns)))
        _ = ax.set_xticklabels(columns2list)
        _ = ax.set_ylabel(y_title)
        #plt.ylim(0,0.3)

        if save_plot==True:
            plot_f=config["tSNR"]["group_results_dir"] + "/mean_brsc_n" + str(len(config["participants_IDs"])) + "_" + y + ".pdf"
            plt.show()
            fig.savefig(plot_f,dpi=300, bbox_inches='tight')
    
