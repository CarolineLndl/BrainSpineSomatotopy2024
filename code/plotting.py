import os
import nibabel as nib
import numpy as np
import seaborn as sns
import glob

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap




from nilearn import plotting
from nilearn import surface

class Plotting:
    '''
    The Plotting class is used to manipulate and visualize maps
    Attributes
    ----------

    '''
    
    

# ======= SPINAL CORD ========
    
    
# ======= Brain ========
class Plot_brain:
    def __init__(self, config):
        self.config = config # Load config info

    def plot_colormap(colors=None,plot_colormap=True,redo=False,verbose=True):
        '''
        Plot your colormap
        Attributes
        ----------
        colors: list of color in html format (i.g. ["#1a04a4",'#0070ff','#07f6e0', "#9bff00",'#e8f703', '#fa9c03', '#ff3a00'])
        '''

        # Create a ListedColormap from the specified colors
        discretized_colormap = ListedColormap(colors)
       
        # Plot a colorbar to visualize the colormap
        if plot_colormap==True:
            plt.figure(figsize=(3,2))
            plt.imshow(np.arange(1, 8).reshape(1, -1), cmap=discretized_colormap, aspect='auto', vmin=1, vmax=7)
            plt.colorbar(ticks=np.arange(1, 8))
            plt.show() 
        
        return discretized_colormap

    def plot_3D(self, i_img=None,hemi_view=["lh","rh"],face_view=["lateral"],vmin=None,vmax=None,threshold=1e-6, mask_img=None,colormap='hot', tag="",output_dir=None, save_results=False):
        '''
        This function help to plot functional 3D maps on a render surface 
        Two nilearn function are used, see details here:
        https://nilearn.github.io/dev/modules/generated/nilearn.surface.vol_to_surf.html
        https://nilearn.github.io/dev/modules/generated/nilearn.plotting.plot_surf_roi.html
        
        to do:
        - Add option to import directly a volume
        - Add new value possibilities for view_per_line
        
        Attributes
        ----------
        i_img <filename>: filename of the input 3D volume image (overlay) should be specify, default: False
        hemi_view <str>: hemisphere of the brain to display, one of the three options should be specified: ["lh","rh"] or ["lh"] or ["rh"]
        face_view <str> : must be a string in: "lateral","medial","dorsal", "ventral","anterior" or "posterior" multiple views can be specified simultaneously (e.g., ["lateral","medial"])
        vmin <float>, optional, Lower bound of the colormap. If None, the min of the image is used.
        vmax <float>, optional, upper bound of the colormap. If None, the min of the image is used. 
        threshold <int> or None optional. If None is given, the image is not thresholded. If a number is given, it is used to threshold the image: values below the threshold are plotted as transparent. If “auto” is given, the threshold is determined magically by analysis of the image. Default=1e-6.
        colormap <str> : specify a colormap for the plotting, default: 'automn'
        tag <str> : specify a tag name for the output plot filenames, default: '' 
        output_dir <directory name>, optional, set the output directory if None, the i_img directory will be used
        save_results <boolean>, optional, set True to save the results
        
        '''
        
        if i_img==None:
            raise Warning("Please provide the filename of the overlay volume (ex: i_img='/my_dir/my_func_img.nii.gz')")
        if output_dir==None:
            output_dir=os.path.dirname(i_img) + "/plots/"
        
        if save_results:
            if not os.path.exists(output_dir):
                os.mkdir(output_dir) # create output folder
                os.mkdir(output_dir + "/tmp/") # create a temporary folder
     
        # 1. Select surface image for background --------------------------------------------------------
        surface_dir = self.config['templates']['surface_dir']
        
        img_surf={}
        for hemi in hemi_view:
            #2. Transform volume into surface image --------------------------------------------------------------------
            # include a mask is there are 0 values that you don't want to include in the interpolation
            img_surf[hemi]=surface.vol_to_surf(i_img,surface_dir+ hemi + ".pial",radius=0, 
                                     interpolation='nearest', kind='line', n_samples=10, mask_img=mask_img, depth=None)
  

            #3. Plot surface image --------------------------------------------------------------------
            side = "left" if hemi == "lh" else "right"
            for face in face_view:
                colorbar=True if face == face_view[-1] else False
                plot=plotting.plot_surf_roi(surface_dir+ hemi +".inflated", roi_map=img_surf[hemi],
                                            cmap=colormap, colorbar=colorbar,mask_img=mask_img,
                                            
                                            hemi=side, view=face,vmin=vmin,vmax=vmax,threshold=int(threshold),
                                            bg_map=surface_dir + hemi +".sulc",darkness=.7)

                if save_results:
                    
                    # Save each plot individually
                    plot.savefig(os.path.join(output_dir + "/tmp/", f'plot_{tag}_{hemi}_{face}.png'),dpi=150)
                    plt.close()

        if save_results:
            # Compute number of columns/rows and prepare subplots accordingly
            view_per_line= len(face_view)
               
            total_rows = ((len(face_view)*len(hemi_view))//view_per_line)
      
            fig, axs = plt.subplots(total_rows, view_per_line, figsize=(10*view_per_line, 8*total_rows))
                    
            for row, hemi in enumerate(hemi_view):
                for col, face in enumerate(face_view):
                    img = plt.imread(os.path.join(output_dir + "/tmp/", f'plot_{tag}_{hemi}_{face}.png'))
                    axs[row,col].imshow(img)
                    axs[row,col].axis('off')
        
            # Save the combined figure
            combined_save_path = os.path.join(output_dir, tag+".pdf")
            plt.savefig(combined_save_path, bbox_inches='tight')
            plt.show()
            
            #remove temporary files
            for hemi in hemi_view:
                for face in face_view:
                    os.remove(os.path.join(output_dir + "/tmp/", f'plot_{tag}_{hemi}_{face}.png'))
    
    
    
    