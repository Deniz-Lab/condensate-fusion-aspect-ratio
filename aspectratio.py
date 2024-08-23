#aspectratio.py
#File imports the DropletFusion class which takes an image stack array and operates on it.
#TO SHARE version of the file

#Written for Chawla, Tom, et al.  by Jenna Tom in Deniz Lab.
#Last updated: July 31, 2024

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import skimage.io
import math
import pandas as pd
import os


from scipy.spatial import distance
from scipy.optimize import curve_fit
from skimage import img_as_float
from skimage.filters import threshold_multiotsu
from skimage.morphology import remove_small_objects, binary_closing
from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import disk

class DropletFusion:
    '''Handles single instance of fusion via stacked image array in put & list of 
    properties recorded. 
    ---
    Parameters:
    image_stack_array = numpy array of image (output of skimage.io.imread())
    properties = categories of regionprops that are reported at end. 
    	Min required: ['area', 'axis_major_length', 'axis_minor_length']
    
    (optional)
    time_int = time interval between frames 
    '''
    def __init__(self,image_stack_array, properties, **kwargs):
        self.im = image_stack_array #numpy array of image (already read in)
        self.properties = properties
        
        self.time_int = kwargs.get('time_int',0.39521)
#         self.pixel_size = kwargs.get('pixel_size', 6.0222)

        #initialize attributes array
        self.values = np.zeros([len(image_stack_array),len(properties)+2])
        self.aspect_ratio = {}

        self.thresh = []
        self.binary = []
        self.clean = []
        self.labels = []
        self.centroid_list = []
        self.figs = {}
    
    
    #Image Pre-processing
    def threshold_image(self):
        for im_array in self.im:
            threshold = skimage.filters.threshold_otsu(im_array)
            self.thresh.append(threshold)
    
    def binary_image(self):
        if len(self.thresh) > 0:
            for im, threshold in zip(self.im, self.thresh):
                self.binary.append(im > threshold)
        else:
            print("No thresholds generated")
            
    #This combines thresholding, binary image formation, and cleaning if prev. functions
    #not invoked previously
    def clean_image(self):
        if len(self.thresh)<1:
            for im_array in self.im:
                image = im_array
                threshold = skimage.filters.threshold_otsu(image)
                binary = image > threshold
                
                self.thresh.append(threshold)
                self.binary.append(binary)
                
                self.clean.append(
                    skimage.morphology.remove_small_holes(
                        skimage.morphology.remove_small_objects(
                            skimage.segmentation.clear_border(binary)))
                )
        else:
            for binary_image in self.binary:
                self.clean.append(
                    skimage.morphology.remove_small_holes(
                        skimage.morphology.remove_small_objects(
                            skimage.segmentation.clear_border(binary_image)))
                )
    
    
    def find_regions(self, i, labels):
        '''Identifies the region to analyze in a single image frame. 
        --
        Parameters: 
        image_index (int): index of image in stack
        labels: labeled region version of image (output of skimage.measure.label(image))'''
        image_index = i
        labels = labels
        
        found_region = False #default to no region found until code identifies a region
        num_labeled_regions = len(skimage.measure.regionprops(labels))
        label_index = "None"
        
        #If only one labeled region labeled by skimage, then default to analyzing that region
        if num_labeled_regions == 1:
            found_region = True
            label_index = 0
        
        #If multiple regions, identify one to use
        elif num_labeled_regions >1:
            
            # Compares to previous image centers and areas if not the first frame
            if image_index != 0:
                area_comp_val = 0
                #k is the index of a previous image. Starts at one image prior and works back
                k = image_index - 1 
                area_index = self.properties.index("area") + 1
                
                #looks for first non-zero previous area 
                #(If segmentation is breaking and no images are labeled, 
                # area will be 0 by default)
                while area_comp_val == 0:
                
                    if self.values[k,area_index] > 0:
                        area_comp_val = self.values[k,area_index]
                    elif k < 1: #stops loop if no prev. images to compare to. Sets k = 0
                        break
                    else:
                        k -= 1
                        
            #Compares area & centroid position. 
            #Area must be within +/- 25% of the previous area (if condensate appears 
            #and is much larger or smaller, should be ignored)
            #Centroid must be within 4
            area_max = 0
            for j, region in enumerate(skimage.measure.regionprops(labels)):

                if image_index > 0 and (region.area > .75*area_comp_val or region.area < 1.25*area_comp_val):
                    centroid_comp_val = self.centroid_list[k]
                    if centroid_comp_val != 0:
                        if distance.cdist([region.centroid], [centroid_comp_val], 'euclidean')[0][0] < 4:
                            label_index = j
                            found_region = True
#                         else:
#                             print(j, k, "centroid condition not met. Dist = {dist}".format(dist = distance.cdist([region.centroid], [centroid_comp_val], 'euclidean')[0][0]))
                elif image_index == 0:
                    if region.area > area_max:
                        area_max = region.area
                        found_region = True
                        label_index = j
        
        if isinstance(label_index, int):
            self.labels += skimage.measure.regionprops(labels)[label_index]
        else:
            self.labels += "None"
        return (found_region, label_index)
    
    def valid_properties(self):
        valid_properties = True
        if "area" not in self.properties or "axis_major_length" not in self.properties or "axis_minor_length" not in self.properties:
            print("Properties List invalid. Requires a minimum of area, axis_major_length, and axis_minor_length to run")
            valid_properties = False
        return valid_properties
    
    def calc_properties(self, attribute_string):
        '''Calculates properties of a particular image_stack and fills .values array with index & properties.
        ---
        Parameters:
        attr_str (str): specifies image stack to analyze. Requires attr_str to be either clean or binary '''
        i = 0
        attr_str = attribute_string.lower()
        
        #specifies the array type analyzed
        if attr_str == "clean" or attr_str == "binary":
            array_stack = getattr(self,attr_str)
        else:
            print("No analysis possible. Specify image type clean or binary.")
            array_stack = np.zeros(shape = (0,0))
        
        #for a single image array from the array stack, 
        #labels image and uses find_regions to select a region to analyze
        for image in array_stack:
            labels = skimage.measure.label(image)
            found_region, label_index = self.find_regions(i, labels)

            if found_region == True:
                region = skimage.measure.regionprops(labels)[label_index]
                
                #limits region area to be large enough
                if region.area >=100:
                    
                    #fills np array of values for given properties with index, properties
                    list_new = [getattr(region,prop) for prop in self.properties]
                    self.values[i,:-1] = [i] + list_new
                    
                    #adds centroid of the region to centroid_list or sets it to zero if none
                    self.centroid_list += [region.centroid]
                else:
                    self.centroid_list += [0]
            else:
                self.centroid_list += [0]
        
            i += 1
    
    def calc_aspect_ratio(self):
        '''Calculates aspect ratio based on major and minor axis lengths'''
        index_major = self.properties.index("axis_major_length")
        index_minor = self.properties.index("axis_minor_length")
        with np.errstate(divide='ignore',invalid='ignore'):
            self.values[:,-1] = self.values[:,index_major+1]/self.values[:,index_minor+1]
        
    def count_no_nan(self,array_col):
        '''Counts the number of values in an array col that are not NaN'''
        return array_col.size - np.isnan(array_col).sum()
    
    def find_zero(self,array):
        '''Approximates the index of self.values where fusion is starting by identifying the first value with 3 subsequent decreasing points'''
        diff_array = np.diff(array[:,-1],1)
        zero_test = False
        i = 0
        num_decrease_points = 4 
        
        #iterates until zero is identified or until it loops through data set
        while zero_test == False:
            
            # If the aspect ratio value in the array is not nan for a given index, look for consecutive decreasing vals 
            if np.isnan(array[i,-1]) == False:
                num_compensate_nan = 0

                comp_array = array[i:i+num_decrease_points+num_compensate_nan, -1]
                num_no_nan = self.count_no_nan(comp_array)
                max_num_comp = array[i:, -1].size - num_compensate_nan - num_decrease_points
                
                #increases index proportionally to number of nan without exceeding number of points remaining in data set
                while num_no_nan < num_decrease_points and max_num_comp > 0:
                    num_compensate_nan += 1
                    max_num_comp -= 1
                    comp_array = array[i:i+num_decrease_points+num_compensate_nan, -1]
                    num_no_nan = self.count_no_nan(comp_array)

                comp_array = comp_array[~np.isnan(comp_array)]
                diff_no_nan_array = np.diff(comp_array)
                
                #checks magnitude of decreases
                if np.all(diff_no_nan_array<0):

                    #Condition that decrease needs to be greater than 12%, 8% and 4% 
                    #for three consecutive points once subtracting 1 from AR AND that the 
                    #starting point has to be higher than AR 1.4 to be a valid zero point
                    
                    #AR[0] > 1.4 condition prevents small magnitude changes 
                    ## around the baseline from being picked up. 
                    
                    #12%, 8%, 4% are hard coded here, but are set based on how this 
                    ## dataset worked. These can be manually adjusted
                    
                    #(AR - 1) is used for the denominator rather than (AR) to account for +1 shift of baseline
                    if comp_array[0] > 1.4 and -1*diff_no_nan_array[0] > (comp_array[0]-1)*.12 and -1*diff_no_nan_array[1] > (comp_array[1]-1)*.08 and -1*diff_no_nan_array[2] > (comp_array[2]-1)*.04:
                        zero_test = True
                        index_zero = i
                        self.aspect_ratio['index_zero'] = index_zero
                        return i
                    else:
                        i+=1
                        if i > len(diff_array)-(num_decrease_points-1):
                            self.aspect_ratio['index_zero'] = np.nan
                            zero_test = True
                else:
                    i+=1
                    if i > len(diff_array)-(num_decrease_points-1):
                        self.aspect_ratio['index_zero'] = np.nan
                        zero_test = True
            else:
                i+=1
                if i > len(diff_array)-(num_decrease_points-1):
                    self.aspect_ratio['index_zero'] = np.nan
                    zero_test = True
                    
    
    def remove_jumps(self,array):
        '''removes single instance jumps > 0.4. Requires that jump invert direction (eg. increase then decrease immediately after or vice versa)'''
        array_no_nan = array[:,-1][~np.isnan(array[:,-1])]
        diff_array_skip_nan = np.diff(array_no_nan)
        index_jump = np.argwhere(np.absolute(diff_array_skip_nan) > 0.4)
        
        #if there are indices found where there are jumps > 0.4, tests for sign inversion. Identifies outlier and sets AR to NaN
        if index_jump.size > 0:
            for i in index_jump:
                if i < diff_array_skip_nan.size-1:
                    if np.sign(diff_array_skip_nan[i]) != np.sign(diff_array_skip_nan[i+1]) and np.absolute(diff_array_skip_nan[i+1]) > 0.4:
                        index_outlier = np.where(array[:,-1] == array_no_nan[i+1])
                        array[index_outlier,-1] = np.nan
                elif i == diff_array_skip_nan.size-1:
                    index_outlier = np.where(array[:,-1] == array_no_nan[i+1])
                    array[index_outlier,-1] = np.nan
        return array
               
    
    def correct_time(self,index_zero, array):
        '''Adjusts time using index col and multiplying by frame by frame time interval attribute'''
        new_array = array.copy()
        
        #Reassigns NaN index_zero to 0 temporarily for the index assignment
        if type(index_zero) != int:
            print("no zero found")
            index_zero = 0

        new_array[:,0] = (new_array[:,0] - new_array[index_zero,0])*self.time_int
        return new_array
    
    def two_param_fit(self, t, a0, tau):
        '''Fitting function that has two parameters: a0 and tau. 
        a0 represents the aspect ratio at time 0, tau is the time constant of decay'''
        return 1 + (a0-1)* np.exp(-t/tau)
    
    def fit_relaxation(self, func):
        '''Fits processed data to a specified function '''
        index_zero = self.find_zero(self.remove_jumps(self.values))
#         print(index_zero)
        full_array = self.correct_time(index_zero, self.values)
        new_array = full_array[index_zero:,].copy()
        new_array = new_array[~np.isnan(new_array).any(axis=1)]
        func = getattr(self,func)
        
        xdata = new_array[:,0]
        y = func(xdata, new_array[0,4], 0.1)
        ydata = new_array[:,-1]

        popt, pcov = curve_fit(func, xdata, ydata, bounds = ([0.,0.],[3.,100.]))
        
        self.aspect_ratio["popt"] = popt
        self.aspect_ratio["pcov"] = pcov
        self.aspect_ratio["data"] = np.stack((xdata,ydata))
        self.aspect_ratio["n_fit"] = len(xdata)
        self.aspect_ratio['majorAxis0'] = full_array[index_zero,self.properties.index("axis_major_length")+1]

    
    def graph_fit(self, func, save_graph = False, file_name = "test_fit", file_type = "svg", surpress_show = False):
        '''Graphs fit calculated and shows plot'''
        func = getattr(self,func)
        if "data" not in self.aspect_ratio.keys():
            print("Data not fit. Run obj.fit_relaxation() first")
        else: 
            index_zero = self.aspect_ratio["index_zero"]
            full_array = self.correct_time(index_zero, self.values)
            
            xdata, ydata = self.aspect_ratio["data"]
            popt = self.aspect_ratio["popt"]
            
            fig = plt.plot(full_array[:,0], full_array[:,-1], 'go', label='all_data')
            plt.plot(xdata, ydata, 'bo', label='fit_data')
            
            plt.plot(xdata, func(xdata, *popt), 'r-',
                 label='fit: a0=%5.3f, tau=%5.3f' % tuple(popt))
            
            plt.xlabel('Time(s)')
            plt.ylabel('Aspect Ratio')
            plt.legend()
            
            if save_graph == True and file_type.lower() in ["svg","png","pdf"]:
                path = os.getcwd() + "/PolyP_Fusion_Figures"
                if os.path.isdir(path) == False:
                    os.mkdir(path)
                export_string = "PolyP_Fusion_Figures/{filename}_fit.{filetype}".format(filename = file_name, filetype = file_type.lower())
                plt.savefig(export_string)
                
            if surpress_show == False:
            	plt.show()
            else:
            	plt.close()
            
    def add_to_fig_dict(self, key_name, figure):
        if key_name not in self.figs.keys():
            self.figs[key_name] = figure
        else:
            overwrite = input("A figure already exists for {name}. Overwrite? (Y/N)".format(name = key_name))
            if overwrite == "Y":
                self.figs[key_name] = figure
                print("Overwritten")
            else:
                print("You typed {response}. The file was not overwritten. If this was not the desired result, re-run and type Y to overwrite".format(response = overwrite))

        

    #Validation/Visualization Functions
    def threshold_fig(self, surpress_show = False):
        if len(self.binary) < 1:
            print("No binary images to threshold")
        else: 
            num_frames = len(self.binary)
            
            im_dim_x = self.im[0].shape[1]
            im_dim_y = self.im[0].shape[0]

            n_col = 3
            n_row = num_frames
            
            fig_height = int(round(2*n_row*im_dim_y/30))
            fig_width = int(round(4*im_dim_x/40 + 4))

            fig, axs = plt.subplots(n_row, n_col, figsize = (fig_width,fig_height))
            fig.tight_layout()
            

            for n, (image, threshold, binary) in enumerate(zip(self.im, self.thresh, self.binary)): 

                axs[n,0].imshow(self.im[n],plt.cm.gray)# add a new subplot iteratively
                axs[n,0].set_title('Original: {im_num}'.format(im_num = n))
                axs[n,0].axis('off')

                axs[n,1].hist(image.ravel(), bins=256)
                axs[n,1].set_title('Histogram')
                axs[n,1].axvline(threshold, color='r')

                axs[n,2].imshow(binary, cmap=plt.cm.gray)
                axs[n,2].set_title('Thresholded')
                axs[n,2].axis('off')
            
            self.add_to_fig_dict("thresholding",fig)
            
            if surpress_show == False:
                plt.show()
            else:
                plt.close()
        
    def image_fig(self,attribute_str, surpress_show = False):
        
        image_set = getattr(self,attribute_str)
        num_images = len(self.im)
        
        if len(image_set) < 1:
            print("Image requested not processed")
        
        else:
            im_dim_x = self.im[0].shape[1]
            im_dim_y = self.im[0].shape[0]
            
            n_col = int(round(5*40/im_dim_x))
            n_row = math.floor(num_images/n_col)+1
            
            fig_width = 15 / 40 * im_dim_x
            fig_height = 12 / 30 * im_dim_y / 4 * n_row
            
            fig = plt.figure(figsize=(fig_width, fig_height))
        
            
            for n in range(num_images): 

                ax = plt.subplot(n_row, n_col, n + 1) 

                ax.imshow(image_set[n],plt.cm.gray)

                ax.set_title(str(n)) # chart formatting
                
            self.add_to_fig_dict(attribute_str, fig)
            
            if surpress_show == False:
                plt.show()
            else:
                plt.close()
    
    def segmentation_overlay(self, attribute_str = "clean", surpress_show = False):
        '''Generates segmentation overlay and shows original and intermediate images,
        as well as the final segmentation.
        
        ---
        Parameters: 
        attribute_str (default = clean): string name of image type to be processed. 
            (Options = "binary", "clean")
        surpress_show (default = False): If True, will not show plot--only will add to .graphs dict
            (Note: surpress show = True, speeds up program run time)'''
        
        
        image_set = getattr(self,attribute_str.lower())
        
        if len(image_set) < 1:
            print("No images of that type ({im_type}) to segment".format(im_type = attribute_str))
        elif len(self.binary) < 1:
            print("No binary images. Segmentation operates best with binary im type")
        else: 
            
            num_frames = len(image_set)
            
            im_dim_x = self.im[0].shape[1]
            im_dim_y = self.im[0].shape[0]

            n_col = 5
            n_row = num_frames
            
            fig_height = int(round((2*n_row*im_dim_y/30)/1.5))
            fig_width = int(round(10*im_dim_x/40 + 4))

            fig, axs = plt.subplots(n_row, n_col, figsize = (fig_width,fig_height))

            fig.tight_layout()
            
            # loop through the length of tickers and keep track of index
            for n, (orig, threshold, binary, image) in enumerate(zip(self.im, self.thresh, self.binary, image_set)): 

                image_label_overlay = skimage.color.label2rgb(image, image=image, bg_label=0)

                axs[n,0].imshow(orig,plt.cm.gray)# add a new subplot iteratively
                axs[n,0].set_title('Original: {im_num}'.format(im_num = n))
                axs[n,0].axis('off')

                axs[n,1].hist(orig.ravel(), bins=256)
                axs[n,1].set_title('Histogram')
                axs[n,1].axvline(threshold, color='r')

                axs[n,2].imshow(binary, cmap=plt.cm.gray)
                axs[n,2].set_title('Thresholded')
                axs[n,2].axis('off')
                
                axs[n,3].imshow(image, cmap=plt.cm.gray)
                axs[n,3].set_title(attribute_str)
                axs[n,3].axis('off')
                
                axs[n,4].imshow(image_label_overlay, cmap=plt.cm.gray)
                axs[n,4].set_title('Fit')
                axs[n,4].axis('off')
                
                labels = skimage.measure.label(image)
                found_region, label_index = self.find_regions(n, labels)
                if found_region == True:
                    region = skimage.measure.regionprops(labels)[label_index]
                    if region.area >=100:

                        orientation = region.orientation

                        y0, x0 = region.centroid

                        x1 = x0 + math.cos(orientation) * 0.5 * region.axis_minor_length
                        y1 = y0 - math.sin(orientation) * 0.5 * region.axis_minor_length
                        x2 = x0 - math.sin(orientation) * 0.5 * region.axis_major_length
                        y2 = y0 - math.cos(orientation) * 0.5 * region.axis_major_length

                        axs[n,4].plot((x0, x1), (y0, y1), '-r', linewidth=2.5)
                        axs[n,4].plot((x0, x2), (y0, y2), '-r', linewidth=2.5)
                        axs[n,4].plot(x0, y0, '.g', markersize=15)

                        deg_ang = -orientation/math.pi*180+90
                        ellipse_fit = mpatches.Ellipse(
                            xy = (x0,y0),
                            width = region.axis_major_length,
                            height = region.axis_minor_length, 
                            angle = deg_ang,
                            fill = None,
                            edgecolor='yellow',
                            linewidth=1,
                        )
                        axs[n,4].add_patch(ellipse_fit)

                        minr, minc, maxr, maxc = region.bbox
                        bx = (minc, maxc, maxc, minc, minc)
                        by = (minr, minr, maxr, maxr, minr)
                        axs[n,4].plot(bx, by, '-b', linewidth=2.5)
            fig_name = "segmentation-" + attribute_str
            
            self.add_to_fig_dict(fig_name,fig)
            
            if surpress_show == False:
                plt.show()
            else:
                plt.close()
        
    
    def export_figs(self,figs_to_export, file_name = "test", file_type = "png"):
    	'''Exports figures in .graphs dict as png, svg, or pdf.
    	
    	Parameters:
    	figs_to_export(str or list): 'name' of figure type to export. 
    		Inputs recognized: 
    		(1) 'all' -> all images in dict exported
    		(2) 'none', 'n', 'no', 'null' -> none (Alternatively though, do not invoke export figs)
    		(3) a specific image str -> 'thresholding', 'segmentation-[analysisVar]','[imageset]'
    		(4) unspecified str -> manual prompt/user approval 
    	file_name (str): name wanted to further identify specific files after saving 
    		(eg: "20240608_22_E_200-249_demo" will precede the specific version of fig exported)
    	file_type (str): specifies save file type (svg, png, pdf)'''
    	export_figs = []
    	if isinstance(figs_to_export,list):
        	for fig_type in figs_to_export:
        		if fig_type.lower() in self.figs.keys():
        			export_figs.append(fig_type.lower())
        		else:
        			print("File does not match internal names. Manually approve by typing Y or N for file(s) to save")
        			for key in self.figs.keys():
        				response_valid = False
        				while response_valid == False:
        					response = input("{name}? (Y/N)".format(name = key))
        					if response.lower() == "y" or response.lower() == "n":
        						response_valid = True
        						if response.lower() == "y":
        							export_figs.append(key)
        					else:
        						print("Response not recognized")
    	elif isinstance(figs_to_export,str):
        	if figs_to_export.lower() == "all":
        		export_figs = self.figs.keys()
        	elif figs_to_export.lower() in ["n","none","no","null"]:
        		print("Fig Exporting Surpressed")
        	elif figs_to_export.lower() in self.figs.keys():
        		export_figs.append(figs_to_export.lower())
        	else:
        		print("File does not match internal names. Manually approve by typing Y or N for file(s) to save")
        		for key in self.figs.keys():
        			response_valid = False
        			while response_valid == False:
        				response = input("{name}? (Y/N)".format(name = key))
        				if response.lower() == "y" or response.lower() == "n":
        					response_valid = True
        					if response.lower() == "y":
        						export_figs.append(key)
        				else:
        					print("Response not recognized")
        
    	for fig_name in export_figs:
        	if file_type.lower() in ["svg","png","pdf"]:
        		path = os.getcwd() + "/PolyP_Fusion_Figures"
        		if os.path.isdir(path) == False:
        			os.mkdir(path)
        	save_name = "PolyP_Fusion_Figures/{filename}_{figname}.{filetype}".format(filename = file_name, figname = fig_name, filetype = file_type.lower())
        	self.figs[fig_name].savefig(save_name)
