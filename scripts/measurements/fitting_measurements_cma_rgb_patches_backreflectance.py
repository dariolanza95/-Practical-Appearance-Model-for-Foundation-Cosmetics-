# Example call: time python ./scripts/measurements/fitting_measurements_cma_rgb_patches.py -cosmetic_type Clarins105 -optimizer CMA-ES 
import os
import argparse
import json
import time

#import cma
import numpy as np
import itertools

from scipy.optimize import minimize
from skimage import io, color
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

# from skimage import data, img_as_float
# from skimage.transform import resize
# from skimage import io, color
# from skimage.metrics import structural_similarity as ssim

def mkdir(folder):
    if not os.path.exists(folder):
        os.mkdir(folder)

class MaterialHelper:

    def __init__(self, phi_i_res = 36, theta_i_res = 18, phi_o_res = 12, theta_o_res = 6, wavelength_res = 3,theta_o_max = 100,theta_o_min = 0,theta_i_min = 0,theta_i_max = 90):
        
        self.verbose = False
        self.phi_i_res = phi_i_res # 5 (72), 10 (36)
        self.theta_i_res = theta_i_res # 5 (36), 10 (18)

        self.phi_o_res = phi_o_res
        self.theta_o_res = theta_o_res
        
        self.wavelength_min = 380
        self.wavelength_max = 830
        self.wavelength_res = wavelength_res

        # Bins per dimension
        self.wavelengths_bins = np.linspace(self.wavelength_min, self.wavelength_max, self.wavelength_res)

        self.phi_i_bins = np.linspace(0,180, self.phi_i_res)
        self.theta_i_bins = np.linspace(theta_i_min, theta_i_max, self.theta_i_res)
        #self.phi_o_bins = np.linspace(0, 360, phi_o_res)
        self.phi_o_bins = np.linspace(-180,0, self.phi_o_res)
        self.theta_o_bins = np.linspace(theta_o_min, theta_o_max, self.theta_o_res)

        #np.set_printoptions(suppress=True)

    # Data binning raw data
    def loadData(self, file_path, output_file,backrefl_data = False):

        raw_data = np.loadtxt(file_path)
        n, m = raw_data.shape

        #print("Debugging Load Data")
        #print(raw_data.shape)
        #print(backrefl_data)
        #input()

        # --------- BSDF evaluation for incident direction and wavelength
        #data = np.zeros((self.phi_i_res+1, self.theta_i_res+1, self.phi_o_res+1, self.theta_o_res+1, self.wavelength_res), dtype=np.dtype(float))
        #data_count = np.ones((self.phi_i_res+1, self.theta_i_res+1, self.phi_o_res+1, self.theta_o_res+1, self.wavelength_res), dtype=np.dtype(int))

        # data = np.full((self.phi_i_res+1, self.theta_i_res+1, self.phi_o_res+1, self.theta_o_res+1, self.wavelength_res), np.zeros,dtype=np.dtype(float))
        # data_count = np.full((self.phi_i_res+1, self.theta_i_res+1, self.phi_o_res+1, self.theta_o_res+1, self.wavelength_res),np.zeros ,dtype=np.dtype(int))
        if self.wavelength_res != 3:
            print("WARNING wavelength_res is NOT 3!!")

        data = np.zeros((self.phi_i_res+1, self.theta_i_res+1, self.phi_o_res+1, self.theta_o_res+1, self.wavelength_res), dtype=np.dtype(float))
        data_count = np.zeros((self.phi_i_res+1, self.theta_i_res+1, self.phi_o_res+1, self.theta_o_res+1, self.wavelength_res), dtype=np.dtype(int))
        # print(data.shape)
        # data = np.zeros((self.phi_i_res, self.theta_i_res, self.phi_o_res, self.theta_o_res, self.wavelength_res), dtype=np.dtype(float))
        # data_count = np.zeros((self.phi_i_res, self.theta_i_res, self.phi_o_res, self.theta_o_res, self.wavelength_res), dtype=np.dtype(int))

        # Loop implementation is really slow, but we can start with this
        for i in range(n):

            phi_i = raw_data[i, 0]
            theta_i = raw_data[i, 1]
            phi_o = raw_data[i, 2]
            theta_o = raw_data[i, 3]
            wavelength = raw_data[i, 4]
            reflectance = raw_data[i, 5]

            phi_i_idx, theta_i_idx, phi_o_idx, theta_o_idx, _ = self.getIndices(phi_i, theta_i, phi_o, theta_o, wavelength, fixed_light = False)

            #custom wavelength_idx
            if i % 3 == 0:
                #red
                wavelength_idx = 0 
            if i % 3 == 1:
                #green
                wavelength_idx = 1
            if i % 3 == 2:
                #blue
                wavelength_idx = 2
            if self.verbose:
                print ("Query phi_i = {0}, theta_i = {1}, phi_o = {2}, theta_o = {3}, wavelength = {4}".format(phi_i, theta_i, phi_o, theta_o, wavelength))
                print ("Query phi_i_dx = {0}, theta_i_idx = {1}, phi_o_idx = {2}, theta_o_idx = {3}, wavelength_idx = {4}".format(phi_i_idx, theta_i_idx, phi_o_idx, theta_o_idx, wavelength_idx))
            # if reflectance > 5:
            #     print(reflectance)
            #     print ("Query phi_i = {0}, theta_i = {1}, phi_o = {2}, theta_o = {3}, wavelength = {4}".format(phi_i, theta_i, phi_o, theta_o, wavelength))
            #     print ("Query phi_i_dx = {0}, theta_i_idx = {1}, phi_o_idx = {2}, theta_o_idx = {3}, wavelength_idx = {4}".format(phi_i_idx, theta_i_idx, phi_o_idx, theta_o_idx, wavelength_idx))
            data[phi_i_idx, theta_i_idx, phi_o_idx, theta_o_idx, wavelength_idx] += reflectance
            data_count[phi_i_idx, theta_i_idx, phi_o_idx, theta_o_idx, wavelength_idx] += 1
        
        data = data / data_count
        # data = np.nan_to_num(data)
        model_values = data
        model_values_pos = model_values[:,:,1,:,:]
        model_values_neg = model_values[:,:,2,:,:]
        model_values_neg = model_values_neg[:,:,1:,:]
        model_values_pos = model_values_pos[:,:,1:,:]
        model_values = np.concatenate(((model_values_neg),model_values_pos),axis=2)
        # working code    
        # helper_bxdf2 = MaterialHelper(phi_i_res = 0, theta_i_res = 91, phi_o_res = 0, theta_o_res = 181, wavelength_res = 81,theta_o_min=-90,theta_o_max=90,theta_i_min=0,theta_i_max=90) # Model sanity check!
        # helper_bxdf2 = MaterialHelper(phi_i_res = 0, theta_i_res = 5, phi_o_res = 0, theta_o_res = 181, wavelength_res = 81,theta_o_min=-90,theta_o_max=90,theta_i_min=30,theta_i_max=70) # Model sanity check!
        #helper_bxdf2 = MaterialHelper(phi_i_res = 0, theta_i_res = 6, phi_o_res = 0, theta_o_res = 181, wavelength_res = 81,theta_o_min=-90,theta_o_max=90,theta_i_min=30,theta_i_max=70) # Back reflectance!

        helper_bxdf2 = MaterialHelper(phi_i_res = 0, theta_i_res = 6, phi_o_res = 0, theta_o_res = 181, wavelength_res = 3,theta_o_min=-90,theta_o_max=90,theta_i_min=20,theta_i_max=70) #Back reflectance!

        #select only incident angles used in measuremnets
        incident_angles = [30,40,50,60,70]
        if backrefl_data:

            print("backrefl_data case")

            incident_angles = [20, 30, 40, 50, 60, 70]

            # final_data = np.full( ( len(incident_angles),180,3),np.nan)
            final_data = np.full( ( len(incident_angles),180,3),np.nan)
            for idx,incident_angle in enumerate(incident_angles):
                _, theta_i_idx, _, _, _ = helper_bxdf2.getIndices(0.0, incident_angle, 0.0, 0.0, 0.0, fixed_light = False)
                print("theta_i_idx {0} vs theta_i {1} ".format(theta_i_idx,incident_angle))
                # pos_vals = model_values[0,theta_i_idx,91:(90+81),:]
                #selecte reflectance values for the positive semiplane (reflection mode) light and camera are on opposite semiplanes
                pos_vals = model_values[0,theta_i_idx,91:(90+91),:]
                #selecte reflectance values for the negative semiplane (backreflection mode), i.e., light and camera are on the same semiplane
                neg_vals = model_values[0,theta_i_idx,0:90,:]
                # pos_vals = model_values[:,:,:,:]
                # pos_vals= data[:,theta_i_idx,91:,:]
                #working but need to be flipped
                #final_data[idx,:,:] = model_values[0,theta_i_idx,:(90*2),:]
                neg_vals = np.flip(neg_vals,axis=0)
                final_data[idx,:,:] = np.concatenate((neg_vals,pos_vals),axis=0)
                # breakpoint()
                
                #final_data[idx,:,:] = model_values[0,theta_i_idx,:(90*2),:]
            #exit()

            #np.save(output_file, final_data)
            #return final_data
        else:
            final_data = np.full( ( len(incident_angles),80,3),np.nan)
            for idx,incident_angle in enumerate(incident_angles):
                _, theta_i_idx, _, _, _ = helper_bxdf2.getIndices(0.0, incident_angle, 0.0, 0.0, 0.0, fixed_light = False)
                print("theta_i_idx {0} vs theta_i {1} ".format(theta_i_idx,incident_angle))
                # pos_vals = model_values[0,theta_i_idx,91:(90+81),:]
                # working code
                # pos_vals = model_values[0,theta_i_idx,91:(90+81),:]
                pos_vals = model_values[0,theta_i_idx,91:(90+81),:]
                
                # pos_vals = model_values[:,:,:,:]
                # pos_vals= data[:,theta_i_idx,91:,:]

                final_data[idx,:,:] = pos_vals
            
            # final_data[idx] = model_values[:,theta_i_idx,:,:]

        #also consider only positive angles

        #final shape of array should be 6x80x3
        np.save(output_file, final_data)


    
        return final_data

    # Data binning indices
    def getIndices(self, phi_i, theta_i, phi_o, theta_o, wavelength, fixed_light = False):

        # Fixed incident ligth
        if fixed_light:
            phi_i_idx = 0
            theta_i_idx = 0
        else:

            #phi_i -= 180

            phi_i_idx = np.digitize(phi_i, self.phi_i_bins)
            theta_i_idx = np.digitize(theta_i, self.theta_i_bins)
        phi_o_idx = np.digitize(phi_o, self.phi_o_bins)
        theta_o_idx = np.digitize(theta_o, self.theta_o_bins)
        wavelength_idx = np.digitize(wavelength, self.wavelengths_bins)

        return phi_i_idx, theta_i_idx, phi_o_idx, theta_o_idx, wavelength_idx

# Old working Code
# def getBSDFData(basename, helper, root_folder):
#     data_file =     os.path.join("./scripts/measurements/fitting/{0}.txt".format(basename))
#     output_file =   os.path.join("./scripts/measurements/fitting/{0}.npy".format(basename))
    
#     # root_folder =os.path.join("scripts","measurements","matching54") 
#     out_folder = os.path.join(root_folder, basename)
def getBSDFData(basename, helper, out_folder,output_file="",data_file="", backrefl_data = False):
    if data_file == "":
        data_file =     os.path.join("./scripts/measurements/fitting/{0}.txt".format(basename))
    if output_file == "":
        output_file =   os.path.join("./scripts/measurements/fitting/{0}.npy".format(basename))
    
    # root_folder =os.path.join("scripts","measurements","matching54") 
    # out_folder = os.path.join(root_folder, basename)
    #mkdir(root_folder)
    #mkdir(out_folder)

    #data = np.load("scripts/measurements/model_bsdf.npy")
    out_path = os.path.join(out_folder,"{0}.png".format(basename))
    helper.verbose = True
    
    data_test = helper.loadData(data_file, output_file, backrefl_data)
    # breakpoint()
    
    #print(getBSDFData)
    #print (data_test.shape)
    #input()

    # print("DATA TEST SHAPE")

    # print(data_test.shape)
    # print(data_test[0,:,:])
    #convert data from wavelengths to rgb
    
    #incident angles
    m = 6
    if backrefl_data:
        n = 115
    else:
    #observation angles
        n = 80
    if output_file.endswith(".npy"):
        data = np.load(output_file)
        # fake_data = np.full((m,n,3),1.0)
        # data = fake_data
        if data.shape != (m,n,3):
            print("WARNING!")
            err_msg = "Please provide a file with shape {m}x{n}x3, where m angle of incidence and n angle of observations x 3 rgb channels."
            print(err_msg)
            print(data.shape)        
    else:
        print("ERROR provided data are not a npy file.")
    return [data,data,out_folder,out_path]

def parseNormalizedReflectanceData(path,output_path,helper):
#first and third dimensions are for phi, which we will ignore.
#use the full resolution, we will ignore all the non-zero angles when plotting
    # refl_tensor = np.zeros((1,181,1,181,81))
    refl_tensor = np.full((1,181,1,181,81),np.nan)

    theta_angles = {}

    with open(path, 'r') as file:
        for line in file:
            data = line.strip().split()  
            #first element gives theta_i, second theta_0, the remaining are the spectral reflectance
            theta_i = int(data[0])
            theta_o = int(data[1])
            _, theta_i_idx, phi_o_idx, theta_o_idx, _ = helper.getIndices(0, theta_i, 0, theta_o, 0)
                

            if theta_i in theta_angles.keys():
                theta_angles[theta_i].append(theta_o)
            else:
                theta_angles[theta_i] = [theta_o]

            spectral_data = data[2:]
            if len(spectral_data) != 81:
                print("WARNING SPECTRAL DATA LENGTH IS NOT 81! It is {0}, this might have implications on how data are processed, be careful".format(spectral_data))
            
            # print("spectral data len: " + str(len(spectral_data)))
            # print("idx {0} vs data {1}".format(theta_i_idx,data[0]))
            # print("idx {0} vs data {1}".format(theta_o_idx,data[1]))
            x = np.array(spectral_data,dtype="float")
            refl_tensor[0,theta_i_idx,0,theta_o_idx,:] = x


    np.save(output_path,refl_tensor)

    return theta_angles

def getMeasuredDataRGB(cosmetic_type = "", backrefl_data = True):
    if cosmetic_type == "":
        cosmetic_type = "Clarins112"
    cosmetic_dict = {}

    # cosmetic_dict["Clarins112"] = "Foundation_112_5W/Clarins112_5W_thick_layer_tab_RGB"
    # cosmetic_dict["Clarins103"] = "Foundation_103N/Clarinsfoundation103N_norm_thick_layer_tab_RGB"
    # cosmetic_dict["Clarins105"] = "Foundation_105/Clarins_105_Nude_coll_norm_RGB"
    # cosmetic_dict["Clarins108"] = "Foundation_108_5W/Clarinsfoundation108_5W_norm_thick_layer_tab_RGB"
    
    cosmetic_dict["Clarins112"] = "Clarins112/Clarins112_5W_thick_layer_tab_CIE_RGB"
    cosmetic_dict["Clarins103"] = "Clarins103/Clarinsfoundation103N_norm_thick_layer_tab_CIE_RGB"
    cosmetic_dict["Clarins105"] = "Clarins105/Clarins_105_Nude_coll_norm_CIE_RGB"
    cosmetic_dict["Clarins108"] = "Clarins108/Clarinsfoundation108_5W_norm_thick_layer_tab_CIE_RGB"
    if backrefl_data:
        measured_data = os.path.join("scripts","validation", "measured_data","linkoping","CIE_RGB_back_reflection",cosmetic_dict[cosmetic_type] + ".npy")
    else:
        measured_data = os.path.join("scripts","validation", "measured_data","linkoping","CIE_RGB",cosmetic_dict[cosmetic_type] + ".npy")

    #print(measured_data)
    #input()

    return np.load(measured_data)
# def getMeasuredDataRGB(cosmetic_type):
# 	if cosmetic_type == "":
# 	    cosmetic_type = "Clarins112Thick"
#     cosmetic_dict = {}
    
# 	cosmetic_dict["Clarins112Thick"] = "Foundation_112_5W/Clarins112_5W_thick_layer_tab_RGB"
#     measured_data = os.path.join("scripts","validation", "measured_data","linkoping",cosmetic_dict[cosmetic_type] + ".npy")
    
#     return np.load(measured_data)


def getMeasuredData(helper,cosmetic_type=""): 
    
    # ---------------------Pick the cosmetic measurement file
    # I will save this in a json file and only reading it here
    if cosmetic_type == "":
        cosmetic_type = "Clarins112Thick"
    
    cosmetic_dict = {}
    cosmetic_dict["Clarins112Thick"] = "Foundation_112_5W/Clarins112_5W_thick_layer_tab"
    cosmetic_dict["Clarins112Thin"] = "Foundation_112_5W/Clarins112_5W_thin_layer_on_silicon_tab"
    cosmetic_dict["Clarins103N"] = "Foundation_103N/Clarinsfoundation103N_norm_thick_layer_tab"
    cosmetic_dict["Clarins108"] = "Foundation_108_5W/Clarinsfoundation108_5W_norm_thick_layer_tab"
    cosmetic_dict["Clarins105"] = "Foundation_105/Clarins_105_Nude_coll_norm"

    cosmetic_dict["Highlight_01"] = "Highlight_Rimmel/Highlight_coll_1_norm"
    cosmetic_dict["Highlight_02"] = "Highlight_Rimmel/Highlight_coll_2_norm"
    cosmetic_dict["Highlight_03"] = "Highlight_Rimmel/Highlight_coll_3_norm"    
    cosmetic_dict["Blusher_01"] = "Blush_Rimmel/Blusher_coll_1_norm"
    cosmetic_dict["Blusher_02"] = "Blush_Rimmel/Blusher_coll_2_norm"
    cosmetic_dict["Blusher_03"] = "Blush_Rimmel/Blusher_coll_3_norm"

    cosmetic_dict["FoamOnly"] = "FoamBlusher/Foam_coll_norm"
    cosmetic_dict["FoamBlusher"] = "FoamBlusher/FoamBlusher_coll_norm"

    # cosmetic_dict["Clarins108_5W"]
    # cosmetic_type = ["Highlight_Rimmel","FoundationClarins","BronzeRimmel","BlushRimmel","Foundation_103N","Foundation_108_5W","Foundation_112_5W"]
    
    cosmetic_table_name = ["Clarins112_5WOnsilicon1400nmLayer_norm","Clarins112_5W_thick_layer_tab","Clarins112_5W_thin_layer_on_silicon_tab",
            "Clarinsfoundation108_5W_norm_thick_layer_tab","Clarinsfoundation103N_norm_thick_layer_tab","RimmelBlushThinLayerOnBlack_norm","Blush_Bulk_coll_norm","HighlightRimmel_thick_layer_tab","Highlight_Bulk_coll_02"]
    
    #change this index for changing the type of material to 
    #cosmetic_type = "Foundation_108_5W"
    
    
    measured_data = os.path.join("scripts","validation", "measured_data","linkoping",cosmetic_dict[cosmetic_type] + ".txt")
    output_path = os.path.join("scripts","validation","measured_data","linkoping",cosmetic_dict[cosmetic_type]+ ".npy")
    print("Cosmetic measurements: {0}".format(measured_data))
        
    # measured_data = os.path.join("scripts","validation", "measured_data","linkoping",cosmetic_type,table_name + ".txt")
    # output_path = os.path.join("scripts","validation","measured_data","linkoping",cosmetic_type,table_name + ".npy")


    theta_angles = parseNormalizedReflectanceData(path=measured_data,output_path=output_path, helper = helper)
    #measure_data_rgb = os.path.join("scripts","validation","measured_data","linkoping",cosmetic_type[idx],cosmetic_table_name[idx] + "_RGB_tab.txt")
    #output_path_rgb= os.path.join("scripts","validation","measured_data","linkoping",cosmetic_type[idx],cosmetic_table_name[idx]+"_RGB_tab.npy")
    #parseRGBData(path=measure_data_rgb,output_path=output_path_rgb)

    # table_name = "Clarinsfoundation108_5W_norm_thick_layer_tab"

    data = np.load(output_path)

    # data_rgb =  np.load(output_path_rgb)
    out_folder = cosmetic_type + "_refl_plot"
    mkdir(out_folder)
    out_path = os.path.join(out_folder,"{0}.png".format(cosmetic_type))
    return [data,np.zeros(1),out_folder,out_path, theta_angles]

# Measure BRDF, the parameters are defined by a json file.
# It takes around 30 seconds per render
def measureBRDF(material_file = "./scripts/measurements/model_parameters.json", parameters = {}):

    #print("Calling measureBRDF")

    # command = "pbrt-v4/build/pbrt_test --test-filter Measurements.MeasureBRDF --nthreads 20"
    command = "pbrt-v4/build/pbrt_test --test-filter Measurements.EvalDiffuse --nthreads 20"

    # Save material parameters
    print("OPTIMIZIDE PARAMS DICT")
    print(parameters)
    #input()

    with open(material_file, 'w', encoding='utf-8') as f:
        json.dump(parameters, f, ensure_ascii=False, indent=4)

    # Execute command
    print("Executing command: {0}".format(command))
    os.system(command)

    #print("Succesfully measuring the BRDF")


def averagedSpectralReflectancePlot(data,helper,theta_i, fig,axs,output_file = "spectral.pdf", markers = 'x',linestyle= '-', plot_format = "pdf",  parameters = [], query = {}, color = (0.5,0,0)):
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"
    #plot first positive theta_i angles (i.e. angles in which the camera and the light are NOT on the same hemisphere.) 
    _, theta_i_idx, _, _, _ = helper.getIndices(0, theta_i, 0, 0, 0)
    
    colorbar = True
    
    # NB: probably it is better to create a box to write the parameters https://matplotlib.org/stable/tutorials/text/text_intro.html
    if parameters:
        #title    = "Spectral BRDF Plot: mean = {0}, thick = {1}, $\\sigma$ = {3}, {2}".format(parameters[0], parameters[2], parameters[3],parameters[4])
        title    = "Spectral BRDF Plot: $\\sigma$ = {0:0.3f}, $D \\%$ = {1:0.2f}, $g_1$ = {2:0.2f}, $g_2$ = {3:0.2f}, $w_g$ = {4:0.2f}, \n $\\mu_\\theta$ = {5:0.2f}, $\\alpha_d$ = {6: 0.2f}, $\\alpha_p$ = {7: 0.2f}, $E_D$ = {8:0.2f}, $E_R$ = {9:0.2f}".format(parameters["mean_size"], \
                 parameters["diff_particles_percs"], parameters["g1"], parameters["g2"], parameters["w_g"], parameters["means_or"], parameters["albedo_diff"], parameters["albedo_flake"], query["data_loss"], query["reg_loss"])
    else:
        title    = "Spectral BRDF Plot"

    #theta_o_indices = query["theta_angles"][theta_i]

    #print(theta_o_indices)
    #print(len(theta_o_indices))

    #first 90 are negative values
    neg_vals = data[:,theta_i_idx,0:90,:]
    pos_vals= data[:,theta_i_idx,91:,:]


    #does this make it sense?
    # avg_refl = np.mean(data,axis=3)
    
    
    # plt.ylim((0.0,1.0))

    if np.all( pos_vals[:] == 0):
        print("for theta_i_idx {0} all values are zero".format(theta_i))
        # print(avg_refl[0,theta_i_idx,0,:])
        return 
        positive_val_mask = avg_refl[0,theta_i_idx,phi_o_idx,:] > 0
    
    # positive_val_mask = avg_refl[0,theta_i_idx,0,:] > 0
    # print(positive_val_mask[0])
    # naffo = pos_vals[:,theta_i_idx,:,:]
    # naffo = naffo.flatten()
    # naffo = naffo[:helper.theta_o_bins.size] 
    # positive_val_mask = naffo[:] > 0
    # avg_refl_plt = naffo[positive_val_mask]

    
    # pos_vals  = pos_vals.flatten()
    # pos_vals = pos_vals[:helper.theta_o_bins.size]
    
    positive_val_mask = pos_vals[:] > 0
    avg_refl_plt = pos_vals[positive_val_mask]
    #print(naffo.shape)
    #print(naffo)


    #avg_refl_plt = naffo
    
    x = helper.theta_o_bins

    # x = helper.theta_o_bins[1:helper.theta_o_bins.size]
# 

    #x = x[positive_val_mask] -1 # There is 1 degree shift (1, 11, ..)
    # x = x[positive_val_mask]

    #print("theta_i = {0}, x = {1}".format(theta_i, x))
    #print("Size of average reflectance: {0}".format(avg_refl.shape))
    #print("Size of x(model), theta_i = {0}: {1}".format(theta_i, x.shape))
    #print("Size of average reflectances (model): {0}".format(avg_refl_plt.shape))

    #input()
    # pos_vals = np.nanmean(np.where(data > 0, data, np.nan), axis=1
    # pos_vals = np.nanmean(pos_vals,axis=2).flatten()

    # pos_vals = np.mean(neg_vals,axis=2).flatten()
    pos_vals = pos_vals[0,:,1].flatten()
    plt.plot(x, pos_vals,linewidth=1.0,color = color,alpha=0.5, linestyle=linestyle,label= r'$\theta_i$ = {0} - fit'.format(theta_i) )
    plt.scatter(x, pos_vals,s=5,c=color,alpha=0.5)
    # neg_vals = np.mean(neg_vals,axis=2).flatten()
    neg_vals = neg_vals[0,:,1].flatten()
    x = np.linspace(-90,0,90)
    plt.plot(x,np.flip(neg_vals),linewidth=1.0,color = color,alpha=0.5, linestyle=linestyle,label= r'$\theta_i$ = {0} - fit'.format(theta_i) )
    plt.scatter(x, np.flip(neg_vals),s=5,c=color,alpha=0.5)
    
    plt.xlabel(r'$\theta_o$', size=14, ha='center', va='top')
    l=plt.ylabel(r'$R$', size=14, ha='right', va='center')
    l.set_rotation(0)
    plt.legend()
    plt.title(title)

    #plt.savefig(output_file)

def plotAveragedSpectralReflectanceMeasurements(data,helper,theta_i,theta_o, fig,axs,output_file = "spectral.pdf",color = (0.5,0,0)):
    
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"
    _, theta_i_idx, _, _, _ = helper.getIndices(0, theta_i, 0, 0, 0)

    avg_refl = np.mean(data,axis=4) #Average spectral reflectance

    #print("Average reflectance: {0}".format(avg_refl.shape))

    #check that there are some data to plot
    if np.all( avg_refl[0,theta_i_idx,0,:] == 0):
        print("for theta_i_idx {0} all values are zero".format(theta_i))
        # print(avg_refl[0,theta_i_idx,0,:])
        return 
        positive_val_mask = avg_refl[0,theta_i_idx, :,:] > 0
    
    naffo = avg_refl[:,theta_i_idx,:,:]
    #naffo = avg_refl[:,theta_i_idx,:, :]
    #naffo = naffo[:, :, theta_o]
    naffo = naffo.flatten()
    # print(naffo)
    # naffo = naffo[1:]

    # Non-zero positions matches the measurement angles (seems to work fine)
    positive_val_mask = naffo[:] > 0
    avg_refl_plt = naffo[positive_val_mask]
    x = helper.theta_o_bins
    x = x[positive_val_mask] -1 # There is 1 degree shift (1, 11, ..)

    #print("theta_i = {0}, x = {1}".format(theta_i, x))
    #print("Size of x, theta_i = {0}: {1}".format(theta_i, x.shape))
    #print("Size of average reflectances: {0}".format(avg_refl_plt.shape))

    plt.plot(x, avg_refl_plt,linewidth=1.0,color = color,linestyle='-',label= r'$\theta_i$ = {0}'.format(theta_i) )
    plt.scatter(x, avg_refl_plt,s=5,c=color)
    
    plt.xlabel(r'$\theta_o$', size=14, ha='center', va='top')
    l=plt.ylabel(r'$R$', size=14, ha='right', va='center')
    l.set_rotation(0)
    plt.legend()

    ax = plt.gca()
    # ax.set_xlim([-10,75]) # We should update this with measurement data
    ax.set_xlim([-10,90]) # We should update this with measurement data


    #input()
    #plt.savefig(output_file)

def plotFitting(basename, output_folder, use_model, model_brdf, use_measurements, data_brdf, helper_bxdf, parameters, query, plot_avg_refl = True):

    # Parameters
    plot_format = "png"
    phi_i = 0
    phi_o = 0
    wavelength = 650
    theta_o = query["theta_o"]

    phi_i = 0
    iterator = 0
    # linestyles = ['-','--',':','-.','-',';',';','-']
    # markers = ['o','x','v','.','s']

    #parameters = ['0', '0', '0', '8', '0', '0.02', '0', '0.95', '1', '0.95', '0.5', '0', '0.5']

    fig, ax = plt.subplots(figsize=(8, 4))
    theta_i_angles = [30,40,50,60,70]
    # theta_i_angles = [30]
    
    if plot_avg_refl:
        for theta_i in theta_i_angles:
            #print(theta_i)

            phi_i_idx, theta_i_idx, phi_o_idx, theta_o_idx, wavelength_idx = helper_bxdf.getIndices(phi_i, theta_i, phi_o, theta_o, wavelength)
            
            if use_model:

                output_file = os.path.join(output_folder,"avg_reflectance_thetaI_{0}_thetaO_{1}.{2}".format(theta_i,theta_o,plot_format))
                
                # avg_refl = np.mean(model_brdf,axis=4)
                
                # #print(model_brdf.shape)
                # #print("Model shape: {0}".format(model_brdf.shape))

                # if np.all( avg_refl[:,theta_i_idx,:,:] == 0):
                #     print("theta i {0} with idx {1} is zero ".format(theta_i,theta_i_idx))
                #     # print(avg_refl[:,theta_i_idx,:,:])
                #     continue 
                t_color_list = ['b','g','r','c','m','y','k']
                linestyles = ['']

                """
                averagedSpectralReflectancePlot(model_brdf, helper=helper_bxdf,theta_i = theta_i,output_file= output_file,
                        plot_format= plot_format,parameters=  parameters,color = t_color_list[iterator],linestyle= '--',
                            axs=ax,
                            fig=fig)
                """
                
                averagedSpectralReflectancePlot(model_brdf, helper=helper_bxdf,theta_i = theta_i,output_file= output_file,
                        plot_format= plot_format,parameters=  parameters, query = query,  color = t_color_list[iterator],linestyle= '--',
                            axs=ax,
                            fig=fig)

            #check if we have backreflectance
            #rgb_color = np.clip(rgb_color,0,1)
            use_color = False
            if use_color:
                rgb_color = rgb_color[rgb_color > 0]
                rgb_color_reshaped = rgb_color.reshape(-1,3)
                
                color_list = [tuple(row) for row in rgb_color_reshaped]
                t_color_list = []
                for c in rgb_color_reshaped:
                    magnitude = np.linalg.norm(c)
                    rgb_color =  c / magnitude
                    t_color_list.append(tuple(rgb_color))
                # print(t_color_list)
            else:
                t_color_list = ['b','g','r','c','m','y','k']
            # rgb_color_avg = np.mean(rgb_color, axis=1)


            if use_measurements:
                output_file = os.path.join(output_folder,"avg_reflectance_thetaI_{0}_thetaO_{1}.{2}".format(theta_i,theta_o,plot_format))
                #print(output_file)
                
                """
                averagedSpectralReflectancePlot(data_brdf, helper=helper_bxdf,theta_i = theta_i,output_file= output_file,
                        plot_format= plot_format,parameters=  parameters,color = t_color_list[iterator],linestyle= '-',
                            axs=ax,
                            fig=fig)
                """

                #print("Measurements shape: {0}".format(data_brdf.shape))

                #helper_measured = MaterialHelper(phi_i_res = 0, theta_i_res = 181, phi_o_res = 1, theta_o_res = 181, wavelength_res = 0,theta_o_min=-90,theta_o_max=90)
                #helper_measured = MaterialHelper(phi_i_res = 0, theta_i_res = 18, phi_o_res = 0, theta_o_res = 180, wavelength_res = 0,theta_o_min=-90,theta_o_max=90)
                #plot RGB patches
                #plotAveragedSpectralReflectanceMeasurements(data_brdf,helper=helper_measured,fig=fig,axs=ax,output_file= output_file,theta_i=theta_i,theta_o=theta_o,color=t_color_list[iterator])
                

                            
            iterator+=1
        
        plt.ylim((0.0,25))
        plt.savefig(os.path.join(output_folder, basename+".png"))
        plt.ylim((0.0,2.5))        
        plt.savefig(os.path.join(output_folder, basename+"_limited.png"))
    
    else:
        print("To be implemented")

    # plt.savefig(out_path)
    #plt.show()
def plotMeasurements(data_bxdf,out_folder):
    img = data_bxdf
    img = np.clip(img, 0.0, 1.0)
    plt.imshow(img, aspect = "auto")
    plt.ticklabel_format(style='plain', axis='x', useOffset=False)
    # plt.xticks(np.arange(0, 80, 5))
    plt.xticks([])

    plt.xlabel(r'$\theta_o$', size=12, ha='center', va='top')
    
    #ax.xlabel(r'Target', size=12, ha='center', va='top')
    plt.yticks([])
    # ax.yticks([])
    plt.ylabel(r'$\theta_i$', size=12, ha='center')
    plt.savefig(out_folder)


def plotCurves(data_bxdf,data_meas,theta_i):
    indexes = {30:0,40:1,50:2,60:3,70:4}
    colors = ['blue','green','orange',"cyan","purple"]
    m_index = indexes[theta_i]
    x = np.linspace(0,80,80)

    plt.plot(x,data_bxdf[m_index,:,0],linestyle = '-', linewidth=1.0,color = colors[m_index],alpha=0.5,label= r'$\theta_i$ = {0} - fit'.format(theta_i) )
    plt.plot(x,data_meas[m_index,:,0],linestyle = '-.', linewidth=1.0,color = colors[m_index],alpha=0.5,label= r'$\theta_i$ = {0} - meas'.format(theta_i) )
    # plt.scatter(x, pos_vals,s=5,c=color,alpha=0.5)
    

def plotData(out_folder,data_values2,data_bxdf,loss_img,loss, parameters,query,out_file = ""):
    
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(12, 9))

    axes[0].set_title("RGB ({0})".format(query["cosmetic_type"]))
    axes[0].set_ylabel("Data", rotation=90, size='large')
    axes[1].set_ylabel("Model", rotation=90, size='large')
    axes[2].set_ylabel("Loss", rotation=90, size='large')

    reference = plotPatch(axes[0],data_values2[:,:,:])
    ours = plotPatch(axes[1],data_bxdf[:,:,:])
    error = plotPatch(axes[2],loss_img[:,:,:],isErrorData=True)

    """
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(12, 9))
    axes[0,0].set_title("RGB ({0})".format(query["cosmetic_type"]))
    axes[0,1].set_title("R")
    axes[0,2].set_title("G")
    axes[0,3].set_title("B")

    axes[0,0].set_ylabel("Data", rotation=90, size='large')
    axes[1,0].set_ylabel("Model", rotation=90, size='large')
    axes[2,0].set_ylabel("Loss", rotation=90, size='large')

    reference = plotPatch(axes[0,0],data_values2[:,:,:])

    plotPatch(axes[0,1],data_values2[:,:,0])
    plotPatch(axes[0,2],data_values2[:,:,1])
    plotPatch(axes[0,3],data_values2[:,:,2])
    ours = plotPatch(axes[1,0],data_bxdf[:,:,:])
    tmp = data_bxdf[:,:,0]
    model_values_LAB = np.mean(tmp)     
    plotPatch(axes[1,1],data_bxdf[:,:,0])
    plotPatch(axes[1,2],data_bxdf[:,:,1])
    plotPatch(axes[1,3],data_bxdf[:,:,2])

    error = plotPatch(axes[2,0],loss_img[:,:,:],isErrorData=True)
    plotPatch(axes[2,1],loss_img[:,:,0])
    plotPatch(axes[2,2],loss_img[:,:,1])
    plotPatch(axes[2,3],loss_img[:,:,2])

    """

    # Ablation studies
    #title  = "BRDF RGB Patch Plot: $\\sigma$ = {0:0.3f}, $D \\%$ = {1:0.2f}, $\\alpha_p$ = ({2: 0.2f}, {3: 0.2f}, {4: 0.2f}), $\\eta$ = {5:0.2f}, $\\rho$ = {6:0.2f}, $E_D$ = {7:0.2f}, $E_R$ = {8:0.2f}".format(parameters["mean_size"], \
    #             parameters["diff_particles_percs"], parameters["albedo_flake_r"], parameters["albedo_flake_g"], parameters["albedo_flake_b"], parameters["ior"], parameters["surf_roughs"],  loss, query["reg_loss"])

    
    title  = "BRDF RGB Patch Plot: $\\sigma$ = {0:0.3f}, $D \\%$ = {1:0.2f}, $g_1$ = {2:0.2f}, $g_2$ = {3:0.2f}, $w_g$ = {4:0.2f}, \n $\\mu_\\theta$ = {5:0.2f}, $\\alpha_d$ = ({6: 0.2f}, {7: 0.2f}, {8: 0.2f}) , $\\alpha_p$ = ({9: 0.2f}, {10: 0.2f}, {11: 0.2f}), $E_D$ = {12:0.2f}, $E_R$ = {13:0.2f}".format(parameters["mean_size"], \
                 parameters["diff_particles_percs"], parameters["g1"], parameters["g2"], parameters["w_g"], parameters["means_or"], parameters["albedo_diff_r"], parameters["albedo_diff_g"], parameters["albedo_diff_b"], parameters["albedo_flake_r"], parameters["albedo_flake_g"], parameters["albedo_flake_b"], loss, query["reg_loss"])

    #plt.title("test")
    plt.suptitle(title)
    #plt.suptitle("Total Loss {0}".format(loss))
    # plt.ticklabel_format(style='plain', axis='x', useOffset=False)
    # plt.xticks(np.arange(0, 80, 5))
    # plt.xlabel(r'$\theta_o$', size=12, ha='center', va='top')
    if out_file == "":
        out_file = os.path.join(out_folder,"log.png")
    
    plt.savefig(out_file)
    plt.clf()
    
    """
    thetas = [30,40,50,60,70]
    for theta_i in thetas:
        plotCurves(data_bxdf,data_meas=data_values2,theta_i=theta_i)
    
    plt.legend()
    plt.title(title)
    out_file_curve = out_file[:-4]
    out_file_curve += "_curve.png"
    


    plt.savefig(out_file_curve)
    plt.clf()

    fig, axes = plt.subplots(figsize=(6,6))
    
    
    left, bottom, width, height = 0.1, 0.6, 0.25, 0.25
    ax_inset = fig.add_axes([left, bottom, width, height])
    plotPatch(axes,data_bxdf[:,:,:])
    # im = plotPatch(ax_inset,loss_img[:,:,:],isErrorData=True,use_ticks=False)
    im = plotError(ax_inset,loss_img)
    plt.colorbar(im,ax=ax_inset)
    out_file_inset = out_file[:-4]
    out_file_inset += "_inset.svg"
    plt.savefig(out_file_inset)
    plt.clf()
    fig, axes = plt.subplots(figsize=(6,6))
    plotPatch(axes,data_values2[:,:,:])

    out_file_inset = out_file[:-4]
    out_file_inset += "_ref.svg"
    plt.savefig(out_file_inset)
    plt.clf()
    # plt.savefig(out_file)
    """
    
from matplotlib.colors import LinearSegmentedColormap
def plotError(ax,img,vmin = None,vmax = None):
    ax.set_xticks([])
    ax.set_yticks([])
    if vmin:
        im = ax.matshow(img[:,:,0], aspect = "auto",vmin=vmin,vmax = vmax)
        
    else:
        im = ax.matshow(img[:,:,0], aspect = "auto")
    ax.set_xticks([])
    ax.set_yticks([])
    return im
def plotPatch(ax,img,isErrorData = False,use_ticks = True):
    mean_val = np.mean(img)
    img = np.clip(img, 0.0, 1.0,dtype=np.float32)
    #desired_shape = (5, 80)
    desired_shape = (6, 120) # Back reflectance
    use_ticks = True

    #print(img.shape)
    #print(np.mean(img))
    if use_ticks:
        custom_y_ticks = [0, 1, 2, 3, 4, 5 ]
        custom_y_tick_labels = ['20', '30', '40', '50', '60','70']

        custom_x_ticks = [0, 40, 119]
        custom_x_tick_labels = ['-40', '0', '80']

        ax.set_xticks(custom_x_ticks)
        ax.set_xticklabels(custom_x_tick_labels)


        # Set the custom y-tick positions and labels
        ax.set_yticks(custom_y_ticks)
        ax.set_yticklabels(custom_y_tick_labels)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    if img.shape == desired_shape:
        #channels
        #print(np.mean(img))
        ax.set_title("Avg val {0:0.3f}".format(mean_val))
        im = ax.imshow(img, aspect = "auto",cmap="viridis")
        # plt.colorbar(im,ax=ax)
    else:
        if isErrorData:
            colors = [(1, 0, 0), (0, 0, 0), (0, 1, 0)]  # Red, White, Green
    
# Define the corresponding values for each color
# Here, we set the color white to represent zero values
            values = [0, 0.5, 1]

# Create the colormap
            cmap = LinearSegmentedColormap.from_list("RedToGreen", list(zip(values, colors)))
            
            # img = img.mean(axis=2)
            im = ax.imshow(img[:,:,0], aspect = "auto",cmap = "RdBu")
            # im = ax.matshow(img[:,:,0], aspect = "auto")

            # plt.show()
            return im
        else:
            ax.imshow(img, aspect = "auto")
    # ax.ticklabel_format(style='plain', axis='x', useOffset=False)
    # ax.xticks(np.arange(0, 80, 5))
    # ax.xlabel(r'$\theta_o$', size=12, ha='center', va='top')
    
    #ax.xlabel(r'Target', size=12, ha='center', va='top')

    #plt.yticks(theta_i)
    # ax.yticks([])
    # ax.ylabel(r'$\theta_i$', size=12, ha='center', va='center_baseline')
    return ax
def plotRGBPatch(img, output_file, name = {0}):

        img = np.clip(img, 0.0, 1.0)

        if name:
            title = "BRDF patch {0}".format(name)
        else:
            title = "BRDF patch"

        plt.title(title)

        plt.rcParams["figure.figsize"] = (10, 8)

        #fig = plt.gcf()
        #fig.set_size_inches(30, 20)

        plt.ticklabel_format(style='plain', axis='x', useOffset=False)
        plt.xticks(np.arange(0, 80, 5))
        plt.xlabel(r'$\theta_o$', size=12, ha='center', va='top')
        plt.imshow(img)

        #plt.yticks(theta_i)
        plt.yticks([])
        plt.ylabel(r'$\theta_i$', size=12, ha='center', va='center_baseline')

        plt.savefig(output_file)

        #matplotlib.image.imsave(output_file, img)
def avoid_data_extrapolation(matrix_values):
    if matrix_values.shape != (6,120,3):
        print("WARNING NOT THE EXPECTED SHAPE, ACTUAL SHAPE IS {0} ".format(matrix_values.shape))
    thetas = [20,30,40,50,60]
    #these are all hardcode with the data we have, it's just ad-hoc formula
    for idx,theta in enumerate(thetas):
        leftmost_known_angle = 30-theta
        offset = 40 + leftmost_known_angle
        matrix_values[idx,0:offset,:] = 0
    return matrix_values

def remove_data_extrapolation(model_values,data_values):
    model_values = np.where(data_values == 0, 0, model_values)
    return model_values    
    
# Objective function to fitting a model.
def objective(x, data_values, parameters, parameter_names, query):


    print("Trying solution: {0}".format(x))

    data_loss_format = query["data_loss_format"]
    reg_loss_format = query["reg_loss_format"]
    reg_w = query["reg_w"]
    R_w = query["R_w"]
    G_w = query["G_w"]
    B_w = query["B_w"]
    checkpoint_iters = query["checkpoint_iters"]
    
    color_space_data = data_loss_format.split("_")[0]
    metric_data = data_loss_format.split("_")[1]

    print("Data loss: Metric {0}, Color Space {1}".format(metric_data, color_space_data))
    #input()

    loss_subset = query["loss_subset"]

    # ---------------- Update material parameters with candidate solution
    for idx, parameter_name in enumerate(parameter_names): 
        parameters[parameter_name] = x[idx]
        # print("idx {0} params name {1} value {2:.2f}".format(idx,parameter_name,x[idx]))

    #print("Parameters = {0}".format(parameters))

    # Parameters helpful to debug the plot construction of the measured model data
    
    # Diffuse particles
    #parameters["mean_size"] = 0.05
    #parameters["diff_particles_percs"] = 0.97
    #parameters["means_or"] = 0.0 

    # Specular microflakes
    #parameters["mean_size"] = 0.02
    #parameters["diff_particles_percs"] = 0.0
    #parameters["means_or"] = 0.0 

    # Ablation studies (Task 2)
    #parameters["albedo_flake_r"] = 1.0
    #parameters["albedo_flake_g"] = 0.98
    #parameters["albedo_flake_b"] = 0.85
    #parameters["surf_roughs"] = 0.02
    #parameters["ior"] = 1.2
    #parameters["mean_size"] = 0.4
    #parameters["n_samples"] = [512]
    

    """
    # Testing parameters directly
    parameters = {
        "name": "clarins108",
        "albedo_flake_r": 0.99,
        "albedo_flake_g": 0.91,
        "albedo_flake_b": 0.99,
        "albedo_diff_r": 0.99,
        "albedo_diff_g": 0.94,
        "albedo_diff_b": 0.90,
        "surf_roughs" : 0.0,
        "thickness": 8.0,
        "refl_levels": 0,
        "means_or": 1.38,
        "stds_or": 0.0,
        "mean_size": 0.149,
        "stds_size": 0.0,
        "ior": 1.0,
        "diff_particles_percs":0.9,
        "g1": 0.18,
        "g2": -0.24,
        "w_g": 0.54,
        "n_samples": [128], 
        "distribution": "ggx"
    }
    """

    # --------------------------------Measure model BRDF
    # Subset dataset options: average spectral reflectance, filter by theta_i, filter smaller values

    material_file = "./scripts/measurements/model_parameters.json"
    measureBRDF(material_file = material_file, parameters = parameters)

    basename = "fitting_model"
 
    #helper_bxdf = MaterialHelper(phi_i_res = 0, theta_i_res = 180, phi_o_res = 0, theta_o_res = 180, wavelength_res = 0,theta_o_min=-90,theta_o_max=90) # Model sanity check!
    #work with 3 wavelengths
    # helper_bxdf = MaterialHelper(phi_i_res = 0, theta_i_res = 91, phi_o_res = 2, theta_o_res = 91, wavelength_res = 3,theta_o_min=0,theta_o_max=90) # Model sanity check!
    #helper_bxdf = MaterialHelper(phi_i_res = 0, theta_i_res = 5, phi_o_res = 2, theta_o_res = 91, wavelength_res = 3,theta_o_min=0,theta_o_max=90,theta_i_min=30,theta_i_max=70) # Model sanity check!
    
    helper_bxdf = MaterialHelper(phi_i_res = 0, theta_i_res = 6, phi_o_res = 2, theta_o_res = 91, wavelength_res = 3,theta_o_min=0,theta_o_max=90,theta_i_min=20,theta_i_max=70) # Model sanity check!
    #helper_bxdf = MaterialHelper(phi_i_res = 0, theta_i_res = 91, phi_o_res = 2, theta_o_res = 181, wavelength_res = 3,theta_o_min=-90,theta_o_max=90,theta_i_min=20,theta_i_max=70) # Backreflectance

    #helper_measured = MaterialHelper(phi_i_res = 0, theta_i_res = 91, phi_o_res = 0, theta_o_res = 181, wavelength_res = 1,theta_o_min=-90,theta_o_max=90)

    #(data_bxdf,data_rgb_bxdf,out_folder,out_path) = getBSDFData(basename, helper_bxdf, query["output_folder"])
    (data_bxdf,data_rgb_bxdf,out_folder,out_path) = getBSDFData(basename, helper_bxdf, query["output_folder"], "", "", True) # Include backreflectance data
    #print("model.shape {0} model.rgb.shape {1}".format(data_bxdf.shape,data_rgb_bxdf.shape))
    model_values = data_bxdf

    print("Model size: {0}".format(model_values.shape))


    #model_values = model_values[:, :120, :]# Back reflectance
    #old working code
    #model_values = model_values[:, 50:170, :]# Back reflectance
    model_values = model_values[:, 50:165, :]# Back reflectance

    #input()

    # NB: it might be helpful to optimize a subset of the data (only the plot?)
    # helper_measured = MaterialHelper(phi_i_res = 0, theta_i_res = 181, phi_o_res = 0, theta_o_res = 181, wavelength_res = 0,theta_o_min=-90,theta_o_max=90)
    #old code
    #helper_measured = MaterialHelper(phi_i_res = 0, theta_i_res = 181, phi_o_res = 0, theta_o_res = 181, wavelength_res = 0,theta_o_min=-90,theta_o_max=90)
    
    #Alina data goes from (ideally) -90 to 90 with a res of 181.
    # helper_measured = MaterialHelper(phi_i_res = 0, theta_i_res = 91, phi_o_res = 0, theta_o_res = 181, wavelength_res = 81,theta_o_min=-90,theta_o_max=90)
    helper_measured = MaterialHelper(phi_i_res = 0, theta_i_res = 91, phi_o_res = 0, theta_o_res = 181, wavelength_res = 1,theta_o_min=-90,theta_o_max=90)

    #print("Measurements size: {0}".format(data_values.shape))
    data_values2 = data_values
    # global NIters
    #----------------------------Data term-------------------------------------------
    NIters = query["NIters"] + 1
    query["NIters"] = NIters
    loss = 0.0

    tmp_data_loss = 0.0
    
    # Setting image format
    #print("Data size")
    #print(data_values.size)
    # breakpoint()
    if data_values.size == (2400):
        #print("baron")
        data_values = data_values[:, 2:].reshape([6, 80, 3])
        #ignore first row because we dont always have data for incident light 20°
        #data_values = data_values2[1:,:,:]
        data_values = data_values[1:,:,:]
    #else:
    #    data_values = data_values[:, 2:].reshape([6, 115, 3]) # Backreflectance
        #data_values = data_values[:, 2:].reshape([5, 80, 3])





    #print("Measurements size 2: {0}".format(data_values.shape))

    #print(data_values2.shape)
    #print(data_bxdf.shape)
    #print(data_values.shape)

    data_values2 = np.copy(data_values)
    #indices = data_values2 > 0.95
    #indices = np.where(data_values2 > 0.95)

    #print("Indices: {0}".format(indices))
    #model_values[:] = data_values[0, 0, :]

    # The data have always non-zero values on the whole 6x120 matrix, since we are interpolating. 
    # This is not correct as we dont have backrefl data for all lines
    # Thus we set to zero the regions in which we know we dont have data
    
    # model_values = avoid_data_extrapolation(model_values)
    # data_values = avoid_data_extrapolation(data_values)
    # data_values2 = avoid_data_extrapolation(data_values2)
    model_values = remove_data_extrapolation(model_values=model_values,data_values=data_values)


    if loss_subset == "diffuse":
        model_values = np.where(data_values2 > 0.95, 0, model_values)
        data_values = np.where(data_values2 > 0.95, 0, data_values2)
        
    # Focus on highlight --> Replace diffuse by 0s
    if loss_subset == "highlight":
        model_values = np.where(data_values2 < 0.95, 0, model_values)
        data_values = np.where(data_values2 < 0.95, 0, data_values2)
    if loss_subset == "cbrt":
        model_values = np.cbrt(model_values)
        data_values =  np.cbrt(data_values2)
    
    model_values_LAB = np.mean(model_values[:,:,0])     
    
    #model_values = data_values - 0.8

    #model_values[indices] = 0
    #data_values[indices] = 0

    if color_space_data == "RGB":
        
        tmp_data_loss = 0.0
        model_R  = model_values.mean(axis=(0,1))[0]
        model_G  = model_values.mean(axis=(0,1))[1]
        model_B  = model_values.mean(axis=(0,1))[2]
        # model_color = color.RGB(model_R,model_G,model_B)

        data_values_R  = data_values.mean(axis=(0,1))[0]
        data_values_G  = data_values.mean(axis=(0,1))[1]
        data_values_B  = data_values.mean(axis=(0,1))[2]
        # data_color = color.RGB(data_values_R,data_values_G,data_values_B)
        # colordiff = color.colordiff(color.RGB(model_color,data_color))
        data_values_LAB = color.rgb2lab(data_values)
        model_values_LAB = color.rgb2lab(model_values)
        # print(colordiff)
        if metric_data == "L1":
            tmp_data_loss = np.mean(np.abs(data_values - model_values))

        if metric_data == "L2":
            tmp_data_loss +=  np.sqrt(np.mean((data_values[:,:,0] - model_values[:,:,0])**2))
            tmp_data_loss +=  np.sqrt(np.mean((data_values[:,:,1] - model_values[:,:,1])**2))
            tmp_data_loss +=  np.sqrt(np.mean((data_values[:,:,2] - model_values[:,:,2])**2))

            # tmp_data_loss += R_w * np.sqrt(np.mean((data_values[:,:,0] - model_values[:,:,0])**2))
            # tmp_data_loss += G_w * np.sqrt(np.mean((data_values[:,:,1] - model_values[:,:,1])**2))
            # tmp_data_loss += B_w * np.sqrt(np.mean((data_values[:,:,2] - model_values[:,:,2])**2))



            #tmp_data_loss += np.sqrt(np.mean((measured_wl.flatten() - model_wl.flatten())**2))

        if metric_data == "L3":
            tmp_data_loss = np.sum(np.abs(model_values - data_values)**3)

        # Mean absolute logarithmic loss (might be helpful for scale fitting)
        if metric_data == "AbsLog":
            tmp_data_loss = np.mean(np.abs(np.log(1 + model_values) - np.log(1 + data_values)))

        if metric_data == "RMSE": 
            np.linalg.norm(model_values -  data_values) / np.sqrt(len(model_values))

        if metric_data == "SSIM":

            #print("Computing SSIM loss")
            #input()

            data_values_ssim = data_values
            model_values_ssim = model_values
            
            #data_values_ssim = resize(data_values, (5, 80))
            #model_values_ssim = resize(model_values, (5, 80))

            #tmp_data_loss += ssim(img, img+2.5)
            tmp_data_loss += ssim(data_values_ssim, model_values_ssim, win_size = 3)

            # Histogram loss

        if metric_data == "Histogram":

            #print("Computing Histogram loss")
            #input()

            data_values3 = data_values.reshape((-1, 3))
            model_values3 = model_values.reshape((-1, 3))

            hist_data, bin_edges = np.histogram(data_values, bins=10)
            hist_model, bin_edges = np.histogram(model_values, bins=10)

            tmp_data_loss = abs(hist_data.sum()-hist_model.sum())
            tmp_data_loss = np.linalg.norm(hist_data-hist_model)

            # Plot histograms
            
            # Random test
            #x = np.random.randn(1000, 3)
            #ax.hist(x, bins='auto', density=True, histtype='bar', stacked=True)

            fig, ax = plt.subplots(nrows=1, ncols=1)
            #ax.hist(data_values3, bins='auto', density=True, histtype='bar', stacked=True, color=["red", "green", "blue"])
            ax.hist(model_values3, bins=10, density=True, histtype='bar', stacked=True, color=["red", "green", "blue"])
            
            ax.set_title("Histogram model ({0})".format(parameters["name"]))

            #fig.savefig("./scripts/measurements/fitting/histogram_data.png")
            fig.savefig("./scripts/measurements/fitting/histogram_model.png")
        
    # LAB space loss --> Ignore perceived lightness

    if color_space_data == "LAB":

        #print("LAB space try")

        # data_values_LAB = color.rgb2lab(data_values.mean(axis=(0,1)))
        # model_values_LAB = color.rgb2lab(model_values.mean(axis=(0,1)))

        """
        # Check if color transformation is truncating the data
        data_values_RGB = color.lab2rgb(data_values_LAB)

        print("Data RGB: [{0}, {1}]".format(data_values.min(), data_values.max()))
        print("Data LAB: [{0}, {1}]".format(data_values_LAB.min(), data_values_LAB.max()))
        print("Data RGB2: [{0}, {1}]".format(data_values_RGB.min(), data_values_RGB.max()))

        input()
        """
        
        tmp_data_loss = 0.0
        
        if metric_data == "L2":
            data_values_LAB = color.rgb2lab(data_values)
            model_values_LAB = color.rgb2lab(model_values)
        
            tmp_data_loss += np.sqrt(np.mean((data_values_LAB[:,:,0] - model_values_LAB[:,:,0])**2))
            tmp_data_loss += np.sqrt(np.mean((data_values_LAB[:,:,1] - model_values_LAB[:,:,1])**2))
            tmp_data_loss += np.sqrt(np.mean((data_values_LAB[:,:,2] - model_values_LAB[:,:,2])**2))
        if metric_data == "colorDistance":
            data_values_LAB_R = np.mean(data_values[:,:,0])
            model_values_LAB_R = np.mean(model_values[:,:,0])         
            data_values_LAB_G =  np.mean(data_values[:,:,1])
            model_values_LAB_G = np.mean(model_values[:,:,1])         
            data_values_LAB_B =  np.mean(data_values[:,:,2])
            model_values_LAB_B = np.mean(model_values[:,:,2])         
            data_values_LAB = color.rgb2lab([data_values_LAB_R,data_values_LAB_G,data_values_LAB_B])
            model_values_LAB = color.rgb2lab([model_values_LAB_R,model_values_LAB_G,model_values_LAB_B])
            # model_values_LAB = color.rgb2lab(model_values)
            tmp_data_loss = color.deltaE_ciede2000(data_values_LAB,model_values_LAB)
            # print("IMPORTANT DATA")
            # print(tmp_data_loss)
            # print(data_values_LAB)
            # print(model_values_LAB)    
    data_loss = tmp_data_loss
    # loss_img = np.sqrt((data_values[:,:,:] - model_values[:,:,:])**2)
    loss_img = data_values[:,:,:] - model_values[:,:,:]
    # print("LOOOSSS IMAGEEE")
    # print(loss_img[:,:,0])
    if np.isnan(data_loss):
        data_loss = np.inf
    
    # ---------------------------Regulartization term-----------------------------
    # TO DO: add all metric from data term?
    #reg_loss = 1.0/np.sum(model_values)

    reg_loss = 0.0
    reg_loss += R_w * np.abs(model_values[:,:,0].sum() - data_values[:,:,0].sum())
    reg_loss += G_w * np.abs(model_values[:,:,1].sum() - data_values[:,:,1].sum())
    reg_loss += B_w * np.abs(model_values[:,:,2].sum() - data_values[:,:,2].sum())
    
    #reg_loss = 1.0/np.linalg.norm(model_values)
    
    loss = data_loss + reg_w * reg_loss
    print("Data loss = {0}".format(loss))
    #print("Regularization loss: {0}".format(reg_loss))
    print("Candidate x = {0}, error = {1}".format(x, loss))

    # ----------------------------Plot the solutions-------------------------------------------------------
    query["data_loss"] = data_loss
    query["reg_loss"] = reg_loss

    basename = "fitting_model"
    #plotData(query["output_folder"],data_values2,data_bxdf,loss_img,loss)
    model_values_LAB = np.mean(model_values[:,:,0])     
    
    plotData(query["output_folder"],data_values2 = data_values, data_bxdf = model_values,loss_img = loss_img,loss = loss, parameters= parameters,query=query)
    
    # plotRGBPatch(data_values2,"target.png")
    # plotRGBPatch(data_bxdf,"optim.png")
    # plotRGBPatch(loss_img,"loss.png")

    # Log only show model, easier to inspect the optimization behaviour
    # Later you can create a video using the following command inside the log folder: ffmpeg -framerate 1 -pattern_type glob -i '*.png' video.mp4 
    # If you want to only see the model predictions set False to data_values case!

    #basename = "fitting_model {:04d}".format(NIters)
    #log_folder = os.path.join(query["output_folder"], "log")
    #plotData(log_folder,data_values, model_values,loss_img,loss, parameters)

    #-----------------------------Checkpoint solutions-----------------------------------------------------
    # Save optimized solution in case the execution is terminated by the operating system
    
    #print("NIters = {0}".format(NIters))
    optimized_parameters = parameters

    if query["HQ"]:
        print("Saving HQ data")
        np.save(os.path.join(query["output_folder"],"best_result_model.npy".format(NIters)),model_values)
        np.save(os.path.join(query["output_folder"],"best_result_ref.npy".format(NIters)),data_values)


    if loss < query["min_loss"]:
        log_folder = os.path.join(query["output_folder"], "log")
        parameters_file = os.path.join(log_folder, "best_loss_{:06d}.json".format(NIters))
        query["min_loss"] = loss
        print("Saving checkpoint: {0}".format(parameters_file))
        #input()
        out_file = os.path.join(log_folder,"best_loss_{:06d}.svg".format(NIters))
        

        plotData(query["output_folder"],data_values, model_values,loss_img,loss, parameters,query,out_file)
        with open(parameters_file, 'w', encoding='utf-8') as f:
            json.dump(optimized_parameters, f, ensure_ascii=False, indent=4)
        np.save(os.path.join(log_folder,"best_loss_model_{:06d}.npy".format(NIters)),model_values)

    if NIters % checkpoint_iters == 0:
        log_folder = os.path.join(query["output_folder"], "log")
        parameters_file = os.path.join(log_folder, "optimized_parameters_{:06d}.json".format(NIters))
        print("Saving checkpoint: {0}".format(parameters_file))
        #input()
        #out_file = os.path.join(log_folder,"opt_params_{:06d}.png".format(NIters))
        #plotData(query["output_folder"],data_values, model_values,loss_img,loss, parameters,query,out_file)

        with open(parameters_file, 'w', encoding='utf-8') as f:
            json.dump(optimized_parameters, f, ensure_ascii=False, indent=4)
        
        numpy_folder  = os.path.join(log_folder,"numpy")
        mkdir(numpy_folder)        
        np.save(os.path.join(numpy_folder,"model_{:06d}.npy".format(NIters)),model_values)
        np.save(os.path.join(numpy_folder,"loss_{:06d}.npy".format(NIters)),data_loss)
    
    if query["HQ"]:
        #if HQ is set then we need only 1 image at really high spp
        exit()
    # exit()

    return loss

def callbackF(Xi):
    
    pass


    """
    time_elapsed = time.time() - start_t

    if time_elapsed > max_time:
        print("Time limit reached")
        return True
    """

# TO DO: Read the parameters from the command line like output folder, cosmetic type and so on
def main():
    parser = argparse.ArgumentParser(description="Fitting model to measurements data")
    parser.add_argument("-optimizer", default="CMA-ES", help="Valid options: CMA-ES, Nelder-Mead, Powell")
    parser.add_argument("-output_folder", default="./scripts/measurements/fitting", help="Path to the output folder")
    parser.add_argument("-cosmetic_type",default="Clarins112" , help="Use this field to plot a measured cosmetic. Valid values are: Clarins112|Clarins112Thin|Clarins103|Clarins108|Clarins105| \
                        Highlight_01|Highlight_02|Highlight_03|Blusher_01|Blusher_02|Blusher_03|")
    parser.add_argument("-initial_parameters", default="./scripts/measurements/fitting/initial_parameters.json", help="Initial parameter file")
    parser.add_argument("-parameter_bounds", default="./scripts/measurements/fitting/parameter_bounds.json", help="Parameter bounds file")
    parser.add_argument("-spp", type = int, help="Resolution used for obtaining samples")
    parser.add_argument("-mode", type = str,default="full", help="Task that should be done. This limits the parameters that are optimized.Legal values are full|task01|task02")
    parser.add_argument("-isHQ", action = "store_true", help="Set this flag to render once the scene (avoid to iterate).")


    args = parser.parse_args()
    
    full_optim = False
    ablation_study_t_1 = True
    ablation_study_t_3 = False
    output_folder = args.output_folder
    cosmetic_type = args.cosmetic_type
    optimizer_type = args.optimizer
    parameters_file = args.initial_parameters
    parameter_bounds_file = args.parameter_bounds
    spp = args.spp
    isHQ = args.isHQ

    start_t = time.time()
    mode = args.mode
    log_folder = os.path.join(output_folder, "log")
    mkdir(output_folder)
    mkdir(log_folder)


    # Material parameters
    parameters = {}
    #parameters_file = "./scripts/measurements/fitting/initial_parameters.json"
    with open(parameters_file, 'r', encoding='utf-8') as f:
        parameters = json.load(f)

    if spp:
        print("spp is defined {0}".format(spp))
        parameters["n_samples"] = [spp]
    
    if mode == "task01":
        # print("****************PASSED************************")
        # Diffusers only
        parameters["diff_particles_percs"] = 1.0
        parameters["w_g"] = 0.5
        parameters["ior"] = 1.0
        parameters["surf_roughs"] = 0.0
    elif mode == "task02":
        # Platelets only
        parameters["diff_particles_percs"] = 0.0
        parameters["w_g"] = 1.0
        parameters["ior"] = 1.0
        parameters["surf_roughs"] = 0.0
    elif mode == "task03":
        #1 HG lobe
        parameters["w_g"] = 1.0
        parameters["ior"] = 1.0
        parameters["surf_roughs"] = 0.0
    elif mode == "task04":
        #Use Full-Model + Surf Interface
        parameters["ior"] = 1.3
        parameters["surf_roughs"] = 0.1

    parameters["stds_or"] = 0.0
    parameters["stds_size"] = 0.0
    parameter_bounds = {}
    #parameter_bounds_file = "./scripts/measurements/fitting/parameter_bounds.json"
    with open(parameter_bounds_file, 'r', encoding='utf-8') as f:
        parameter_bounds = json.load(f)


    # BRDF query, especially helpful for fitting configuration and plots
    # Color spaces: RGB, LAB
    # Metrics: L1, L2, L3, AbsLog, SSIM, Histogram, ...
    # Data loss options: <color_space>_<metric>
    # Reg loss options: <color_space>_<metric>
    # Loss_subset options: full, highlight, diffuse
    # Average method (R, G, B have the same weights)
    # Luminosity methods (R = 0.3, G = 0.59, B = 0.11)
    # "loss_subset": "cbrt", 
    query = {
        "phi_i_res": 0,
        "theta_i_res": 101,
        "phi_o_res" : 2,
        "theta_o_res": 101,
        "wavelength_res": 0,
        "output_folder": output_folder,
        "theta_i" : [],
        "theta_o" : [],
        "loss_subset": "", 
        "data_loss_format": "RGB_L1",
        "reg_loss_format": "colorDistance",
        "reg_w": 0.0,
        "R_w" : 0.3,
        "G_w" : 0.59,
        "B_w" : 0.11,
        "checkpoint_iters": 1,
        "NIters": 0,
        "HQ": False,
        "cosmetic_type":"cosmetic_type",
        "min_loss": np.inf
    }
    if isHQ:
        query["HQ"] = True    
    #"loss_subset": "diffuse", 
    #"data_loss_format": "LAB_colorDistance",

    # ----------------------------Reading data from measurements

    #helper_measured = MaterialHelper(phi_i_res = 0, theta_i_res = 181, phi_o_res = 1, theta_o_res = 181, wavelength_res = 0,theta_o_min=-90,theta_o_max=90)
    helper_measured = MaterialHelper(phi_i_res = 0, theta_i_res = 91, phi_o_res = 0, theta_o_res = 181, wavelength_res = 1,theta_o_min=-90,theta_o_max=90)
    # cosmetic_type = "Clarins112Thick"
    #oldy code
    # (data_measured,data_rgb_measured,out_folder,out_path, theta_angles) = getMeasuredData(helper_measured,cosmetic_type)
    #for rgb patches
    data_measured = getMeasuredDataRGB(cosmetic_type, True)
    print(data_measured.shape)
    print(data_measured.size)
    theta_angles = {}
    reflectance_data = data_measured

    print("Measured data {0}".format(data_measured.shape))

    #print(np.unique(data_measured[:, 1]))
    #input()

    query["theta_angles"] = theta_angles

    angles_file = "./scripts/measurements/angles.json"
    with open(angles_file, 'w', encoding='utf-8') as f:
        json.dump(theta_angles, f, ensure_ascii=False, indent=4)

    #theta_o_min = np.array(theta_o_values).min()
    #theta_o_max = np.array(theta_o_values).max()


    #angles_file = "./scripts/measurements/angles.txt"

    #saveDirections(angles_file, theta_i_values, theta_o_values, [], [])

        #helper_bxdf = MaterialHelper(phi_i_res = 0, theta_i_res = 101, phi_o_res = 2, theta_o_res = 101, wavelength_res = 0)
    # helper_bxdf = MaterialHelper(phi_i_res = 0, theta_i_res = 180, phi_o_res = 0, theta_o_res = 180, wavelength_res = 0,theta_o_min=-90,theta_o_max=90)

    # (data_bxdf,data_rgb_bxdf,out_folder,out_path) = getBSDFData(basename, helper_bxdf,output_folder)
    #     #print("data.shape {0} data_rgb.shape {1}".format(data_bxdf.shape,data_rgb_bxdf.shape))

    # reflectance_data = data_bxdf

    # basename = "solution2"
    # plotFitting(basename, output_folder, False, [], True, reflectance_data, helper_bxdf, parameters, query, plot_avg_refl = True)
    #     #plotFitting(basename, output_folder, True, reflectance_data, False, [], helper_bxdf, parameters, query, plot_avg_refl = True)

        

    #     # ----------------------------Reading data from measurements

    #     #helper_measured = MaterialHelper(phi_i_res = 0, theta_i_res = 181, phi_o_res = 1, theta_o_res = 181, wavelength_res = 0,theta_o_min=-90,theta_o_max=90)
    #     helper_measured = MaterialHelper(phi_i_res = 0, theta_i_res = 91, phi_o_res = 0, theta_o_res = 181, wavelength_res = 1,theta_o_min=-90,theta_o_max=90)
    #     # cosmetic_type = "Clarins112Thick"
    #     #oldy code
    #     # (data_measured,data_rgb_measured,out_folder,out_path, theta_angles) = getMeasuredData(helper_measured,cosmetic_type)
    #     #for rgb patches
    #     data_measured = getMeasuredDataRGB(cosmetic_type)
    #     theta_angles = {}
    #     reflectance_data = data_measured

    #     query["theta_angles"] = theta_angles

    #     angles_file = "./scripts/measurements/angles.json"
    #     with open(angles_file, 'w', encoding='utf-8') as f:
    #         json.dump(theta_angles, f, ensure_ascii=False, indent=4)



    #     #theta_o_min = np.array(theta_o_values).min()
    #     #theta_o_max = np.array(theta_o_values).max()


    #     #angles_file = "./scripts/measurements/angles.txt"

    #     #saveDirections(angles_file, theta_i_values, theta_o_values, [], [])

    #     plot_meas = False
    #     if plot_meas:
    #         out_folder = "./rgb_{0}.png".format(cosmetic_type) 
    #         data_values = reflectance_data
    #         if data_values.size == (2400):
    #             data_values2 = data_values[:, 2:].reshape([6, 80, 3])
    #             #ignore first row because we dont always have data for incident light 20°
    #             data_values2 = data_values2[1:,:,:]
    #         else:
    #             data_values2 = data_values[:, 2:].reshape([5, 80, 3])
    #         plotMeasurements(data_values2,out_folder)
    #     # exit()
    #     basename = "solution"
    #     # plotRGBPatch(data_values2,"target.png")
    #     # plotRGBPatch(data_bxdf,"optim.png")

    #     # plotFitting(basename, output_folder, False, [], True, reflectance_data, helper_measured, parameters, query, plot_avg_refl = True)

    #     print("Reference plot is ready")

    # ----------------------------Brute Force Optimization----------------------------

    if optimizer_type == "BruteForce":
        # TO DO: extend this strategy to the rest of the parameters? 
        # TO DO: implement low-discrepancy quasi-random (Halton) sampling

        print("Trying the optimizer = {0}".format(optimizer_type))

        tolerance = 1e-2
        best_loss = int(1e8)
        best_solution = np.array([0.0, 0.0, 0.0])

        #parameter_names = ["mean_size", "diff_particles_percs", "g1", "g2", "w_g", "means_or", "albedo_flake_r","albedo_flake_g","albedo_flake_b",\
        #                 "albedo_diff_r","albedo_diff_g","albedo_diff_b"] # Full optimization
        parameter_names = ["albedo_diff_r","albedo_diff_g","albedo_diff_b"]

        albedo_diff_r = np.linspace(0.0, 1.0, 100)
        albedo_diff_g = np.linspace(0.0, 1.0, 100)
        albedo_diff_b = np.linspace(0.0, 1.0, 100)

        parameter_values = np.array([x for x in itertools.product(albedo_diff_r, albedo_diff_g, albedo_diff_b)])

        for idx, value in enumerate(parameter_values):
            print("Trying parameters: {0}".format(value))

            loss = objective(value, reflectance_data, parameters, parameter_names, query)

            if loss < best_loss:

                best_loss = loss
                best_solution = value

            if loss < tolerance:
                break
                
        # Save optimized solution
        for idx, parameter_name in enumerate(parameter_names): 
            optimized_parameters[parameter_name] = x_opt[idx]

        parameters_file = os.path.join(output_folder, "optimized_parameters.json")

        with open(parameters_file, 'w', encoding='utf-8') as f:
            json.dump(optimized_parameters, f, ensure_ascii=False, indent=4)

        print("Succesfully fitting to reflectance data")

    #-----------------------------CMA-ES Optimization---------------------------------

    if optimizer_type == "CMA-ES":
        
        parameter_names = ["mean_size", "diff_particles_percs", "g1", "g2", "w_g", "means_or", "albedo_flake_r","albedo_flake_g","albedo_flake_b","albedo_diff_r","albedo_diff_g","albedo_diff_b"] # Full optimization
        #parameter_names = ["albedo_diff_r","albedo_diff_g","albedo_diff_b"] # Only albedo
        
        # Ablation studies
        parameter_names = ["mean_size", "ior", "surf_roughs", "albedo_flake_r","albedo_flake_g","albedo_flake_b"]
        #parameter_names = ["surf_roughs"]

        x_init = []
        min_bounds = []
        max_bounds = []

        for parameter_name in parameter_names:
            x_init.append(parameters[parameter_name])

            val_range = parameter_bounds[parameter_name]
            min_bounds.append(val_range[0])
            max_bounds.append(val_range[1])


        bounds = [min_bounds, max_bounds]
        #bounds = [ [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]

        print("Parameter initialization: {0}".format(x_init))
        print("Parameter bounds: {0}".format(bounds))
        print("USING CMA-ES")
        options = {}
        options["bounds"] = bounds
        options["ftarget"] = 1e-4
        #options['maxiter'] = 10 #maximum number of iterations (10 --> 5 min, 100 --> 40 min)
        #options['popsize'] = 10 #number of new solution per iteration


        optimized_parameters = parameters

        x_opt, es = cma.fmin2(objective_function=objective, x0=x_init, sigma0=0.5, 
                                args=(reflectance_data, parameters, parameter_names, query), 
                                options = options, eval_initial_x = True)

        # Save optimized solution
        for idx, parameter_name in enumerate(parameter_names): 
            optimized_parameters[parameter_name] = x_opt[idx]

        parameters_file = os.path.join(output_folder, "optimized_parameters.json")
        
        with open(parameters_file, 'w', encoding='utf-8') as f:
            json.dump(optimized_parameters, f, ensure_ascii=False, indent=4)

        print("Succesfully fitting to reflectance data")
    scipy_optimizers = ["Nelder-Mead", "Powell"]

    if optimizer_type in scipy_optimizers:

        print("Trying the optimizer = {0}".format(optimizer_type))

        tolerance = 1e-2 # Fine detail: 1e-6
        tolerance = 1e-4 # Fine detail: 1e-6
        maxIters = 1000000
        if mode == "full":
            parameter_names = ["mean_size", "diff_particles_percs", "g1", "g2", "w_g", "means_or", "albedo_flake_r","albedo_flake_g","albedo_flake_b",\
                        "albedo_diff_r","albedo_diff_g","albedo_diff_b","stds_or"] # Full optimization
        
        elif mode == "task01":
        #ablation study - Task 1
            print("************ WARNING IN ABLATION STUDY TASK 1 MODE!! ************")
            #diffusers only
            parameter_names = ["g1","g2","w_g", "albedo_diff_r","albedo_diff_g","albedo_diff_b"] # Full optimization
            # parameter_names = ["ior", "surf_roughs", "g1","albedo_diff_r","albedo_diff_g","albedo_diff_b"] # Full optimization
        elif mode == "task02":
            print("************ WARNING IN ABLATION STUDY TASK 2 MODE!! ************")            
            parameter_names = ["mean_size", "means_or", "albedo_flake_r","albedo_flake_g","albedo_flake_b"] # Full optimization
            # parameter_names = ["mean_size", "diff_particles_percs", "g1", "means_or", "albedo_flake_r","albedo_flake_g","albedo_flake_b","albedo_diff_r","albedo_diff_g","albedo_diff_b","stds_or"] # Full optimization
        #parameter_names = ["albedo_diff_r","albedo_diff_g","albedo_diff_b"] # Only albedo

        elif mode == "task03":
            print("************ WARNING IN ABLATION STUDY TASK 3 MODE!! ************")            
            parameter_names = ["mean_size", "means_or", "albedo_flake_r","albedo_flake_g","albedo_flake_b","g1","albedo_diff_r","albedo_diff_g","albedo_diff_b"] # Full optimization
        elif mode == "task04":
            print("************ WARNING IN ABLATION STUDY TASK 4 MODE!! ************")            
            parameter_names = ["ior","surf_roughs", "mean_size", "means_or", "albedo_flake_r","albedo_flake_g","albedo_flake_b","g1","albedo_diff_r","albedo_diff_g","albedo_diff_b"] # Full optimization
        elif mode == "color":
            print("************ WARNING IN COLOR MATCHING MODE!! ************")            
            parameter_names = ["albedo_flake_r","albedo_flake_g","albedo_flake_b","albedo_diff_r","albedo_diff_g","albedo_diff_b"] # Full optimization
        elif mode == "color02":
            print("************ WARNING IN COLOR MATCHING MODE!! ************")            
            parameter_names = ["albedo_diff_r","albedo_diff_g","albedo_diff_b"] # Full optimization
        
        x_init = []
        bounds = []
        for parameter_name in parameter_names:
            x_init.append(parameters[parameter_name])

            val_range = parameter_bounds[parameter_name]
            bounds.append((val_range[0], val_range[1]))

        print("Parameter initialization: {0}".format(x_init))
        print("Parameter bounds: {0}".format(bounds))

        optimized_parameters = parameters

        optimizer = optimizer_type # Options: "Nelder-Mead", "Powell"

        res = minimize(objective,
                    x_init,
                    method=optimizer,
                    args=(reflectance_data, parameters, parameter_names, query),
                    tol=tolerance,
                    options={"maxiter": maxIters, 'disp': True,"maxfev":maxIters},
                    bounds = bounds,
                    callback=callbackF
                    #callback = MinimizeStopper(1E-3)
                    )

        # Save optimized solution
        for idx, parameter_name in enumerate(parameter_names): 
            optimized_parameters[parameter_name] = res.x[idx]

        parameters_file = os.path.join(output_folder, "optimized_parameters.json")

        with open(parameters_file, 'w', encoding='utf-8') as f:
            json.dump(optimized_parameters, f, ensure_ascii=False, indent=4)


        print("Solution {0}".format(res))


if __name__ == "__main__":
    main()

