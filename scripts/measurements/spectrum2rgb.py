# Helper to transform Spectral Reflectance to RGB reflectance
# Example Spectrum --> sRGB: python ./scripts/measurements/spectrum2rgb.py -input_file ./scripts/measurements/Blush_Bulk_coll_norm.txt -output_file ./scripts/measurements/Blush_Bulk_coll_RGB.txt -illuminant D65
# Example Spectrum --> CIE RGB: python ./scripts/measurements/spectrum2rgb.py -input_file ./scripts/measurements/Blush_Bulk_coll_norm.txt -output_file ./scripts/measurements/Blush_Bulk_coll_RGB.txt -colour_space CIE\ RGB
# Example CIE RGB -> sRGB: python ./scripts/measurements/spectrum2rgb.py -input_file ./scripts/measurements/best_loss_model_000507.npy -output_file ./scripts/measurements/best_loss_model_000507_sRGB.npy -illuminant D65
# Tutorial: https://colour.readthedocs.io/en/v0.3.11/tutorial.html#convert-to-tristimulus-values

import numpy as np
import argparse

np.set_printoptions(suppress=True)

import colour
import colour.plotting
#import colour.plotting.colorimetry

import scipy.misc
import matplotlib.image
import matplotlib.pyplot as plt

#print(colour.__version__)
#from colour.plotting import *
#colour.plotting.colorimetry.colour_plotting_defaults()
#visible_spectrum_plot()

class SpectrumHelper:

    def __init__(self, lambda_min = 380, lambda_max = 780, lambda_res = 81, transformation = "RGB", illuminant = "D65", colour_space = "CIE RGB"):

        self.wavelengths = np.linspace(lambda_min, lambda_max, num = lambda_res)
        self.illuminant = illuminant
        self.transformation = transformation
        self.colour_space = colour_space
        self.use_back_reflection = True

    def readData(self, input_file = "./scripts/measurements/Clarins112_5W_thin_layer_on_silicon_norm_tab.txt"):

        # Format: theta_i theta_o spectrum reflectance (81 wavelengths)
        
        if input_file.endswith('.npy'):
            self.reflectance_data = np.load(input_file)
        else:
            self.reflectance_data = np.loadtxt(input_file)

    # NB: It should be a properp data binning indices
    def getIndices(self, theta_i):

        theta_i_values = [20, 30, 40, 50, 60, 70]

        for idx, value in enumerate(theta_i_values):

            if value == theta_i:
                return idx
       

    def plotRGBPatch(self, img, output_file, name = "", min_theta_o = 0, max_theta_o = 80):

        img = np.clip(img, 0.0, 1.0)

        if name:
            title = "BRDF patch {0}".format(name)
        else:
            title = "BRDF patch"

        #plt.title(title)

        plt.rcParams["figure.figsize"] = (10, 8)

        #fig = plt.gcf()
        #fig.set_size_inches(30, 20)

        #print(min_theta_o)
        #print(max_theta_o)
        #print(np.arange(min_theta_o, max_theta_o, 5))
        #input()
        max_theta_o = 75
        min_theta_o = -40
        x_labels = np.arange(int(min_theta_o), int(max_theta_o), 20)
        print(x_labels)
        #x_labels = list(map(str, x_labels))

        #print(x_labels) 
        #input()

        # Original code
        #plt.ticklabel_format(style='plain', axis='x', useOffset=False)
        #plt.xticks(np.arange(0, 80, 5))
        #plt.xticks(np.arange(min_theta_o, max_theta_o, 10))
        #plt.yticks([])
        #plt.ylabel(r'$\theta_i$', size=12, ha='center', va='center_baseline')
        #plt.imshow(img)

        # Backlighting code
        # code for plotting the correct xticks and yticks
        plt.xticks(np.arange(0, max_theta_o + abs(min_theta_o), 10), x_labels,size=24)
        plt.xlabel(r'$\theta_o$', size=34, ha='center', va='top')
        plt.yticks([0, 1, 2, 3, 4, 5], ["20", "30", "40", "50", "60", "70"],size=24)
        plt.ylabel(r'$\theta_i$', size=34, ha='center', va='center_baseline', rotation=90)
        # plt.imshow(img)
        
        #remove xticks and yticks for final results
        # plt.gca().set_xticks([])

        # plt.xticks([])
        # plt.yticks([])
        plt.imshow(img, aspect="auto")

        #plt.yticks(theta_i)

        plt.savefig(output_file )
        

        #matplotlib.image.imsave(output_file, img)

    def toRGB(self, output_file, interpolate = True, name = ""):

        #print(output_file)
        #input()

        n, m = self.reflectance_data.shape



        # Format: theta_i theta_o RGB reflectance
        reflectance_rgb = np.zeros([n, 5], dtype="float32")

        theta_i = np.unique(self.reflectance_data[:, 0])
        initial_theta_i_len = theta_i.shape[0]
        if theta_i.shape[0] == 5:
            theta_i = np.insert(theta_i,0,20)
            
        
        min_theta_o = np.min(self.reflectance_data[:, 1])
        max_theta_o = np.max(self.reflectance_data[:, 1])

        if not self.use_back_reflection:
            min_theta_o = 0
            max_theta_o = 80
        else:
            #fixed number of cols
            min_theta_o = -40
            max_theta_o = 75

        theta_o = np.arange(min_theta_o, max_theta_o, 1)

        #print(theta_o)
        #input()

        #theta_o = np.arange(min_theta_o, max_theta_o, 1)

        m_theta_o_n = max_theta_o - min_theta_o 

        # new_n = theta_i.shape[0] * theta_o.shape[0]

        #fix all matrices to be have the same rows and columns        
        new_n = 6 * ( max_theta_o - min_theta_o)

        #print(theta_i.shape)
        #print(theta_o.shape)
        #print(new_n)
        #input()

        reflectance_rgb = np.zeros([n, 5], dtype="float32")
        final_reflectance = np.zeros([new_n, 5], dtype="float32")

        #print("theta_o = [{0}, {1}]".format(min_theta_o, max_theta_o))
        #print("theta_o = {0}".format(theta_o))

        angles = {}

        #print("theta_i = {0}".format(theta_i))

        #print(self.reflectance_data[14, :])
        #print(theta_i.shape)
        #input()

        for i in range(n):

            spd_data = self.reflectance_data[i, 2:]
            sample_spd_data = {}

            for idx, wavelength in enumerate(self.wavelengths):
                # print(wavelength)

                key = int(wavelength)

                sample_spd_data[wavelength] = spd_data[idx]

            spd = colour.SpectralDistribution(sample_spd_data, name='Sample')
            #print(repr(spd))

            # Plotting the sample spectral power distribution.
            #colour.plotting.single_spd_plot(spd)

            # ----------------Convert to Tristimulus Values
            cmfs = colour.MSDS_CMFS['CIE 1931 2 Degree Standard Observer']
            illuminant = colour.SDS_ILLUMINANTS[self.illuminant]

            # Calculating the sample spectral power distribution *CIE XYZ* tristimulus values.
            XYZ = colour.sd_to_XYZ(spd, cmfs, illuminant)
            #print(XYZ)

            # The output domain of *colour.spectral_to_XYZ* is [0, 100] and the input
            # domain of *colour.XYZ_to_sRGB* is [0, 1]. We need to take it in account and
            # rescale the input *CIE XYZ* colourspace matrix.
            
            if self.transformation == "SpectrumtosRGB":
                RGB = colour.XYZ_to_sRGB(XYZ / 100)
            else:
                #NB: If we use the CIE RGB space, then we do not need the illuminant
                #NB: The last parameter can be set up to "Bradford"

                #XYZ = np.array([0.21638819, 0.12570000, 0.03847493])
                #illuminant = np.array([0.34570, 0.35850])
                self.colour_space = "CIE RGB"
                
                #print(self.colour_space)
                #input()
                RGB = colour.XYZ_to_RGB(XYZ / 100, self.colour_space, None, None)

            #print("XYZ = {0}, RGB = {1}".format(XYZ.shape, RGB.shape))

            reflectance_rgb[i, 0] = self.reflectance_data[i, 0]
            reflectance_rgb[i, 1] = self.reflectance_data[i, 1]
            reflectance_rgb[i, 2] = RGB[0]
            reflectance_rgb[i, 3] = RGB[1]
            reflectance_rgb[i, 4] = RGB[2]

            if not interpolate:

                theta_o_n = theta_o.shape[0]
                n = theta_i.shape[0] * theta_o.shape[0]

                #new_idx = self.getIndices(self.reflectance_data[i, 0]) * theta_o_n + self.reflectance_data[i, 1] - min_theta_o - 1
                new_idx = self.getIndices(self.reflectance_data[i, 0]) * theta_o_n + int(self.reflectance_data[i, 1]) - int(min_theta_o)
                
                # Ugly fix: first column and last column were wrong with the computation above
                if new_idx % theta_o_n == 0 and self.reflectance_data[i, 1] != min_theta_o:
                    new_idx -= 1

                """
                print(int(min_theta_o))

                if (self.reflectance_data[i, 1] == -40.0):
                    print("test")
                    print(new_idx)
                    print(self.reflectance_data[i, 0])
                    print(self.reflectance_data[i, 1])
                    print(self.getIndices(self.reflectance_data[i, 0]))
                    #new_idx += 1
                    input()
                """

                #if new_idx == 720:
                #    continue

                #print(theta_o_n)
                #print(min_theta_o)

                final_reflectance[new_idx, 0] = self.reflectance_data[i, 0]
                final_reflectance[new_idx, 1] = self.reflectance_data[i, 1]
                final_reflectance[new_idx, 2] = RGB[0]
                final_reflectance[new_idx, 3] = RGB[1]
                final_reflectance[new_idx, 4] = RGB[2]

        #np.savetxt(output_file, reflectance_rgb, delimiter='\t', fmt='%f')

        # Only interpolate where they are values
        # NB: Maybe better to interpolate directly for final reflectance (sparse measurements)
        if interpolate:
            theta_o_n = theta_o.shape[0]
            #theta_o_n = m_theta_o_n
            for idx, angle in enumerate(theta_i):
                
                aux = reflectance_rgb[reflectance_rgb[:, 0] == angle]
                if idx == 0 and initial_theta_i_len == 5:
                    continue
                # if aux[:,1].size == 0:
                    # continue
                
                new_theta_o_min = int(np.min(aux[:, 1]))
                new_theta_o_max = int(np.max(aux[:, 1]))
                
                new_theta_o_max = int(min(new_theta_o_max,75))
                # new_theta_o_min = int(min(new_theta_o_max,-40))
                new_theta_o = np.arange(new_theta_o_min, new_theta_o_max)

                #print(aux[:, 1])
                #print(new_theta_o_n)
                #input()

                # Old code
                """
                R = np.interp(theta_o, aux[:, 1], aux[:, 2])
                G = np.interp(theta_o, aux[:, 1], aux[:, 3])
                B = np.interp(theta_o, aux[:, 1], aux[:, 4])

                final_reflectance[idx * theta_o_n: (idx+1)*theta_o_n, 0] = angle
                final_reflectance[idx * theta_o_n: (idx+1)*theta_o_n, 1] = theta_o
                final_reflectance[idx * theta_o_n: (idx+1)*theta_o_n, 2] = R
                final_reflectance[idx * theta_o_n: (idx+1)*theta_o_n, 3] = G
                final_reflectance[idx * theta_o_n: (idx+1)*theta_o_n, 4] = B
                """

                # New code --> Most back reflectance entries are black (upper back diagonal)
                # NB: We should exlude theta_o angles outside of the data range [theta_o_min, theta_o_max] for a given theta_i
                R = np.interp(new_theta_o, aux[:, 1], aux[:, 2])
                G = np.interp(new_theta_o, aux[:, 1], aux[:, 3])
                B = np.interp(new_theta_o, aux[:, 1], aux[:, 4])
                #print("[{0}, {1}]".format(new_theta_o_min, new_theta_o_max))
                #print(R.shape)
                #print(G.shape)
                #print(B.shape)

                theta_o_min_idx = int(new_theta_o_min) - int(min_theta_o)
                theta_o_max_idx = int(new_theta_o_max) - int(min_theta_o)

                #print("[{0}, {1}]".format(theta_o_min_idx, theta_o_max_idx))
                #input()

                #final_reflectance[idx * new_theta_o_min: idx * new_theta_o_max, 0] = angle
                #final_reflectance[idx * new_theta_o_min: idx * new_theta_o_max, 1] = theta_o


                final_reflectance[idx * theta_o_n + theta_o_min_idx: idx * theta_o_n + theta_o_max_idx, 2] = R
                final_reflectance[idx * theta_o_n + theta_o_min_idx: idx * theta_o_n + theta_o_max_idx, 3] = G
                final_reflectance[idx * theta_o_n + theta_o_min_idx: idx * theta_o_n + theta_o_max_idx, 4] = B

                #print(np.interp(theta_o, aux[:, 1], aux[:, 2]))
                #print(aux)

                #break

            
            if output_file.endswith('.npy'):
                # print(final_reflectance.shape)
#                np.save(output_file, final_reflectance)
                n = theta_i.shape[0]
                m = theta_o.shape[0]
                img = final_reflectance[:, 2:].reshape([n, m, 3])
                np.save(output_file, img)
            else: 
                np.savetxt(output_file, final_reflectance, delimiter='\t', fmt='%f')
        else:

            #print(final_reflectance.shape)
            #input()

            if output_file.endswith('.npy'):
                np.save(output_file, reflectance_rgb)
            else: 
                np.savetxt(output_file, reflectance_rgb, delimiter='\t', fmt='%f')

        # Plot reflectance RGB
        n = theta_i.shape[0]
        m = theta_o.shape[0]
        # breakpoint()
        img = final_reflectance[:, 2:].reshape([n, m, 3])

        print("Img shape {0}".format(img.shape))

        ext = output_file[-4: ]

        print("Extension {0}".format(ext))

        #output_file = output_file.replace(ext, '.pdf')
        output_file = output_file.replace(ext, '.png')
        output_file = output_file.replace(ext, '.svg')

        print(output_file)

        #if self.use_back_reflection:
        #    max_theta_o = 120

        self.plotRGBPatch(img, output_file, name, min_theta_o, max_theta_o)

        return reflectance_rgb

    def RGBtoSRGB(self, output_file):
        
        print("testing RGBtosRGB")
        print(self.reflectance_data.shape)
        n, m, c = self.reflectance_data.shape

        #print("Size of the image: {0}x{1}".format(n, m))

        RGB = self.reflectance_data
        #print("RGB {0}".format(RGB.shape))

        XYZ = colour.RGB_to_XYZ(RGB, self.colour_space, None, None)
        sRGB = colour.XYZ_to_sRGB(XYZ)

        #print("sRGB {0}".format(sRGB.shape))

        # Save raw data
        if output_file.endswith('.npy'):
            np.save(output_file, sRGB)
        else: 
            np.savetxt(output_file, sRGB, delimiter='\t', fmt='%f')

        # Plot reflectance sRGB
        """
        # NB: we need this if the format is theta_i, theta_o, R, G, B
        n = theta_i.shape[0]
        m = theta_o.shape[0]

        img = final_reflectance[:, 2:].reshape([n, m, 3])

        """

        ext = output_file[-4: ]

        #print("Extension {0}".format(ext))

        #output_file = output_file.replace(ext, '.pdf')
        #output_file = output_file.replace(ext, '.png')
        output_file = output_file.replace(ext, '.svg')

        #print(output_file)

        self.plotRGBPatch(sRGB, output_file, "Cosmetic Model sRGB")

        return sRGB

    def sRGBtoRGB(self, output_file):
        
        print("testing sRGBtoRGB")

        #print(self.reflectance_data.shape)
        n, m, c = self.reflectance_data.shape

        #print("Size of the image: {0}x{1}".format(n, m))

        sRGB = self.reflectance_data
        #print("RGB {0}".format(RGB.shape))

        XYZ = colour.sRGB_to_XYZ(sRGB)
        RGB = colour.XYZ_to_RGB(XYZ, self.colour_space, None, None)

        #print("sRGB {0}".format(sRGB.shape))

        # Save raw data
        if output_file.endswith('.npy'):
            np.save(output_file, RGB)
        else: 
            np.savetxt(output_file, RGB, delimiter='\t', fmt='%f')

        # Plot reflectance sRGB
        #n = theta_i.shape[0]
        #m = theta_o.shape[0]

        #img = final_reflectance[:, 2:].reshape([n, m, 3])

        ext = output_file[-4: ]

        #print("Extension {0}".format(ext))

        #output_file = output_file.replace(ext, '.pdf')
        output_file = output_file.replace(ext, '.png')

        #print(output_file)

        self.plotRGBPatch(RGB, output_file, "Cosmetic Model sRGB")

        return RGB

parser = argparse.ArgumentParser(description="Fitting model to measurements data")
#parser.add_argument("-transformation", default="RGB", help="Valid options: RGB | sRGB")
parser.add_argument("-input_file", default="./scripts/measurements/Blush_Bulk_coll_norm.txt", help="Input file (Spectrum)")
parser.add_argument("-output_file", default="./scripts/measurements/Blush_Bulk_coll_RGB.txt", help="Output file (RGB transformed)")
parser.add_argument("-illuminant", default="D65", help="Factory illuminant to plot. Options: D65, D50")
parser.add_argument("-colour_space", default="CIE RGB", help="Factory illuminant to plot. Options: CIEG RGB, D50")
parser.add_argument("-transformation", default="SpectrumtosRGB", help="Transformation reflectance type. Options: SpectrumtosRGB, SpectrumtoRGB, RGBtosRGB, sRGBtoRGB")
parser.add_argument("-interpolate", action='store_true')

args = parser.parse_args()

# Colour space availables
# print(sorted(colour.RGB_COLOURSPACES))

input_file = args.input_file
output_file = args.output_file
illuminant = args.illuminant
colour_space = args.colour_space
transformation = args.transformation # Options: RGBtoSRGB
interpolate = args.interpolate
name = ""
#name = "Clarins Foundation 103N Thick Layer Tab"
#name = "Clarins Foundation 108_5W Thick Layer Tab"
#name = "Clarins Foundation 112_5W Thin Layer Tab"


spectrum_helper = SpectrumHelper(transformation = transformation, illuminant = illuminant, colour_space = colour_space)

#file = "./scripts/measurements/Clarins112_5W_thin_layer_on_silicon_norm_tab.txt"
#file = "./scripts/measurements/Highlight_Bulk_coll_norm.txt"
#file = "./scripts/measurements/RimmelBlushThinLayerOnBlack_norm.txt"
#input_file = "./scripts/measurements/Blush_Bulk_coll_norm.txt"
spectrum_helper.readData(input_file)

if transformation == "RGBtosRGB":
    reflectance_rgb = spectrum_helper.RGBtoSRGB(output_file)
else:
    if transformation == "sRGBtoRGB":
        reflectance_rgb = spectrum_helper.sRGBtoRGB(output_file)
    else:
        #output_file = file = "./scripts/measurements/Blush_Bulk_coll_RGB.txt"
        reflectance_rgb = spectrum_helper.toRGB(output_file, interpolate, name)