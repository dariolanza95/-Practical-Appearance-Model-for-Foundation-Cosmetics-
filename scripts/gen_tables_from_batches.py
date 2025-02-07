import os
import shutil
import glob
import sys
from pathlib import Path
import argparse
import itertools

import json

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
class AppearanceExplorer:
    def __init__(self, verbose, spp, n_threads, nodebind, width, height):
        self.spp = spp
        self.n_threads = n_threads
        self.nodebind = nodebind
        self.width = width
        self.height = height
        self.verbose = verbose
    def createGridPlot(self, input_folder, output_path, format_file, row_name, column_name, row_label, column_label, row_values, column_values):
        # The key idea is plotting the renders on a "table"

        parameters = list(itertools.product(row_values, column_values))

        # Filling the column labels

        cols = []

        for element in column_values:
            cols.append("{0} = {1}".format(column_label, element))
            
        # Filling the row labels
        rows = []

        for element in row_values:
            rows.append("{0} = {1}".format(row_label, element))

        #define grid of plots
        # A key question is how can we pick the fig_size to obtain a reasonable output image

        squeeze = False

        if len(cols) == 1 or len(rows) == 1:
            print("We have only an array, we should squeeze")
            squeeze = True

        fig, axs = plt.subplots(figsize=(8,8), nrows=len(rows), ncols=len(cols), sharex=True, sharey=True, gridspec_kw = {'wspace':0, 'hspace':0}, squeeze = squeeze)

        # Setting the table titles
        for ax, col in zip(axs[0][:], cols):
            ax.set_title(col, size='small')

        for ax, row in zip(axs[:,0], rows):
            ax.set_ylabel(row, rotation=90, size='small')
        
        for i, row_value in enumerate(row_values):

            for j, column_value in enumerate(column_values):

                #axs[i][j].set_axis_off() # turn off the axis
                axs[i][j].set_xticklabels([])
                axs[i][j].set_yticklabels([])
                axs[i][j].tick_params(left=False, bottom=False)

                file = os.path.join(input_folder, "render_{0}_{1}_{2}_{3}.png".format(column_name, column_value, row_name, row_value))
                # file = os.path.join(input_folder, "rendering_{0}_size_{1}_conc_{2}_af_{3}_ad_{4}_g1_{5}_g2_{6}_wg_{7}.png".format(column_name, column_value, row_name, row_value))
                # file = os.path.join("scripts","rendering_batteries", input_folder, "rendering_0.0_size_0.1_conc_{0}_af_0.9_ad_{1}_g1_0.0_g2_0.0_wg_0.5.png".format( row_value,column_value))
                
                # 
                img = mpimg.imread(file)
                axs[i][j].imshow(img)
    
        #fig.savefig(output_path, format="format", bbox_inches="tight", facecolor='white')
        fig.savefig(output_path, format=format_file, bbox_inches="tight")



# Argument parameters of the program (general setup)
def main():
    parser = argparse.ArgumentParser(description="Given a valid scene, this script will render frames and create a final grid plot where two parameters are variable and the rest are fixed")
    parser.add_argument('--scene', type=str, default='./scenes/sphere_test/sphere_elliptical_video.xml', help='Mitsuba scene to be rendered')
    parser.add_argument("-spp", "--numSamples", type=int,default = 64, help="set the number of samples per pixel")
    parser.add_argument("-t", "--threads", type=int,default = 8, help="set the number of threads to be used")
    parser.add_argument("-width", "--width", type=int, default=512, help="width of the image")
    parser.add_argument("-height", "--height", type=int, default=512, help="height of the image")
    parser.add_argument("-out", "--output_folder", default="./results", help="output folder with the final renders")   
    parser.add_argument("--nodebind", "-node", action="store_true", help="use NUMACTL option -cpunodebind=0 ")

    # Argument parameters of the program (specific experiments) 
    parser.add_argument('--conf_file', type=str, default='./scripts/rendering_batteries/conf_file.json', help='Specific configuration file with all the information to create the final plots')

    args = parser.parse_args()

    # output_folder = args.output_folder
    output_folder = "./results"
    output_path = Path(output_folder)

    scene_name = args.scene
    #experiment_name = args.experiment
    verbose = True
    render = True
    subplots = True
    grid_plot = False
    plot_labels = True

    # Performing a particular experiment given the json configuration file
    configuration_folder = ""

    explorer = AppearanceExplorer(verbose, args.numSamples, args.threads, args.nodebind, args.width, args.height)

    #input_conf_path = "./scripts/appearance_exploration/thin_film.json"
    input_conf_path = args.conf_file
    # explorer.performExperiment(input_conf_path, output_folder)

    file = open (input_conf_path, "r")  
    conf = json.loads(file.read())
    file.close()

            # Extract the specific parameters
    row_values = conf["row_values"]
    column_values = conf["column_values"]

    print(row_values)
    print(column_values)

    row_name = conf["row_name"]
    column_name = conf["column_name"]

    row_label = conf["row_label"]
    column_label = conf["column_label"]

    experiment_name = conf["experiment_name"]
    plot_type = conf["plot_type"]
    output_plot_path = conf["output_plot_path"]
    #make_renders = conf["make_renders"]
    format_file = conf["format_file"]
    render_folder = conf["render_folder"]
    temp_pa = os.path.join("scripts","rendering_batteries", output_plot_path)
    if len(row_values) > 1 and len(column_values) > 1:
        explorer.createGridPlot(render_folder, temp_pa, format_file, row_name, column_name, row_label, column_label, row_values, column_values)

main()