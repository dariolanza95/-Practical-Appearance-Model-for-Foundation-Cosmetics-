import os
import argparse
def mkdir(folder):
    if not os.path.isdir(folder):
        os.mkdir(folder)

def teaserImgs(spp,outfolder):
    scenes = ["bare","bot","all"]
    mkdir(outfolder)
    for scene in scenes:
        outfile = os.path.join(outfolder,"teaser_{0}_v009.exr".format(scene)) 
        command = "pbrt-v4/build/pbrt  scenes/emily/emily_teaser_v009_{0}.pbrt -spp {1} --outfile {2}".format(scene,spp,outfile) 
        os.system(command=command)
        
def apperanace_matching(spp,outfolder):
    cosmetics = ["105","108"]
    mkdir(outfolder)
    for cosmetic in cosmetics:
        outfile = os.path.join(outfolder,"app_matching_clarins{0}.exr".format(cosmetic)) 
        command = "pbrt-v4/build/pbrt  scenes/emily/Emily_results_scene_forehead_v03_C{0}.pbrt -spp {1} --outfile {2}".format(cosmetic,spp,outfile) 
        os.system(command=command)

parser = argparse.ArgumentParser(description="Given a valid scene, this script will render frames and create a final grid plot where two parameters are variable and the rest are fixed")

parser.add_argument("--task", type=str, default="teaser", help="select the task that the script will perform (teaser | appMatching)")
parser.add_argument("--spp", help="spp used for the task")
parser.add_argument("--out_folder", default=".", help="output folder where to store the results")

args = parser.parse_args()
task = args.task
spp = args.spp
outfolder = args.out_folder
if task == "teaser":     
    if not spp:
        spp = 8192
    teaserImgs(spp,outfolder)
if task == 'app_matching':
    if not spp:
        spp = 1024
    apperanace_matching(spp,outfolder)