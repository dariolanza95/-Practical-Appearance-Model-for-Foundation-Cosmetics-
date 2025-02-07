import os
def replaceContent(line_numbers_to_change,lines,new_contents,out_file_path):
    for i, line_number in enumerate(line_numbers_to_change):
        if line_number < len(lines):
            print(new_contents[i])
            lines[line_number] = new_contents[i] + '\n'
    #  = "../../pbrt-v4/build/pbrt simple_sphere_edited_test.pbrt"
    
    with open(out_file_path, 'w') as file:
        file.writelines(lines)
    return out_file_path

#this script is very dumb, but was the fastest thing that came to my mind
def editSceneTeaser(scene_dict,lines):
    line_numbers_to_change = [425-1,426-1,]
    new_content_for_lines = []
    #flakes for top layer definition at line 425
    new_content_for_lines.append( "\"rgb albedo_flakes_top \" [{0} {1} {2}]".format(scene_dict["albedo_flake_top"][0],scene_dict["albedo_flake_top"][1],scene_dict["albedo_flake_top"][2]))
    #spherical diffusers for top layer  definition at line 426
    new_content_for_lines.append( "\"rgb albedo_diff_top \" [{0} {1} {2}]".format(scene_dict["albedo_diff_top"][0],scene_dict["albedo_diff_top"][1],scene_dict["albedo_diff_top"][2]))
    
    return replaceContent(line_numbers_to_change=line_numbers_to_change,lines=lines,new_contents=new_content_for_lines, out_file_path=os.path.join("scenes","emily","emily_teaser_v009_all_expl.pbrt"))

def editSceneArm(scene_dict,lines):
    line_numbers_to_change = [24-1, 72-1, 173-1, 216-1,226-1,228-1,229-1,231-1,232-1,233-1,234-1,236-1,237-1,238-1,239-1,240-1]
    new_content_for_lines = []
    
    #change the cropwindow 24
    new_content_for_lines.append("\"float cropwindow\" [{0} {1} {2} {3}]".format( scene_dict["cropwindow"][0],scene_dict["cropwindow"][1],scene_dict["cropwindow"][2],scene_dict["cropwindow"][3]))
    #change the skin map 72
    new_content_for_lines.append("\"string filename\" [\"{0}\"]".format(scene_dict["albedo_map"]))
    #cosmetic_map; line 173
    new_content_for_lines.append("\"string filename\" [\"{0}\"]".format(scene_dict["cosmetic_map"]))

    #section for changing the material
    #etaSkin; line 216
    new_content_for_lines.append("\"float etaSkin\" [{0}]".format(scene_dict["etaSkin"]))
    #thickness; line 226
    new_content_for_lines.append("\"float thick_scale\" [{0}]".format(scene_dict["thickness"]))
    #albedo diffuse line 228
    new_content_for_lines.append( "\"rgb albedo_diff \" [{0} {1} {2}]".format(scene_dict["albedo_diff"][0],scene_dict["albedo_diff"][1],scene_dict["albedo_diff"][2]))
    #albedo flake line 229
    new_content_for_lines.append( "\"rgb albedo_flakes \" [{0} {1} {2}]".format(scene_dict["albedo_flake"][0],scene_dict["albedo_flake"][1],scene_dict["albedo_flake"][2]))
    #particle size; line 231
    new_content_for_lines.append("\"float sigma\" [{0}]".format(scene_dict["sigma"]))
    #particles orientations; line 232
    new_content_for_lines.append("\"float mean_p_orientation\" [{0}]".format(scene_dict["mean"]))
    #particles std size; line 233
    new_content_for_lines.append("\"float stdSize\" [{0}]".format(scene_dict["stdSize"]))
    #particles std_or orientations; line 234
    new_content_for_lines.append("\"float std_or\" [{0}]".format(scene_dict["std_or"]))
    #particles concentrations; line 236
    new_content_for_lines.append("\"float diff_particles_perc\" [{0}]".format(scene_dict["concentration"]))
    
    #diffuser forward; line 237
    #new_content_for_lines.append("\"float eta\" [{0}]".format(scene_dict["ior"]))
    #diffuser forward; line 238
    #new_content_for_lines.append("\"float roughness\" [{0}]".format(scene_dict["surf_rough"]))

    
    #diffuser forward; line 237
    new_content_for_lines.append("\"float g1\" [{0}]".format(scene_dict["g1"]))
    #diffuser backward; line 238
    new_content_for_lines.append("\"float g2\" [{0}]".format(scene_dict["g2"]))
    #diffuser forward/backward weight; line 239
    new_content_for_lines.append("\"float w_g\" [{0}]".format(scene_dict["w_g"]))
    #maxDepth; line 240
    new_content_for_lines.append("\"integer maxdepth\" [{0}]".format(int(scene_dict["maxdepth"])))
    return replaceContent(line_numbers_to_change=line_numbers_to_change,lines=lines,new_contents=new_content_for_lines, out_file_path=os.path.join("scenes","arm","arm_scene_edited.pbrt"))

def editSceneEmily(scene_dict,lines):
    line_numbers_to_change = [216-1,226-1,228-1,229-1,231-1,232-1,233-1,234-1,236-1,237-1,238-1,239-1,240-1]
    # line_numbers_to_change = [216,226,228,229,231,232,233,234,236,237,238,239,240]
    new_content_for_lines = []
    #etaSkin; line 216
    new_content_for_lines.append("\"float etaSkin\" [{0}]".format(scene_dict["etaSkin"]))
    #thickness; line 226
    new_content_for_lines.append("\"float thick_scale\" [{0}]".format(scene_dict["thick_scale"]))
    #albedo diffuse line 228
    new_content_for_lines.append( "\"rgb albedo_diff \" [{0} {1} {2}]".format(scene_dict["albedo_diff"][0],scene_dict["albedo_diff"][1],scene_dict["albedo_diff"][2]))
    #albedo flake line 229
    new_content_for_lines.append( "\"rgb albedo_flakes \" [{0} {1} {2}]".format(scene_dict["albedo_flake"][0],scene_dict["albedo_flake"][1],scene_dict["albedo_flake"][2]))
    #particle size; line 231
    new_content_for_lines.append("\"float sigma\" [{0}]".format(scene_dict["sigma"]))
    #particles orientations; line 232
    new_content_for_lines.append("\"float mean_p_orientation\" [{0}]".format(scene_dict["mean"]))
    #particles std size; line 233
    new_content_for_lines.append("\"float stdSize\" [{0}]".format(scene_dict["stdSize"]))
    #particles std_or orientations; line 234
    new_content_for_lines.append("\"float std_or\" [{0}]".format(scene_dict["std_or"]))
    #particles concentrations; line 236
    new_content_for_lines.append("\"float diff_particles_perc\" [{0}]".format(scene_dict["diff_particles_perc"]))
    
    #diffuser forward; line 237
    #new_content_for_lines.append("\"float eta\" [{0}]".format(scene_dict["ior"]))
    #diffuser forward; line 238
    #new_content_for_lines.append("\"float roughness\" [{0}]".format(scene_dict["surf_rough"]))

    
    #diffuser forward; line 237
    new_content_for_lines.append("\"float g1\" [{0}]".format(scene_dict["g1"]))
    #diffuser backward; line 238
    new_content_for_lines.append("\"float g2\" [{0}]".format(scene_dict["g2"]))
    #diffuser forward/backward weight; line 239
    new_content_for_lines.append("\"float w_g\" [{0}]".format(scene_dict["w_g"]))
    #maxDepth; line 240
    new_content_for_lines.append("\"integer maxdepth\" [{0}]".format(int(scene_dict["maxdepth"])))
    return replaceContent(line_numbers_to_change=line_numbers_to_change,lines=lines,new_contents=new_content_for_lines, out_file_path=os.path.join("scenes","emily","Emily_results_scene_edited.pbrt"))
# def editFile(mean,incident,s,conc,lines):
def editFile(scene_dict,lines):

# line_numbers_to_change = [8-1,34-1,50-1,51-1,55-1]
    # line_numbers_to_change = [10-1,46-1,47-1,51-1]
    line_numbers_to_change = [42-1,44-1,45-1,49-1,50-1,54-1,55-1,56-1,57-1,58-1,59-1,63-1]
    new_content_for_lines = []
    #it's important the order otherwise things wont align
    #name
    # new_content_for_lines.append("\"string filename\" \"simple_sphere_mean_{0}_size_{1}_conc_{2}_af_{3}_ad_{4}_g1_{5}_g2_{6}_wg_{7}.png\"".format(mean,s,conc,albedo_flake,albedo_diff,g1,g2,w_g))
    #Incident light Rotation
    # new_content_for_lines.append("Rotate {0} 1 0 0".format(incident))
    #thickness; line 42
    new_content_for_lines.append("\"float thickness\" [{0}]".format(scene_dict["thickness"]))
    #albedo diffuse line 44
    new_content_for_lines.append( "\"rgb albedo_diff \" [{0} {1} {2}]".format(scene_dict["albedo_diff"][0],scene_dict["albedo_diff"][1],scene_dict["albedo_diff"][2]))
    #albedo flake line 45
    new_content_for_lines.append( "\"rgb albedo_flakes \" [{0} {1} {2}]".format(scene_dict["albedo_flake"][0],scene_dict["albedo_flake"][1],scene_dict["albedo_flake"][2]))
    #particle size; line 49
    new_content_for_lines.append("\"float sigma\" [{0}]".format(scene_dict["sigma"]))
    #particles orientations; line 50
    new_content_for_lines.append("\"float mean\" [{0}]".format(scene_dict["mean"]))
    #particles concentrations; line 54
    new_content_for_lines.append("\"float diff_particles_perc\" [{0}]".format(scene_dict["concentration"]))
    #diffuser forward; line 55
    new_content_for_lines.append("\"float eta\" [{0}]".format(scene_dict["ior"]))
    #diffuser forward; line 56
    new_content_for_lines.append("\"float roughness\" [{0}]".format(scene_dict["surf_rough"]))

    
    #diffuser forward; line 57
    new_content_for_lines.append("\"float g1\" [{0}]".format(scene_dict["g1"]))
    #diffuser backward; line 58
    new_content_for_lines.append("\"float g2\" [{0}]".format(scene_dict["g2"]))
    #diffuser forward/backward weight; line 59
    new_content_for_lines.append("\"float w_g\" [{0}]".format(scene_dict["w_g"]))
    #beneath layer color forward/backward weight; line 63

    new_content_for_lines.append("\"rgb reflectance\" [{0} {1} {2}]".format(scene_dict["botLayer"][0],scene_dict["botLayer"][1],scene_dict["botLayer"][2]))



    out_file_path = os.path.join("scenes","pbrt","test_ball_edited.pbrt")

    return replaceContent(line_numbers_to_change=line_numbers_to_change,lines=lines,new_contents=new_content_for_lines,out_file_path=out_file_path)

def get_scene_dictionary(params_range):
        final_dict = {}

        final_dict["mean"] = params_range["mean"][0]

        if "thickness" in params_range.keys():
            final_dict["thickness"] = params_range["thickness"][0]
        else:
            final_dict["thickness"] = 8.0
        
        if "surf_rough" in params_range.keys():
            final_dict["surf_rough"] = params_range["surf_rough"][0]
        else:
            final_dict["surf_rough"] = 0.0
        
        if "ior" in params_range.keys():
            final_dict["ior"] = params_range["ior"][0]
        else:
            final_dict["ior"] = 0.0

        if "botLayer" in params_range.keys():
            final_dict["botLayer"] = params_range["botLayer"][0]
        else:
            final_dict["botLayer"] = [ 0,0,0]

        if "etaSkin" in params_range.keys():
            final_dict["etaSkin"] = params_range["etaSkin"][0]
        else:
            final_dict["etaSkin"] = 1.3
        if "stdSize" in params_range.keys():
            final_dict["stdSize"] = params_range["stdSize"][0]
        else:
            final_dict["stdSize"] = 0.0
        
        if "std_or" in params_range.keys():
            final_dict["std_or"] = params_range["std_or"][0]
        else:
            final_dict["std_or"] = 0.0
        if "maxdepth" in params_range.keys():
            final_dict["maxdepth"] = params_range["maxdepth"][0]
        else:
            final_dict["maxdepth"] = 512
        
        if "albedo_map" in params_range.keys():
            final_dict["albedo_map"] = params_range["albedo_map"][0]
        else:
            final_dict["albedo_map"] = "./hand_textures/albedo.tga"
        

        final_dict["incident_or"] = params_range["incident_or"][0]
        final_dict["concentration"] = params_range["concentration"][0]
        final_dict["sigma"] = params_range["sigma"][0]

        final_dict["g1"] = params_range["g1"][0]
        final_dict["g2"] = params_range["g2"][0]
        final_dict["w_g"] = params_range["w_g"][0]
        final_dict["albedo_flake"] = params_range["albedo_flake"][0]
        final_dict["albedo_diff"] = params_range["albedo_diff"][0]
        final_dict["experiment_name"] = params_range["experiment_name"][0]

    
        # Modify the specific lines
        row_name = params_range["row_name"]
        column_name = params_range["column_name"]
        if column_name == "color":
            print("skipping")
            # final_dict[column_name] = params_range[column_name]
        else:
            final_dict[column_name] = params_range[column_name]

        final_dict[row_name] = params_range[row_name]
        return final_dict

import os,re
def renderingRoutine(column_name,column_value,row_name,row_value,experiment_name,out_folder = "",scene_path="scenes/pbrt/test_ball_edited.pbrt",nodebind=-1,arm_idx=-1):
    import subprocess
    # Define the shell command you want to run
    if column_name == "color" and arm_idx != -1:
        pattern = r"skin_type_[IVX]+"
        match = re.search(pattern, row_value)
        print(row_value)
        print(match)
        skin_type = match[0]
        out_name = "render_{0}_{1}_{2}_{3}.exr".format(column_name, arm_idx, row_name, skin_type)
    elif column_name == "color" and arm_idx == -1:
        "Warning!! You haven't passed the arm_idx, this means that all the rendering for the same column value will have the same name."
    else:
        out_name = "render_{0}_{1}_{2}_{3}.exr".format(column_name, column_value, row_name, row_value)

    # out_name = "render_{0}_{1}_{2}_{3}.png".format(column_name, column_value, row_name, row_value)
                

    # out_name = "rendering_{0}_size_{1}_conc_{2}_af_{3}_ad_{4}_g1_{5}_g2_{6}_wg_{7}.png".format(mean,s,conc,albedo_flake[0][0],albedo_diff[0][0],g1,g2,w_g)

    if out_folder != "":
        if not os.path.exists(os.path.join("scripts","rendering_batteries",out_folder,experiment_name)):
            os.mkdir(os.path.join("scripts","rendering_batteries",out_folder,experiment_name))
        out_path = os.path.join("scripts","rendering_batteries",out_folder,experiment_name,out_name)
    else:
        if not os.path.exists(os.path.join("scripts","rendering_batteries",experiment_name)):
            os.mkdir(os.path.join("scripts","rendering_batteries",experiment_name))
        out_path = os.path.join("scripts","rendering_batteries",experiment_name,out_name)
    if nodebind == -1:
        command = "./pbrt-v4/build/pbrt {1} --outfile {0}".format(out_path,scene_path)  # Replace with the command you want to execute
    else:
        command = "numactl --cpunodebind {2} ./pbrt-v4/build/pbrt {1} --outfile {0}".format(out_path,scene_path,nodebind)  # Replace with the command you want to execute
    print("**Shell command**")
    print(command)
    # Run the shell command
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return out_path
def gen_cosmetic_pos_dict():
    #these are hardcoded value found by trial and error
    cosmetic_pos_dict = {}
    
    #cosmetic_pos_dict[0] = {"cropwindow": [0.55, 0.80, 0.42, 0.63] , "cosmetic_map": "./hand_textures/cosmetic_map_v04_pos_01.png"}
    #cosmetic_pos_dict[1] = {"cropwindow": [0.80, 0.85, 0.42, 0.63] , "cosmetic_map": "./hand_textures/cosmetic_map_v04_pos_02.png"}
    #cosmetic_pos_dict[2] = {"cropwindow": [0.85, 0.98, 0.42, 0.63] , "cosmetic_map": "./hand_textures/cosmetic_map_v04_pos_03.png"}
    cosmetic_pos_dict[0] = {"cropwindow": [0.57, 0.70, 0.45, 0.58] , "cosmetic_map": "./hand_textures/cosmetic_map_v06_pos_01.png"}
    cosmetic_pos_dict[1] = {"cropwindow": [0.70, 0.76, 0.45, 0.58] , "cosmetic_map": "./hand_textures/cosmetic_map_v06_pos_01.png"}
    cosmetic_pos_dict[2] = {"cropwindow": [0.76, 0.83, 0.45, 0.58] , "cosmetic_map": "./hand_textures/cosmetic_map_v06_pos_01.png"}
    cosmetic_pos_dict[3] = {"cropwindow": [0.83, 0.95, 0.45, 0.58] , "cosmetic_map": "./hand_textures/cosmetic_map_v06_pos_01.png"}
    
    
    return cosmetic_pos_dict
#default
#"float cropwindow" [0.55 0.98 0.42 0.63]

# valid for position 01
#"float cropwindow" [0.55 0.80 0.42 0.63]

# valid for position 02
#"float cropwindow" [0.80 0.85 0.42 0.63]


# valid for position 03
#"float cropwindow" [0.85 0.93 0.42 0.63]


from PIL import Image
import cv2
import OpenEXR,numpy as np
import Imath
def readExr(path):
    exr_file = OpenEXR.InputFile(path)
    header = exr_file.header()
    dw = header['dataWindow']
    size = (dw.max.x - dw.min.x + 1, dw.max.y - dw.min.y + 1)
    print(path)
    # Read the RGB channels
    # breakpoint()

    R = np.frombuffer(exr_file.channel('R',Imath.PixelType(Imath.PixelType.FLOAT)), dtype=np.float32).reshape(size[1], size[0])
    G = np.frombuffer(exr_file.channel('G',Imath.PixelType(Imath.PixelType.FLOAT)), dtype=np.float32).reshape(size[1], size[0])
    B = np.frombuffer(exr_file.channel('B',Imath.PixelType(Imath.PixelType.FLOAT)), dtype=np.float32).reshape(size[1], size[0])
    # Combine channels into a single image
    image = cv2.merge([B, G, R])
    return image
def writeExr(output_path,stitched_image):
    stitched_image = stitched_image.astype(np.float32)
    exr = OpenEXR.OutputFile(output_path, OpenEXR.Header(stitched_image.shape[1], stitched_image.shape[0]))
    exr.writePixels({'R': stitched_image[:, :, 2].tobytes(), 'G': stitched_image[:, :, 1].tobytes(), 'B': stitched_image[:, :, 0].tobytes()})
    exr.close()

import natsort
def stitch_v_imgs(out_folder,ext='exr'):
    images = []
    
    row_paths = glob.glob(os.path.join(out_folder)  +  "/row_*.{0}".format(ext))
    row_paths = natsort.natsorted(row_paths)
    for path in row_paths:
        if os.path.exists(path):
            ext = path.split('.')[-1].lower()
            if ext == '.png':
                img = cv2.imread(path)
            elif ext == 'exr':
                print("readingExr")
                img = readExr(path)
            images.append(img) 
    stitched_image = cv2.vconcat(images)
    if ext == '.png':
        cv2.imwrite(os.path.join(os.path.dirname(out_folder),"final_img.png"), stitched_image)
    elif ext == 'exr':
        print("readingExr")
        img = writeExr(os.path.join(os.path.dirname(out_folder),"final_img.exr"),stitched_image)

def stitch_imgs(out_paths,output_path):
    #take the path images and stich them together to create a single image
    #images = [Image.open(path) for path in out_paths if os.path.exists(path)]
    #images = [cv2.imread(path) for path in out_paths if os.path.exists(path)]
    images = []
    for path in out_paths:
        if os.path.exists(path):
            ext = path.split('.')[-1].lower()
            if ext == '.png':
                img = cv2.imread(path)
            elif ext == 'exr':
                print("readingExr")
                img = readExr(path)
            images.append(img) 
    # Stitch images horizontally
    stitched_image = cv2.hconcat(images)
    if ext == '.png':
        cv2.imwrite(output_path, stitched_image)
    elif ext == 'exr':
        print("readingExr")
        img = writeExr(output_path,stitched_image)

    return stitched_image

def params_exploration_arm(params_range_path = "scripts/test_params_range.json",out_folder = "",nodebind=-1):
    import os
    scene_path = os.path.join("scenes","arm","arm_scene_default.pbrt")
    with open(scene_path, 'r') as file:
        lines = file.readlines()

    import json
    out_paths = []
    print(params_range_path)
    with open(params_range_path, 'r') as file:
            params_range = json.load(file)
    # print(params_range)
    scene_dict = get_scene_dictionary(params_range)
    file_path = "test_param_ranges.json"
    experiment_name = params_range["experiment_name"]
    row_name = params_range["row_name"]
    column_name = params_range["column_name"]
    row_values = params_range["row_values"]
    column_values = params_range["column_values"]
    cosmetic_pos_dict = gen_cosmetic_pos_dict()
    for row_idx,row_value in enumerate(row_values):
        for col_idx,column_value in enumerate(column_values):
            # params_range[column_name]
            row_name = params_range["row_name"]
            column_name = params_range["column_name"]
            if column_name == "color":
                scene_dict["albedo_diff"] = column_value[0]
                scene_dict["albedo_flake"] = column_value[1]
            else:    
                scene_dict[column_name] = column_value
            scene_dict[row_name] = row_value
            scene_dict["cropwindow"] = cosmetic_pos_dict[col_idx]["cropwindow"]
            scene_dict["cosmetic_map"] =  cosmetic_pos_dict[col_idx]["cosmetic_map"]
            print("******* NEW ITERATION *******")
            print(" col_name {0} col val {1}".format(column_name,column_value) )
            print(" row_name {0} row val {1}".format(row_name,row_value) )
            print(scene_dict["cosmetic_map"])
            
            # edited_scene_path = editSceneArm(scene_dict,lines=lines)            
            # out_path = renderingRoutine(column_name,column_value,row_name,row_value,experiment_name,out_folder,scene_path= "scenes/arm/arm_scene_edited.pbrt", nodebind=nodebind,arm_idx = col_idx)
            # out_paths.append(out_path)
        ext = 'exr'
        # row_output_path = os.path.join("scripts","rendering_batteries",out_folder,experiment_name,"row_{0}.{1}".format(row_idx,ext))
        # stitch_imgs(out_paths,row_output_path)
        out_paths.clear()
    img_out_path = os.path.join("scripts","rendering_batteries",out_folder,experiment_name)    
    stitch_v_imgs(img_out_path)
def params_exploration(params_range_path = "scripts/test_params_range.json",out_folder = "",nodebind=-1):
    import os
    scene_path = os.path.join("scenes","pbrt","test_ball.pbrt")
    with open(scene_path, 'r') as file:
        lines = file.readlines()

    import json
    print(params_range_path)
    with open(params_range_path, 'r') as file:
            params_range = json.load(file)
    # print(params_range)
    scene_dict = get_scene_dictionary(params_range)
    file_path = "test_param_ranges.json"
    experiment_name = params_range["experiment_name"]
    row_name = params_range["row_name"]
    column_name = params_range["column_name"]
    row_values = params_range["row_values"]
    column_values = params_range["column_values"]
    for column_value in column_values:
        for row_value in row_values:
            params_range[column_name]
            
            row_name = params_range["row_name"]
            column_name = params_range["column_name"]

            scene_dict[column_name] = column_value
            scene_dict[row_name] = row_value
            print("******* NEW ITERATION *******")
            print(" row_name {0} row val {1}".format(row_name,row_value) )
            print(" col_name {0} col val {1}".format(column_name,column_value) )
            
            
            edited_scene_path = editFile(scene_dict,lines=lines)
            
            renderingRoutine(column_name,column_value,row_name,row_value,experiment_name,out_folder,nodebind=nodebind)
    
    # for mean in means:
    #         for conc in concentrations:
    #             for s in sigmas:
    #                 for g1 in g1s:
    #                     for g2 in g2s:
    #                         for w_g in w_gs:
    #                             for af in albedo_flake:
    #                                 for ad in albedo_diff:
    #                                     edited_scene_path = editFile(mean=mean,g1=g1,g2=g2,w_g=w_g,albedo_flake=af,albedo_diff=ad,conc=conc,s=s,lines=lines)
    #                                     # rendering_routine(mean,s,conc,albedo_flake,albedo_diff,g1,g2,w_g,experiment_name)

                # render(edited_scene_path)
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
    def createGridPlot(self, input_folder, output_path, format_file, row_name, column_name, row_label, column_label, row_values, column_values,m_title,batch_folder):
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
        
        plt.suptitle(m_title)
        
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
                
                file = os.path.join(batch_folder,input_folder, "render_{0}_{1}_{2}_{3}.png".format(column_name, column_value, row_name, row_value))
                print(file)
                img = mpimg.imread(os.path.join("scripts","rendering_batteries",file))
                axs[i][j].imshow(img)
    
        #fig.savefig(output_path, format="format", bbox_inches="tight", facecolor='white')
        # fig.savefig(output_path, format=format_file, bbox_inches="tight")
        print(output_path)        
        fig.savefig(output_path, format=format_file, bbox_inches="tight")
        


# Argument parameters of the program (general setup)
import json
import itertools
def genTable(config_path,batch_folder):
    explorer = AppearanceExplorer(False, 128, 8, True, 256, 256)
    # explorer = AppearanceExplorer(verbose, args.numSamples, args.threads, args.nodebind, args.width, args.height)

    #input_conf_path = "./scripts/appearance_exploration/thin_film.json"
    input_conf_path = config_path
    # explorer.performExperiment(input_conf_path, output_folder)
    print(input_conf_path)
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
    # output_plot_path = conf["output_plot_path"]
    output_plot_path = os.path.join("./scripts/rendering_batteries/results/",experiment_name)
    print(output_plot_path)
    
#    make_renders = conf["make_renders"]
    format_file = conf["format_file"]
    format_file = "png"
    render_folder = experiment_name
    
    if "means_or" in conf.keys():
        means_or = conf["means_or"]
    else:
        means_or = 0.0
    albedo_diff =conf["albedo_diff"][0]
    albedo_flake = conf["albedo_flake"][0]
    print(albedo_flake)
    title    = "Spectral BRDF Plot: $\\sigma$ = {0:0.3f}, $D \\%$ = {1:0.2f}, $t$ = {12:0.2f} \
$g_1$ = {2:0.2f}, $g_2$ = {3:0.2f}, $w_g$ = {4:0.2f}, \n $\\mu_\\theta$ = {5:0.2f}, $\\alpha_d$ = {6: 0.2f} {7: 0.2f} {8: 0.2f}, $\\alpha_p$ = {9: 0.2f} {10: 0.2f} {11: 0.2f}".format(conf["sigma"][0], conf["concentration"][0], conf["g1"][0], conf["g2"][0], conf["w_g"][0], means_or,albedo_diff[0],albedo_diff[1],albedo_diff[2] ,albedo_flake[0],albedo_flake[1],albedo_flake[2],conf["thickness"][0] )
    
    if len(row_values) > 1 and len(column_values) > 1:
        explorer.createGridPlot(render_folder, output_plot_path, format_file, row_name, column_name, row_label, column_label, row_values, column_values,title,batch_folder)


def convertOptParamsToMaterialParams(param_dict):
    """Function to transform from the optical parameters saved as json to the values used as input in the pbrt file scene"""
    out_dict = {}
    out_dict["etaSkin"] = 1.3
    out_dict["albedo_diff"] = [param_dict["albedo_diff_r"],param_dict["albedo_diff_g"],param_dict["albedo_diff_b"]]
    out_dict["albedo_flake"] = [param_dict["albedo_flake_r"],param_dict["albedo_flake_g"],param_dict["albedo_flake_b"]]
    
    out_dict["thick_scale"] = param_dict["thickness"]
    
    #platelets
    out_dict["mean"] = param_dict["means_or"]
    out_dict["std_or"] = param_dict["stds_or"]
    
    out_dict["sigma"] = param_dict["mean_size"]
    out_dict["stdSize"] = param_dict["stds_size"]

    out_dict["diff_particles_perc"] = param_dict["diff_particles_percs"]
    
    #diffusers
    out_dict["g1"] = param_dict["g1"]
    out_dict["g2"] = param_dict["g2"]
    out_dict["w_g"] = param_dict["w_g"]
    
    out_dict["maxdepth"] = param_dict["n_samples"]
    return out_dict

from itertools import product
def generate_cartesian_product(lists):
    """
    Generate the Cartesian product of multiple lists.
    
    Parameters:
    - lists: A list of input lists.
    
    Returns:
    A list of tuples representing the Cartesian product.
    """
    return list(product(*lists))

def paramsExplorationTeaser(out_folder,scene_path):
    scene_path = os.path.join("scenes","emily","emily_teaser_v009_all.pbrt")
    with open(scene_path, 'r') as file:
        lines = file.readlines()
    import numpy as np
    listR = np.linspace(0.8,0.99,4)
    listG = np.linspace(0.2,0.5,5)
    listB = np.linspace(0.2,0.5,5)

    cartesian_product_list = generate_cartesian_product([listR, listG, listB])
    cartesian_product_list = [list(t) for t in cartesian_product_list]
    print(cartesian_product_list)

    for _list in cartesian_product_list:
        dict_samples_opt_prop = {}
        dict_samples_opt_prop["albedo_flake_top"] = _list
        dict_samples_opt_prop["albedo_diff_top"] = _list
        edited_scene_path = editSceneTeaser(dict_samples_opt_prop,lines=lines)
        column_name = "flakes"
        column_value = "R_{0:.2f}_G_{1:0.2f}_B_{1:0.2f}".format( dict_samples_opt_prop["albedo_flake_top"][0],dict_samples_opt_prop["albedo_flake_top"][1],dict_samples_opt_prop["albedo_flake_top"][2])
        row_name = "diff"
        row_value = "R_{0}_G_{1}_B_{2}".format(dict_samples_opt_prop["albedo_diff_top"][0],dict_samples_opt_prop["albedo_diff_top"][1],dict_samples_opt_prop["albedo_diff_top"][2])
        experiment_name = "teaser_v009_all_explor_v02"
        renderingRoutine(column_name,column_value,row_name,row_value,experiment_name,out_folder,edited_scene_path)



def paramsExplorationSkin(out_folder,scene_path):
    """function used to edit and render the various scenes with the input parameters."""
    
    # scene_path = os.path.join("scenes","emily","emily_sss_setup.pbrt")
    # scene_path = os.path.join("scenes","emily","emily_ablation_study_skin_v02.pbrt")
    with open(scene_path, 'r') as file:
        lines = file.readlines()

    dict_samples_opt_prop = {}
    import json
    foundations = ["Clarins103","Clarins105","Clarins108","Clarins112"]
    for foundation in foundations:
        path = os.path.join("optimized_foundations",foundation,"best_result.json")
        # path = os.path.join("optimized_foundations",foundation,"best_result.json")
        with open(path, 'r') as file:
            opt_params = json.load(file)
        sample_opt_properties_estimated = convertOptParamsToMaterialParams(opt_params)
        sample_opt_properties_estimated['maxdepth'] = 1024
        sample_opt_properties_estimated['nsamples'] = 16
        print(sample_opt_properties_estimated)
        dict_samples_opt_prop[foundation] = sample_opt_properties_estimated
    

    # thickness_levels = [16]
    thickness_levels = [0.0,0.25,0.5,0.75,1.0]
    for foundation in foundations:
        for thickness in thickness_levels:
            dict_samples_opt_prop[foundation]['thick_scale'] = thickness
            edited_scene_path = editSceneEmily(dict_samples_opt_prop[foundation],lines=lines)
            column_name = "Clarins"
            column_value = foundation[-3:]
            row_name = "thick"
            row_value = thickness
            experiment_name = "Fig07"
            renderingRoutine(column_name,column_value,row_name,row_value,experiment_name,out_folder,edited_scene_path)
    exit()
    params_range = {}
    # print(params_range)
    scene_dict = get_scene_dictionary(params_range)
    file_path = "test_param_ranges.json"
    experiment_name = params_range["experiment_name"]
    row_name = params_range["row_name"]
    column_name = params_range["column_name"]
    row_values = params_range["row_values"]
    column_values = params_range["column_values"]
    for column_value in column_values:
        for row_value in row_values:
            params_range[column_name]
            
            row_name = params_range["row_name"]
            column_name = params_range["column_name"]

            scene_dict[column_name] = column_value
            scene_dict[row_name] = row_value
            print("******* NEW ITERATION *******")
            print(" row_name {0} row val {1}".format(row_name,row_value) )
            print(" col_name {0} col val {1}".format(column_name,column_value) )
            
            
            edited_scene_path = editFile(scene_dict,lines=lines)
            
            renderingRoutine(column_name,column_value,row_name,row_value,experiment_name,out_folder)

import os,glob,argparse
def main():
    parser = argparse.ArgumentParser(description="Fitting model to measurements data")
#parser.add_argument("-transformation", default="RGB", help="Valid options: RGB | sRGB")
    parser.add_argument("-scene", default="emily", help="Type of scene that should be rendered. Legal values are (emily|knob)")
    parser.add_argument("-nodebind", type= int , default="-1", help="Nodebind that should be used (i.e. 0,1). If -1 is used then all possible nodes will be used. Default -1")
    parser.add_argument("-batch",  default="batch_13", help="If knob is selected, then select from which folder the configuration files should be taken.")
    
    args = parser.parse_args()

# Colour space availables
# print(sorted(colour.RGB_COLOURSPACES))
    nodebind = args.nodebind
    scene = args.scene
    if scene == "emily_frontal":

        # scene_path = os.path.join("scenes","emily","emily_ablation_study_skin_v02.pbrt ") 
        scene_path = os.path.join("scenes","emily","emily_ablation_study_skin.pbrt") 
        out_folder = ""
        paramsExplorationSkin(out_folder=out_folder,scene_path=scene_path) 
    elif scene == "emily_side":
        scene_path = os.path.join("scenes","emily","emily_sss_setup.pbrt") 
        out_folder = ""
        paramsExplorationSkin(out_folder=out_folder,scene_path=scene_path) 

    elif scene == "teaser":

        scene_path = os.path.join("scenes","emily","emily_teaser_v009_all.pbrt") 
        out_folder = "teaser_exploration" 
        
        paramsExplorationTeaser(out_folder=out_folder,scene_path=scene_path)
    elif scene == 'knob':    
        batch_folder = args.batch
        config_paths = glob.glob(os.path.join("scripts","rendering_batteries", batch_folder)  +  "/*.json")

        for config_path in config_paths:
            bn = os.path.basename(config_path) 
            if bn == "surfRoug_vs_ior_con_0.8_g1_0_pSize_0.1" + ".json":
                print("skipping")
                continue
            
            params_exploration(config_path,batch_folder,nodebind=nodebind)
            genTable(config_path=config_path,batch_folder= batch_folder)
    elif scene == "arm":
        batch_folder = args.batch
        config_paths = glob.glob(os.path.join("scripts","rendering_batteries", batch_folder)  +  "/*.json")
        print(config_paths)
        for config_path in config_paths:
            bn = os.path.basename(config_path) 
            if bn == "surfRoug_vs_ior_con_0.8_g1_0_pSize_0.1" + ".json":
                print("skipping")
                continue
            
            params_exploration_arm(config_path,batch_folder,nodebind=nodebind)
            
        
        # for config_name in config_names:
            # config_path = os.path.join("scripts","rendering_batteries","batch_03",config_name+".json") 
            # params_exploration(config_path)
            # gen_table(config_path=config_path)

main()




# with open(file_path, 'w') as file:
    # file.writelines(lines)
# Save the edited content back to the file
