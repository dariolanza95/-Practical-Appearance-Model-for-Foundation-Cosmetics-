import os
cosmetics = ["Clarins112","Clarins108","Clarins103","Clarins105"]
# cosmetics = ["Clarins108"]

tasks = ["full"]
# tasks = ["task01","task0"]
# tasks = ["task01","task02"]
HQ = False
for cosmetic in cosmetics:
    for task in tasks:
        res = 128
        spp_hq = 8192
        if task == "full":
            if HQ:
                res = spp_hq
                command = "python scripts/measurements/fitting_measurements_cma_rgb_patches_backreflectance.py -mode full -isHQ -optimizer \"Powell\" " + \
                "-cosmetic_type {0} -spp {1} -parameter_bounds parameter_bound/parameter_bounds.json -output_folder resultsOursFinal04/{0} -initial_parameters resultsOursFinal04/{0}/best_result.json".format(cosmetic,res)
            else:
                command = "python scripts/measurements/fitting_measurements_cma_rgb_patches_backreflectance.py -mode full -optimizer \"Powell\" " + \
                "-cosmetic_type {0} -spp {1} -parameter_bounds parameter_bound/parameter_bounds.json -output_folder resultsOursFinal06/{0} -parameter_bounds parameter_bound/parameter_bounds.json -initial_parameters optimized_foundations/{0}/best_result.json".format(cosmetic,res)
                # "-cosmetic_type {0} -spp {1} -parameter_bounds parameter_bound/parameter_bounds.json -output_folder resultsOursFinal06/{0} -initial_parameters resultsOurs05/{0}/best_result.json".format(cosmetic,res)
            os.system(command)
        elif task == "task01":
            if HQ:
                spp = spp_hq
                command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -optimizer \"Powell\" -mode HQ " + \
            "-cosmetic_type {0} -spp {1} -output_folder resultsAblationT01/{0}/HQ -initial_parameters resultsAblationT01/{0}/best_result.json".format(cosmetic,res)
            else:
                command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -optimizer \"Powell\" -mode task01 " + \
            "-cosmetic_type {0} -spp {1} -output_folder resultsAblationT01/{0} -initial_parameters resultsOursFinal0/{0}/best_result.json -parameter_bounds parameter_bound/parameter_bounds.json".format(cosmetic,res)
            os.system(command)

        elif task == "task02":
            if HQ:
                spp = spp_hq
                command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -optimizer \"Powell\" -mode HQ " + \
            "-cosmetic_type {0} -spp {1} -output_folder resultsAblationT01/{0}/HQ -initial_parameters resultsAblationT02/{0}/best_result.json".format(cosmetic,res)
            else:
                command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -optimizer \"Powell\" -mode task02 " + \
            "-cosmetic_type {0} -spp {1} -output_folder resultsAblationT02/{0} -initial_parameters optimized_foundations/{0}/best_result.json -parameter_bounds parameter_bound/parameter_bounds.json".format(cosmetic,res)
            os.system(command)
        elif task == "task03":
            if HQ:
                spp = spp_hq
                command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -optimizer \"Powell\" -mode HQ " + \
            "-cosmetic_type {0} -spp {1} -output_folder resultsAblationT01/{0}/HQ -initial_parameters resultsAblationT03/{0}/best_result.json".format(cosmetic,res)
            else:
                command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -optimizer \"Powell\" -mode task03 " + \
            "-cosmetic_type {0} -spp {1} -output_folder resultsAblationT03/{0} -initial_parameters optimized_foundations/{0}/best_result.json -parameter_bounds parameter_bound/parameter_bounds.json".format(cosmetic,res)
            os.system(command)
        elif task == "task04":
            if HQ:
                spp = 2048
                command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -optimizer \"Powell\" -mode HQ " + \
            "-cosmetic_type {0} -spp {1} -output_folder resultsAblationT01/{0}/HQ -initial_parameters resultsAblationT04/{0}/best_result.json".format(cosmetic,res)
            else:
                command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -optimizer \"Powell\" -mode task04 " + \
            "-cosmetic_type {0} -spp {1} -output_folder resultsAblationT04/{0} -initial_parameters optimized_foundations/{0}/best_result.json -parameter_bounds parameter_bound/parameter_bounds_T04.json".format(cosmetic,res)
            os.system(command)
        elif task == "color":
            tasks = ["01","02","03"]
            for task in tasks:
                command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -cosmetic_type {0} -optimizer Powell -mode color -output_folder resultsAblationT{1}/{0}/ -initial_parameters resultsAblationT{1}/{0}/best_result.json -spp {2} \
                    -parameter_bounds parameter_bound/parameter_bounds.json".format(cosmetic,task,res)
                os.system(command)

            # command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -mode color02 -cosmetic_type {0} -optimizer Powell -output_folder resultsOursFinal/{0} -initial_parameters resultsOurs/{0}/log/best_loss_000001.json -spp {1} -parameter_bounds resultsOurs/Clarins105/parameter_bounds.json".format(cosmetic,res)
            os.system(command)
        elif task == "rescale":
            tasks = ["01","02","03"]
            for task in tasks:
                command = "python scripts/measurements/replot_fitting_data.py -cosmetic_type {0} -output_folder resultsAblationT{1}/{0}/log/ ".format(cosmetic,task)
                os.system(command)
    


#Trash
"""
if task == "ours-full":
            command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -mode full -optimizer \"Powell\" " + \
                "-cosmetic_type {0} -spp {1} -output_folder resultsOursFinal02/{0} -initial_parameters resultsOurs/{0}/log/best_loss_000001.json".format(cosmetic,res)
            os.system(command)
        elif task == "task02":
            command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -optimizer \"Powell\" -mode task02 " + \
            "-cosmetic_type {0} -spp {1} -output_folder resultsAblationT02/{0} -initial_parameters AblationTask02/{0}/log/best_loss_000001.json".format(cosmetic,res)
            os.system(command)
        elif task == "task03":
            command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -optimizer \"Powell\" -mode task03 " + \
            "-cosmetic_type {0} -spp {1} -output_folder resultsAblationT03/{0} -initial_parameters AblationTask03/{0}/log/best_loss_000001.json".format(cosmetic,res)
            os.system(command)
        elif task == "task04":
            command = "python scripts/measurements/fitting_measurements_cma_rgb_patches.py -optimizer \"Powell\" -mode task04 " + \
            "-cosmetic_type {0} -spp {1} -output_folder resultsAblationT04/{0} -initial_parameters AblationTask04/{0}/log/best_loss_000001.json".format(cosmetic,res)
            os.system(command)
"""