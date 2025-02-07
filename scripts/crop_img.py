from PIL import Image
import os,glob
import argparse

# def crop_image(input_image_path, output_image_path, crop_area):
#     """
#     Crop an image and save it.

#     :param input_image_path: Path to the input image.
#     :param output_image_path: Path to save the cropped image.
#     :param crop_area: Tuple specifying the crop rectangle (left, upper, right, lower).
#     """
#     print(input_image_path)
#     try:
#         with Image.open(input_image_path) as img:
#             cropped_img = img.crop(crop_area)
#             cropped_img.save(output_image_path)
#             print(f"Cropped image saved to {output_image_path}")
#     except Exception as e:
#         print(f"An error occurred: {e}")
# if __name__ == "__main__":
    
#     parser = argparse.ArgumentParser(description="Fitting model to measurements data")
#     parser.add_argument("input_dir",  help="folder where the pngs file that needs to be cropped are. This dir is where the cropped image will be saved")
#     parser.add_argument("crop_area", type=tuple help="folder where the pngs file that needs to be cropped are. This dir is where the cropped image will be saved")
    
#     args = parser.parse_args()

#     # Example crop area (left, upper, right, lower)
#     # Modify these values as needed
#     crop_area = (100, 100, 200, 200)
#     input_dir = args.input_dir
#     png_files = glob.glob(os.path.join(input_dir, "*.png"))
#     for image in png_files:
#         basename = os.path.splitext(os.path.basename(png_files))[0]
#         output_path = os.path.join(input_dir,basename + "_cropped.png")
#         crop_image(image, output_path, crop_area)



# import os
# import glob
# import argparse
# from PIL import Image

# def collect_png_files(directory):
#     """
#     Collect all PNG files in the specified directory and return their base names without extensions.

#     :param directory: Path to the directory to search for PNG files.
#     :return: List of paths to PNG files.
#     """
#     return glob.glob(os.path.join(directory, "*.png"))

# def crop_image(input_image_path, output_image_path, crop_area):
#     """
#     Crop an image and save it.

#     :param input_image_path: Path to the input image.
#     :param output_image_path: Path to save the cropped image.
#     :param crop_area: Tuple specifying the crop rectangle (left, upper, right, lower).
#     """
#     try:
#         with Image.open(input_image_path) as img:
#             cropped_img = img.crop(crop_area)
#             cropped_img.save(output_image_path)
#             print(f"Cropped image saved to {output_image_path}")
#     except Exception as e:
#         print(f"An error occurred: {e}")

def crop_image(input_image_path, output_image_path, crop_area):
    """
    Crop an image and save it.

    :param input_image_path: Path to the input image.
    :param output_image_path: Path to save the cropped image.
    :param crop_area: Tuple specifying the crop rectangle (left, upper, right, lower).
    """
    try:
        with Image.open(input_image_path) as img:
            cropped_img = img.crop(crop_area)
            cropped_img.save(output_image_path)
            print(f"Cropped image saved to {output_image_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

def collect_png_files(directory):
    """
    Collect all PNG files in the specified directory and return their base names without extensions.

    :param directory: Path to the directory to search for PNG files.
    :return: List of base names of PNG files without extensions.
    """
    png_files = glob.glob(os.path.join(directory, "*.png"))
    return png_files

def tuple_type(strings):
    strings = strings.replace("(", "").replace(")", "")
    mapped_int = map(int, strings.split(","))
    return tuple(mapped_int)

def parse_crop_area(crop_area_str):
    """
    Parse a crop area string to a tuple.

    :param crop_area_str: String specifying the crop area.
    :return: Tuple specifying the crop area (left, upper, right, lower).
    """
    try:
        values = tuple(map(int, crop_area_str.strip('()').split(',')))
        if len(values) != 4:
            raise ValueError
        return values
    except:
        raise argparse.ArgumentTypeError("Crop area must be a tuple of four integers (left, upper, right, lower)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collect PNG files and optionally crop them.")
    parser.add_argument('directory', type=str, help="Path to the directory to search for PNG files.")
    parser.add_argument('config', type=str, help="Based on the config type different crop regions will generated. Options are: Emily_frontal, Emily_side, custom")
    
    # parser.add_argument('--crop', type=parse_crop_area,  help="Crop area as a tuple of four integers (left, upper, right, lower).")
    #frontal emily (225,333,671,565)
    #sss emily (612.181,1057,412)
    args = parser.parse_args()
    directory = args.directory
    import json

    # Specify the path to the JSON file
    t = os.path.dirname(__file__)
    file_path = os.path.join(t,"crop_img_settings.json")
    # Open and load the JSON file
    with open(file_path, 'r') as file:
        data = json.load(file)
    crop_area = data

    if not os.path.isdir(directory):
        print(f"The path {directory} is not a valid directory.")
    else:
        png_files = collect_png_files(directory)
        
        if png_files:
            print("PNG files found in the directory:")
            for file in png_files:            
                print(file)
                if crop_area:
                    base_name = os.path.splitext(os.path.basename(file))[0]
                    output_path = os.path.join(directory, f"{base_name}_cropped.png")
                    crop_image(file, output_path, crop_area)
        else:
            print("No PNG files found in the directory.")
