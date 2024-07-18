from PIL import Image
import os
import glob
import numpy as np
from paraview.simple import *

def white_base(file_path, n):
    """ 
    This function replaces the last n lines of an image with white lines. 
    In the context of a granuheap photo, this function replaces the granuheap base with white lines. 
    
    file_path (string) : path of the file to modify (e.g. '/home/user/work/.../your_file.png')
    n (int) : number of lines to change white.
    """

    image = Image.open(file_path).convert('L')
    width, height = image.size
    if n > height:
        raise ValueError("The number of lines to be deleted is greater than the image height.")
    new_height = height - n
    new_image = Image.new('L', (width, height), color=255) 
    for x in range(width):
        for y in range(new_height):
            pixel = image.getpixel((x, y))
            new_image.putpixel((x, y), pixel)
    new_image.save(file_path)

    print( f"white_base : {n} lines was change to white from the bottom of {file_path}. ")
    return 

def scaling_paraview(width_exp, pvd_path, pvd_name, out_path, scale):
    """
    This function is used to scale simulation images to the same scale as experimental photos. 
    This function verify if the number of pixel in width of both image is the same.
    The function returns whether the images are at the same scale (True/False), as well as the height of the support in the simulations (float).

    width_exp (float) : width of the experiemental support.
    pvd_path (string) : path to your simulation directory (e.g. '/home/user/work/.../yourdirectory/')
    pvd_name (string) : name of your simulation (e.g. 'granuheap')
    out_path (string) : name of your output directory (e.g. 'output')
    scale (float) : scale number in paraview macro
    """

    path_input_file = pvd_path + out_path + '/'
    filename1 = pvd_name + '.pvd'
    filename2 = pvd_name + '.solid_surface.01.pvd'
    nb_image = 32
    ImageResolution =[400, 600]
    path_output_file = pvd_path + 'map/'

    if not os.path.exists(path_output_file):
        # Créer le dossier
        os.makedirs(path_output_file)
        print(f"Directory '{path_output_file}' created.")


    ############################
    #   PARAVIEW MACRO START   #
    ############################

    #### import the simple module from the paraview
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'PVD Reader'
    outpvd = PVDReader(registrationName=filename1, FileName=path_input_file+filename1)

    # get animation scene
    animationScene1 = GetAnimationScene()
    timeKeeper1 = GetTimeKeeper()
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # create a new 'PVD Reader'
    outsolid_object01pvd = PVDReader(registrationName=filename2, FileName=path_input_file+filename2)

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    outpvdDisplay = Show(outpvd, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    outpvdDisplay.Representation = 'Surface'
    outpvdDisplay.ColorArrayName = [None, '']
    outpvdDisplay.SelectTCoordArray = 'None'
    outpvdDisplay.SelectNormalArray = 'None'
    outpvdDisplay.SelectTangentArray = 'None'
    outpvdDisplay.OSPRayScaleArray = 'ID'
    outpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    outpvdDisplay.SelectOrientationVectors = 'None'
    outpvdDisplay.ScaleFactor = 0.0006039923056960106
    outpvdDisplay.SelectScaleArray = 'ID'
    outpvdDisplay.GlyphType = 'Arrow'
    outpvdDisplay.GlyphTableIndexArray = 'ID'
    outpvdDisplay.GaussianRadius = 3.019961528480053e-05
    outpvdDisplay.SetScaleArray = ['POINTS', 'ID']
    outpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    outpvdDisplay.OpacityArray = ['POINTS', 'ID']
    outpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    outpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    outpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    outpvdDisplay.ScalarOpacityUnitDistance = 0.0005043115016812344
    outpvdDisplay.OpacityArrayName = ['POINTS', 'ID']
    outpvdDisplay.SelectInputVectors = ['POINTS', 'fem_force']
    outpvdDisplay.WriteLog = ''
    outpvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 7055.0, 1.0, 0.5, 0.0]
    outpvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 7055.0, 1.0, 0.5, 0.0]

    # reset view to fit data
    renderView1.ResetCamera(False)

    # show data in view
    outsolid_object01pvdDisplay = Show(outsolid_object01pvd, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    outsolid_object01pvdDisplay.Representation = 'Surface'
    outsolid_object01pvdDisplay.ColorArrayName = [None, '']
    outsolid_object01pvdDisplay.SelectTCoordArray = 'None'
    outsolid_object01pvdDisplay.SelectNormalArray = 'None'
    outsolid_object01pvdDisplay.SelectTangentArray = 'None'
    outsolid_object01pvdDisplay.OSPRayScaleArray = 'displacement'
    outsolid_object01pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    outsolid_object01pvdDisplay.SelectOrientationVectors = 'None'
    outsolid_object01pvdDisplay.ScaleFactor = 0.0009994862601161003
    outsolid_object01pvdDisplay.SelectScaleArray = 'None'
    outsolid_object01pvdDisplay.GlyphType = 'Arrow'
    outsolid_object01pvdDisplay.GlyphTableIndexArray = 'None'
    outsolid_object01pvdDisplay.GaussianRadius = 4.997431300580502e-05
    outsolid_object01pvdDisplay.SetScaleArray = ['POINTS', 'displacement']
    outsolid_object01pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    outsolid_object01pvdDisplay.OpacityArray = ['POINTS', 'displacement']
    outsolid_object01pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    outsolid_object01pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    outsolid_object01pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    outsolid_object01pvdDisplay.ScalarOpacityUnitDistance = 0.001127148890411883
    outsolid_object01pvdDisplay.OpacityArrayName = ['POINTS', 'displacement']
    outsolid_object01pvdDisplay.SelectInputVectors = ['POINTS', 'displacement']
    outsolid_object01pvdDisplay.WriteLog = ''
    outsolid_object01pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
    outsolid_object01pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source to outpvd
    SetActiveSource(outpvd)
    outpvdDisplay.SetRepresentationType('Point Gaussian')
    outpvdDisplay.GaussianRadius = 0.5
    outpvdDisplay.ScaleByArray = 1
    outpvdDisplay.SetScaleArray = ['POINTS', 'diameter']
    outpvdDisplay.UseScaleFunction = 0
    outpvdDisplay.DiffuseColor = [0.0, 0.0, 0.0]

    # set active source to outsolid_object01pvd
    SetActiveSource(outsolid_object01pvd)
    outsolid_object01pvdDisplay.DiffuseColor = [0.0, 0.0, 0.0]

    # scene manipulation
    animationScene1.GoToLast()
    renderView1.ResetActiveCameraToPositiveZ()
    renderView1.AdjustRoll(-90.0)
    renderView1.InteractionMode = '2D'
    renderView1.ResetCamera(True)

    layout1 = GetLayout()
    layout1.SetSize(400, 600)

    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [-0.0060681740287691355, -2.51084566116333e-06, 0.029958768572906826]
    renderView1.CameraFocalPoint = [-0.0060681740287691355, -2.51084566116333e-06, 1.2052478268742561e-05]
    renderView1.CameraViewUp = [1.0, 2.220446049250313e-16, 0.0]
    renderView1.OrientationAxesVisibility = 0
    renderView1.CameraParallelScale = scale
    Hide(outpvd, renderView1)

    SaveScreenshot(path_output_file + "1.png", renderView1, ImageResolution=ImageResolution)

    Delete(outpvd)
    del outpvd
    ############################
    #    PARAVIEW MACRO END    #
    ############################

    img = Image.open(path_output_file + "1.png").convert('L')
    pixels = img.load()

    # Get the width of the support 
    width = 0
    j=img.size[1]-1
    while width ==0:
        for i in range(img.size[0]):
            if pixels[i,j] < 128:
                width += 1
        j -= 1 
        if j<0:
        	raise ValueError("White image")

    #Get the height of the support 
    height = 0
    for i in range(img.size[0]):
        black_px = 0
        for j in range(img.size[1]):
        	if pixels[i, j] < 128:
        		black_px += 1
        if black_px > height:
        	height = black_px
        	
    if width == width_exp:
        to_print = f"scaling_paraview : Sucess - width = { str(width)} black px, Height = {height} px, Scale = {scale}"
        return True, height, width, to_print
    else:
        to_print = f"scaling_paraview : Failed - width = { str(width)} black px, Expected width = {width_exp} px, Scale = {scale} "
        return False, height, width, to_print

def picture_paraview(pvd_path, pvd_name, out_path, scale):
    """
    This function is used to take 16 pictures of the simulation at the final state following a 180 arc around the heap. 
    In the context of a granuheap simulation, this function imitate the process of taking picture at the end of granuheap experiment.

    pvd_path (string) : path to your simulation directory (e.g. '/home/user/work/.../yourdirectory/')
    pvd_name (string) : name of your simulation (e.g. 'granuheap')
    out_path (string) : name of your output directory (e.g. 'output')
    scale (float) : scale number in paraview macro
    """

    path_input_file = pvd_path + out_path + '/'
    filename1 = pvd_name + '.pvd'
    filename2 = pvd_name + '.solid_surface.01.pvd'
    nb_image = 32
    ImageResolution =[400, 600]
    path_output_file = pvd_path + 'map/'

    if not os.path.exists(path_output_file):
        # Créer le dossier
        os.makedirs(path_output_file)
        print(f"Dossier '{path_output_file}' créé.")

    ############################
    #   PARAVIEW MACRO START   #
    ############################

    #### import the simple module from the paraview
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'PVD Reader'
    outpvd = PVDReader(registrationName=filename1, FileName=path_input_file+filename1)

    # get animation scene
    animationScene1 = GetAnimationScene()
    timeKeeper1 = GetTimeKeeper()
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # create a new 'PVD Reader'
    outsolid_object01pvd = PVDReader(registrationName=filename2, FileName=path_input_file+filename2)

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    outpvdDisplay = Show(outpvd, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    outpvdDisplay.Representation = 'Surface'
    outpvdDisplay.ColorArrayName = [None, '']
    outpvdDisplay.SelectTCoordArray = 'None'
    outpvdDisplay.SelectNormalArray = 'None'
    outpvdDisplay.SelectTangentArray = 'None'
    outpvdDisplay.OSPRayScaleArray = 'ID'
    outpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    outpvdDisplay.SelectOrientationVectors = 'None'
    outpvdDisplay.ScaleFactor = 0.0006039923056960106
    outpvdDisplay.SelectScaleArray = 'ID'
    outpvdDisplay.GlyphType = 'Arrow'
    outpvdDisplay.GlyphTableIndexArray = 'ID'
    outpvdDisplay.GaussianRadius = 3.019961528480053e-05
    outpvdDisplay.SetScaleArray = ['POINTS', 'ID']
    outpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    outpvdDisplay.OpacityArray = ['POINTS', 'ID']
    outpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    outpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    outpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    outpvdDisplay.ScalarOpacityUnitDistance = 0.0005043115016812344
    outpvdDisplay.OpacityArrayName = ['POINTS', 'ID']
    outpvdDisplay.SelectInputVectors = ['POINTS', 'fem_force']
    outpvdDisplay.WriteLog = ''
    outpvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 7055.0, 1.0, 0.5, 0.0]
    outpvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 7055.0, 1.0, 0.5, 0.0]

    # reset view to fit data
    renderView1.ResetCamera(False)

    # show data in view
    outsolid_object01pvdDisplay = Show(outsolid_object01pvd, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    outsolid_object01pvdDisplay.Representation = 'Surface'
    outsolid_object01pvdDisplay.ColorArrayName = [None, '']
    outsolid_object01pvdDisplay.SelectTCoordArray = 'None'
    outsolid_object01pvdDisplay.SelectNormalArray = 'None'
    outsolid_object01pvdDisplay.SelectTangentArray = 'None'
    outsolid_object01pvdDisplay.OSPRayScaleArray = 'displacement'
    outsolid_object01pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    outsolid_object01pvdDisplay.SelectOrientationVectors = 'None'
    outsolid_object01pvdDisplay.ScaleFactor = 0.0009994862601161003
    outsolid_object01pvdDisplay.SelectScaleArray = 'None'
    outsolid_object01pvdDisplay.GlyphType = 'Arrow'
    outsolid_object01pvdDisplay.GlyphTableIndexArray = 'None'
    outsolid_object01pvdDisplay.GaussianRadius = 4.997431300580502e-05
    outsolid_object01pvdDisplay.SetScaleArray = ['POINTS', 'displacement']
    outsolid_object01pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    outsolid_object01pvdDisplay.OpacityArray = ['POINTS', 'displacement']
    outsolid_object01pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    outsolid_object01pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    outsolid_object01pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    outsolid_object01pvdDisplay.ScalarOpacityUnitDistance = 0.001127148890411883
    outsolid_object01pvdDisplay.OpacityArrayName = ['POINTS', 'displacement']
    outsolid_object01pvdDisplay.SelectInputVectors = ['POINTS', 'displacement']
    outsolid_object01pvdDisplay.WriteLog = ''
    outsolid_object01pvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
    outsolid_object01pvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source to outpvd
    SetActiveSource(outpvd)
    outpvdDisplay.SetRepresentationType('Point Gaussian')
    outpvdDisplay.GaussianRadius = 0.5
    outpvdDisplay.ScaleByArray = 1
    outpvdDisplay.SetScaleArray = ['POINTS', 'diameter']
    outpvdDisplay.UseScaleFunction = 0
    outpvdDisplay.DiffuseColor = [0.0, 0.0, 0.0]

    # set active source to outsolid_object01pvd
    SetActiveSource(outsolid_object01pvd)
    outsolid_object01pvdDisplay.DiffuseColor = [0.0, 0.0, 0.0]

    # scene manipulation
    animationScene1.GoToLast()
    renderView1.ResetActiveCameraToPositiveZ()
    renderView1.AdjustRoll(-90.0)
    renderView1.InteractionMode = '2D'
    renderView1.ResetCamera(True)

    layout1 = GetLayout()
    layout1.SetSize(400, 600)

    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [-0.0060681740287691355, -2.51084566116333e-06, 0.029958768572906826]
    renderView1.CameraFocalPoint = [-0.0060681740287691355, -2.51084566116333e-06, 1.2052478268742561e-05]
    renderView1.CameraViewUp = [1.0, 2.220446049250313e-16, 0.0]
    renderView1.OrientationAxesVisibility = 0
    renderView1.CameraParallelScale = scale

    SaveScreenshot(path_output_file + "1.png", renderView1, ImageResolution=ImageResolution)
    img_sim = Image.open(path_output_file + "1.png").convert('L')

    for i in range(32):
        renderView1.AdjustAzimuth(180/32)
        SaveScreenshot(path_output_file + str(i+1) + ".png", renderView1, ImageResolution=ImageResolution)
    
    Delete(outpvd)
    del outpvd

    ############################
    #    PARAVIEW MACRO END    #
    ############################

    print(f"picture_paraview : 32 images taken - {pvd_path}")

def map_processing(input_files_list, output_path, height, height_exp):
    """
    This function take a list of image paths and create a single image that reflects the relative presence of black across the set of input image. 

    input_files_list (array) :  list of images' path that will create the output image. (use the find_files function)
    output_path (string) : path where will be saved the output image (e.g. '/map.png')
    height (int) : number of pixels in height of the numerical support (use the scaling_paraview function)
    height_exp (int) : number of pixels in height of the experimental support
    """

    img = Image.open(input_files_list[0]).convert('L')
    pixels_map = np.full((img.size[0], img.size[1]), 255)
    nb_image = len(input_files_list)
    div_nb_image = 1/nb_image

    # map processing
    for i in range(nb_image):
        img = Image.open(input_files_list[i]).convert('L')
        pixels = img.load()
        for j in range(img.size[1]):
            for i in range(img.size[0]):
                if pixels[i,j] < 100:
                    pixels_map[i,j] -= (255*div_nb_image)
    
    # counting the number of white lines at the bottom
    white_line = True
    nb_white_lines = 0
    for j in reversed(range(img.size[1])):
        for i in reversed(range(img.size[0])):
            if pixels[i,j] < 100:
                white_line = False
        if white_line:
            nb_white_lines += 1 

    # adjust the position of the heap in the picture to match the experimental resutls by adding and deleting lines.
    rows_to_delete = nb_white_lines - height_exp + height
    rows_to_delete_list = np.arange(img.size[1]-rows_to_delete, img.size[1])
    pixels_map = np.delete(pixels_map, rows_to_delete_list, axis=1)
    row_to_add = np.full((img.size[0], rows_to_delete), 255)
    pixels_map = np.hstack((row_to_add, pixels_map))

        
    new_image_pixels = np.array(pixels_map, dtype=np.uint8)
    new_image = Image.fromarray(new_image_pixels)
    new_image = new_image.rotate(-90, expand=True)
    new_image.save(output_path)
    white_base(output_path, height_exp)
    
    print( f"map_processing : Map of {nb_image} images saved at {output_path}")
    return output_path


def find_files(root_dir):
    """
    This function return all the files in a directory and return a list of them (array).

    root_dir (string) : path to the directory where the images are saved
    """
    found_files = []
    
    # Parcourir le répertoire racine et ses sous-répertoires
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            # Construire le chemin absolu du fichier trouvé
            file_path = os.path.join(dirpath, filename)
            # Ajouter le chemin absolu à la liste des fichiers trouvés
            found_files.append(file_path)
    
    print(f'find_files_with_name : {len(found_files)} paths returned in a list.')
    return found_files


def difference_images(image_path1, image_path2, output_path, border_size=5):
    """
    This function take two gray scaled images and create a new image by substracting image 2 from image 1. 
    If the results is positive, bleu is use. If the results is negative, red is use. 
    The image is cropped and saved

    image_path1 (string) : path to the first image (e.g. 'image1.png')
    image_path2 (string) : path to the second image (e.g. 'image2.png')
    output_path (string) : path where the resulting image will be saved. (e.g. 'image_difference.png')
    border_size (int) : border size in pixels for the cropped image. 

    """
    image1 = Image.open(image_path1).convert('RGB')
    pixels1 = image1.load()
    image2 = Image.open(image_path2).convert('RGB')
    pixels2 = image2.load()

    # Create the new image 
    pixels_map = np.zeros((image1.size[0], image1.size[1], 3), dtype=np.uint8)
    for i in range(image1.size[0]):
        for j in range(image1.size[1]):
            diff = pixels1[i,j][0] - pixels2[i,j][0]
            if diff >= 0 :
                pixels_map[i, j] = (255, 255 - diff, 255 - diff)
            elif diff < 0 :
                pixels_map[i, j] = (255 - abs(diff), 255 - abs(diff), 255)
    
    # Crop image 
    width, height = image1.size

    left_limit = 0
    right_limit = width - 1
    for x in range(width):
        if any(pixels_map[x, y, 0] != 255 for y in range(height)):
            left_limit = max(0, x - border_size)
            break
    for x in range(width - 1, -1, -1):
        if any(pixels_map[x, y, 0] != 255 for y in range(height)):
            right_limit = min(width - 1, x + border_size)
            break
    
    top_limit = 0
    bottom_limit = height - 1
    for y in range(height):
        if any(pixels_map[:, y, 0] != 255):
            top_limit = max(0, y - border_size)
            break
    for y in range(height - 1, -1, -1):
        if any(pixels_map[:, y, 0] != 255):
            bottom_limit = min(height - 1, y + border_size)
            break

    cropped_pixels_map = pixels_map[left_limit:right_limit+1, top_limit:bottom_limit+1, :]
    
    # save the new image
    new_image = np.array(cropped_pixels_map, dtype=np.uint8)
    new_image = Image.fromarray(new_image)
    new_image = new_image.rotate(-90, expand=True)
    new_image.save(output_path)

    print(f" difference_image : the subtraction between {image_path1} and {image_path2} have been saved at {output_path}")
    return 

def RMSE(image_exp, image_num):
    """
    This function calculate the root mean squared error of a numerical image. 

    image_exp (string) : path to the experimental result (e.g. 'experimental_results.png')
    image_num (string) : path to the numerical result (e.g. 'map.png')
    """
    image_exp = Image.open(image_exp).convert('L')
    pixels_exp = image_exp.load()
    image_num = Image.open(image_num).convert('L')
    pixels_num = image_num.load()
    sum = 0
    for i in range(image_exp.size[0]):
        for j in range(image_exp.size[1]):
            sum += abs(pixels_exp[i,j] - pixels_num[i,j])**2
    
    rmse = np.sqrt(sum/(image_exp.size[0]*image_exp.size[1]))

    print(f" RMSE : the Root Mean Squared Error  is {rmse}.")
    return rmse

