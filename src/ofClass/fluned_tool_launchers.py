import os
import subprocess


def launch_volume_func_object(path):
    """
    this function launch the utilities that calculates the volumes
    - this field is needed only for the activation rate
    interpolation
    """

    print("calculating volumes  ...")
    orig_folder = os.getcwd()
    os.chdir(path)
    launch_volumes = "postProcess -func writeCellVolumes".split()
    with open("log", "a", encoding="utf-8") as outfile:
        subprocess.Popen(launch_volumes, stdout=outfile).wait()
    os.chdir(orig_folder)

    return None


def launch_centroid_func_object(path):
    """this function launch the utilities that calculates the volume and
    the centroids - this fields are needed only for the activation rate
    interpolation"""

    print("calculating centroids ...")
    origFolder = os.getcwd()
    os.chdir(path)
    launchCentroids = "postProcess -func writeCellCentres".split()
    with open("log", "a") as outfile:
        _ = subprocess.Popen(launchCentroids, stdout=outfile).wait()
    os.chdir(origFolder)

    return


def launch_grad_func_object(path):
    """
    launches the utility to calculate T gradient
    """
    print("launching gradient function object")
    orig_folder = os.getcwd()
    os.chdir(path)
    launch_grad = "postProcess -func 'grad(T)'".split()
    with open("log", "a", encoding="utf-8") as out_file:
        subprocess.Popen(launch_grad, stdout=out_file).wait()
    os.chdir(orig_folder)

    return


def generate_vtk(path):
    """
    launches the utility to create a vtk file
    """

    print("launching FoamToVTK utility")
    orig_folder = os.getcwd()
    os.chdir(path)
    cmd_str_1 = "foamToVTK -latestTime -noFaceZones -noFunctionObjects "
    cmd_str_3 = ' -excludePatches (".*")'
    launch_f2vtk = cmd_str_1.split()
    launch_f2vtk.append("-fields")
    launch_f2vtk.append("(T Ta Td)")
    launch_f2vtk.extend(cmd_str_3.split())
    # print (launch_f2vtk)
    with open("log", "a", encoding="utf-8") as outfile:
        subprocess.Popen(launch_f2vtk, stdout=outfile).wait()
    os.chdir(orig_folder)

    return
