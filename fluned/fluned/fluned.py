# -*- coding: utf-8 -*-
import re
import sys
import os
import vtk
from vtk.util import numpy_support as VN
import pathlib
import shutil
import numpy as np
import subprocess
import argparse
import h5py
from fluned.fluned_h5_utils import get_dataset_keys
from fluned.fluned_h5_utils import get_h5_dataset
from fluned.fluned_h5_utils import get_h5_dataset_multi
from fluned.fluned_h5_utils import get_h5_path_dataset
from fluned.fluned_h5_utils import extract_multiblock
from fluned.fluned_bin_utils import get_fluent_binarray_double
from fluned.fluned_bin_utils import get_fluent_parse_headers
from fluned.fluned_bin_utils import get_fluent_parse_regions
from fluned.of_class import SimulationOF



def create_input_template():
    """
    this function generates an input template file
    """

    input_name = "inputTemplate"
    current_folder = os.getcwd()

    input_path = os.path.join(current_folder,input_name)

    template_text = """CASE  FLUNED_01_DEFAULT_N16
TIME_TREATMENT  steadyState #steadyState or transient supported
#ACTIVATION_FILE  {}
#ACTIVATION_DATASET    "Value - Total"
#ACTIVATION_CONST 1e16
#ACTIVATION_NORMALIZATION 0 #leave zero if no normalization is required
#DECAY_CONSTANT 0.1661825   #N17 decay const
DECAY_CONSTANT 0.09721559   #N16 decay const
INLET_CONC 1e10
MOLECULAR_DIFFUSION 2e-09
SCHMIDT_NUMBER   0.7
CFD_PATH      "{}"
CFD_TYPE    OpenFoam       # OpenFoam, fluent-h5-multi types supported
#FLUENT_FLUID_REGION_NAME     region_name
"""

    with open(input_path,'w',encoding='utf-8') as fw:
        fw.write(template_text.format(current_folder,current_folder))

    return 0


def sample_coordinates_VTK(vtk_file, dataset_name, coordinates):
    """
    this function reads the vtk and sample the reaction rates
    """

    checkFile = os.path.isfile(vtk_file)
    if not checkFile:
        print ("ERROR activation file not found")
        sys.exit()


    #read the vtk file with an unstructured grid
    reader = vtk.vtkStructuredGridReader()
    reader.SetFileName(vtk_file)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()



    #define probe

    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(coordinates))

    #print ("parsed  0/{:d}".format(len(coordinates)))
    #step = pow(10,int(math.log10(len(coordinates)))-1)

    for i,val in enumerate(coordinates):
        points.SetPoint(i,val[0],val[1],val[2])
   #     if (i+1) % step == 0:
   #         print ("parsed  {:d}/{:d}".format(i+1,len(coordinates)))



    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)


    #Perform the interpolation
    probeFilter =vtk.vtkProbeFilter()
    probeFilter.SetSourceData(data)
    probeFilter.SetInputData(polydata)
    probeFilter.Update()

    vtkArray = probeFilter.GetOutput().GetPointData().GetArray(dataset_name)

    reacRates = VN.vtk_to_numpy(vtkArray)

    return  reacRates


#def sampleCoordinatesVTK1(vtkFile, datasetName, coordinates):
#    """this function reads the vtk and sample the reaction rates"""
#
#    checkFile = os.path.isfile(vtkFile)
#    if not checkFile:
#        print ("ERROR activation file not found")
#        sys.exit()
#    print ("Sampling Reaction Rate file... points ")
#
#    #read the vtk file with an unstructured grid
#    reader = vtk.vtkStructuredGridReader()
#    reader.SetFileName(vtkFile)
#    reader.ReadAllVectorsOn()
#    reader.ReadAllScalarsOn()
#    reader.Update()
#    data = reader.GetOutput()
#
#
#    bounds =list(data.GetBounds())
#    dims =list(data.GetDimensions())
#    print (data)
#    print (bounds)
#    print (bounds)
#
#
#
#    #define probe
#
#    points = vtk.vtkPoints()
#    points.SetNumberOfPoints(len(coordinates))
#
#    for i,val in enumerate(coordinates):
#        points.SetPoint(i,val[0],val[1],val[2])
#
#
#
#    polydata = vtk.vtkPolyData()
#    polydata.SetPoints(points)
#
#
#    #Perform the interpolation
#    interpolator = vtk.vtkPointInterpolator()
#    #interpolator.SetInputData(polydata)
#    interpolator.SetInputConnection(polydata.GetPointData().GetOutputPort())
#    #probeFilter.SetSourceData(data)
#    #probeFilter.SetInputData(polydata)
#    #probeFilter.Update()
#    #vtkArray = probeFilter.GetOutput().GetPointData().GetArray(datasetName)
#
#    #reacRates = VN.vtk_to_numpy(vtkArray)
#
#    interpolator.Update()
#    sys.exit()
#
#
#    sys.exit()
#
#
#    return  reacRates

def heronFormula(points):

    dist1 = point_distance(points[0],points[1])
    dist2 = point_distance(points[1],points[2])
    dist3 = point_distance(points[2],points[0])

    semip = (dist1+dist2+dist3)/2

    area = pow((semip*(semip-dist1)*(semip-dist2)*(semip-dist3)),0.5)
    return area

def calculateEdgeCentroid(vector):

    sumX = 0
    sumY = 0
    sumZ = 0

    for vertex in vector:
        sumX += vertex[0]
        sumY += vertex[1]
        sumZ += vertex[2]

    nPoints = len(vector)

    centroid = [sumX/nPoints, sumY/nPoints, sumZ/nPoints]

    return centroid

def point_distance(point1,point2):
    """
    this function calculates the distance between two points in 3D and 2D

    Args:
        point1 (list): coordinates of the first point
        point2 (list): coordinates of the second point

    """

    dist = pow((sum([pow((point1[i] - point2[i]),2) for i in [0,1,2]])),0.5)

    return dist

def calculateArea(points):
    """
    this function calculate the area of the element surface - these can
    be tetra with a triangular face or wedge with a quadrilateral face
    """

    if len(points) == 3:
        area = heronFormula(points)
    elif len(points) == 4:
        centroid =  calculateEdgeCentroid(points)
        #print (points)
        area = 0
        area1 = heronFormula([points[0],points[1],centroid])
        #print("area1", area1)
        area2 = heronFormula([points[1],points[2],centroid])
        #print("area2", area2)
        area3 = heronFormula([points[2],points[3],centroid])
        #print("area3", area3)
        area4 = heronFormula([points[3],points[0],centroid])
        #print("area4", area4)
        area = area1 + area2 + area3 + area4
    else:
        centroid =  calculateEdgeCentroid(points)
        area = 0
        for i in range(len(points)-1):
            area += heronFormula([points[i],points[i+1], centroid])

        area += heronFormula([points[-1],points[0], centroid])

    return area



def checkInt(str):
    try:
        int(str)
        return True
    except ValueError:
        return False

def read_input_file(path):
    """this function reads the user input file"""

    casePath = re.compile ("^\s*case .{1,}?(?=^\s*case|\Z)",
                           re.MULTILINE | re.DOTALL | re.IGNORECASE)
    casesVec = []
    parameters = [ 'case','time_treatment',  'activation_const',
                   'activation_dataset','activation_dataset_error',
                    'activation_file', 'activation_normalization',
                   'inlet_conc','decay_constant', 'cfd_path',
                   'molecular_diffusion','schmidt_number','div_scheme',
                    'cfd_type', 'fluent_fluid_region_name']
    try:
        fin = open(path,'r',encoding="utf8", errors='ignore')
    except IOError:
        print("couldn't open file")
        sys.exit()
    with fin:
        textBlock = fin.read()

    casesBlocks = casePath.findall(textBlock)

    for case in casesBlocks:
        parametersDic = {}
        caseLines = case.splitlines()
        for line in caseLines:
            if len(line.strip()) == 0:
                continue
            if '"' in line:
                args = line.strip().split('"')
                key = args[0].strip().lower()
                if key in parameters:
                    parametersDic[key] = args[1]
            else:
                args = line.strip().split()
                if args[0].lower() in parameters:
                    parametersDic[args[0].lower()] = args[1]

        casesVec.append(parametersDic)

    #print (casesVec)

    return casesVec

class FlunedCase:
    """
    FLUNED simulation class
    """
    def __init__(self, argDic):
        """initialize case and create FLUNED case folder"""

        self.case = argDic['case']
        if 'decay_constant' in argDic:
            self.decay_constant = float(argDic['decay_constant'])
        else:
            self.decay_constant = 0
        if 'activation_const' in argDic:
            self.activation_const = float(argDic['activation_const'])
        else:
            self.activation_const = 0
        if 'activation_normalization' in argDic:
            self.activation_normalization = float(
                argDic['activation_normalization'])
        else:
            self.activation_normalization = 0
        if 'activation_dataset' in argDic:
            self.activation_dataset = argDic['activation_dataset']
        else:
            self.activation_dataset = ''

        if 'activation_dataset_error' in argDic:
            self.activation_dataset_error = argDic['activation_dataset_error']
        else:
            self.activation_dataset_error = ''

        if 'activation_file' in argDic:
            self.activation_file = os.path.normcase(
                    argDic['activation_file'])
        else:
            self.activation_file = ''

        if 'cfd_type' not in argDic:
            print ("ERROR: cfd type not specified")
            sys.exit()
        else:
            self.cfd_type = argDic['cfd_type'].lower()

        if 'time_treatment' not in argDic:
            print ("ERROR: type of time treatment not specified")
            sys.exit()
        elif (argDic['time_treatment'].lower() not in
                ['steadystate','transient']):
            print ("ERROR: type of time treatment not recognized")
            sys.exit()
        else:
            self.time_treatment = argDic['time_treatment'].lower()

        if self.cfd_type in ['fluent-h5-multi','fluent-multi']:
            if 'fluent_fluid_region_name' not in argDic:
                print ("ERROR: name of the fluid region to extract not")
                print ("specified! use parameter FLUENT_FLUID_REGION_NAME")
                sys.exit()
            else:
                self.fluent_fluid_region_name = (
                    argDic['fluent_fluid_region_name'])

        self.molecular_diffusion = float(argDic['molecular_diffusion'])
        self.schmidt_number = float(argDic['schmidt_number'])
        self.inlet_conc = float(argDic['inlet_conc'])
        self.cfd_path = os.path.normcase(argDic['cfd_path'])
        self.cfd_simulation = SimulationOF(self.cfd_path)


        if not os.path.isdir(self.cfd_path):
            raise OSError(f"Folder not found: {self.cfd_path}")

        self.fluned_path = os.path.join(self.cfd_path,self.case)

        #self.fluned_simulation = SimulationOF(self.fluned_path)

        if not os.path.exists(self.fluned_path):
            os.mkdir(self.fluned_path)

        self.num_internal_cells = 0


        return

    def generate_zero_ta(self):
        """
        this function generate the Ta file at t=0
        """


        zero_folder = os.path.join(self.fluned_path,'0')
        zero_ta_path = os.path.join(zero_folder,'Ta')


        ta_header_text="""
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Ta;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

"""
        end_text = """

}
"""

        with open(zero_ta_path,'w',encoding='utf-8') as fw:
            fw.write(ta_header_text)

            for face in self.faces:
                fw.write("    " + face['faceID'] + '\n    {\n')
                if face['type'] == 'wall':
                    fw.write("        type            zeroGradient;\n")
                elif face['type'] == 'outlet':
                    fw.write("        type            zeroGradient;\n")
                elif face['type'] == 'inlet':
                    fw.write("        type            fixedValue;\n")
                    vals_string="        value          uniform 0;\n"
                    fw.write(vals_string)
                else:
                    print ("ERROR face type not recognized")
                    print (face['type'])
                    sys.exit()

                fw.write("    }\n" )

            fw.write(end_text)

        return

    def generate_zero_tr(self):
        """
        this function generate the Tr file at t=0
        """


        zero_folder = os.path.join(self.fluned_path,'0')
        zero_tr_path = os.path.join(zero_folder,'Tr')


        tr_header_text="""
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Tr;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

"""
        end_text = """

}
"""

        with open(zero_tr_path,'w',encoding='utf-8') as fw:
            fw.write(tr_header_text)

            for face in self.faces:
                fw.write("    " + face['faceID'] + '\n    {\n')
                if face['type'] == 'wall':
                    fw.write("        type            zeroGradient;\n")
                elif face['type'] == 'outlet':
                    fw.write("        type            zeroGradient;\n")
                elif face['type'] == 'inlet':
                    fw.write("        type            fixedValue;\n")
                    if self.time_treatment == 'steadystate':
                        vals_string="        value          uniform 0;\n"
                    elif self.time_treatment == 'transient':
                        vals_string="        value          uniform 1;\n"
                    fw.write(vals_string)
                else:
                    print ("ERROR face type not recognized")
                    print (face['type'])
                    sys.exit()

                fw.write("    }\n" )

            fw.write(end_text)

        return

    def generate_zero_td(self):
        """
        this function generate the Td file at t=0
        """


        zero_folder = os.path.join(self.fluned_path,'0')
        zero_td_path = os.path.join(zero_folder,'Td')


        td_header_text="""
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Td;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

"""
        end_text = """

}

"""

        with open(zero_td_path,'w',encoding='utf-8') as fw:
            fw.write(td_header_text)

            for face in self.faces:
                fw.write("    " + face['faceID'] + '\n    {\n')
                if face['type'] == 'wall':
                    fw.write("        type            zeroGradient;\n")
                elif face['type'] == 'outlet':
                    fw.write("        type            zeroGradient;\n")
                elif face['type'] == 'inlet':
                    fw.write("        type            fixedValue;\n")
                    vals_string="        value         uniform {};\n"
                    fw.write(vals_string.format(self.inlet_conc))
                else:
                    print ("ERROR face type not recognized")
                    print (face['type'])
                    sys.exit()

                fw.write("    }\n" )

            fw.write(end_text)

        return

    def generate_zero_t(self):
        """
        this function generate the T file at t=0
        """


        zero_folder = os.path.join(self.fluned_path,'0')
        zero_t_path = os.path.join(zero_folder,'T')


        t_header_text="""
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      T;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

"""
        end_text = """

}

"""

        with open(zero_t_path,'w',encoding='utf-8') as fw:
            fw.write(t_header_text)

            for face in self.faces:
                fw.write("    " + face['faceID'] + '\n    {\n')
                if face['type'] == 'wall':
                    fw.write("        type            zeroGradient;\n")
                elif face['type'] == 'outlet':
                    fw.write("        type            zeroGradient;\n")
                elif face['type'] == 'inlet':
                    fw.write("        type            fixedValue;\n")
                    vals_string="        value         uniform {};\n"
                    fw.write(vals_string.format(self.inlet_conc))
                else:
                    print ("ERROR face type not recognized")
                    print (face['type'])
                    sys.exit()

                fw.write("    }\n" )

            fw.write(end_text)

        return

    def generateConstantFiles(self):
        """This function creates the files in the constant folder"""
        transportPropertiesText = """
/*-------------------------------*- C++ -*----------------------------------\\
| =========                |                                                 |
| \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration    | Version:  v2012                                 |
|   \\  /    A nd          | Website:  www.openfoam.com                      |
|    \\/     M anipulation |                                                 |
\*--------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      transportProperties;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


transportModel  Newtonian;


DT             DT [0 2 -1 0 0 0 0] {};

lambda         lambda [0 0 -1 0 0 0 0] {};

Sct            Sct [ 0 0 0 0 0 0 0 ] {};

// ************************************************************************ //
"""

        constantFolder = os.path.join(self.fluned_path,'constant')
        transpPropPath = os.path.join(constantFolder,'transportProperties')


        with open(transpPropPath,'w') as fw:
            fw.write(transportPropertiesText.format(
                self.molecular_diffusion,
                self.decay_constant,
                self.schmidt_number))

        return


    def create_case_folders_fluent(self):
        """
        This function creates the folders for the FLUENT simulation converted
        into the OpenFOAM format
        """


        case_folder = self.cfd_path

        # T=0 folder
        zero_folder = os.path.join(case_folder,'0')
        if not os.path.exists(zero_folder):
            os.mkdir(zero_folder)

        # T=1 folder
        one_folder = os.path.join(case_folder,'1')
        if not os.path.exists(one_folder):
            os.mkdir(one_folder)

        # constant folder
        const_folder = os.path.join(case_folder,'constant')
        if not os.path.exists(const_folder):
            os.mkdir(const_folder)

        # polyMesh folder
        poly_folder = os.path.join(case_folder,'constant','polyMesh')
        if not os.path.exists(poly_folder):
            os.mkdir(poly_folder)

        # system folder
        sys_folder = os.path.join(case_folder,'system')
        if not os.path.exists(sys_folder):
            os.mkdir(sys_folder)

        # case foam file
        case_file = os.path.join(case_folder,'case.foam')
        pathlib.Path(case_file).touch()

        return

    def create_case_folder(self):
        """
        This function creates an OpenFoam Case - it calls some
        function for the definition of the specific files
        """

        # case Folder
        case_folder = self.fluned_path
        if not os.path.exists(case_folder):
            os.mkdir(case_folder)

        # T=0 folder
        zero_folder = os.path.join(case_folder,'0')
        if not os.path.exists(zero_folder):
            os.mkdir(zero_folder)

        # constant folder
        const_folder = os.path.join(case_folder,'constant')
        if not os.path.exists(const_folder):
            os.mkdir(const_folder)

        # polyMesh folder
        poly_folder = os.path.join(case_folder,'constant','polyMesh')
        if not os.path.exists(poly_folder):
            os.mkdir(poly_folder)

        # system folder
        sys_folder = os.path.join(case_folder,'system')
        if not os.path.exists(sys_folder):
            os.mkdir(sys_folder)

        # case foam file
        case_file = os.path.join(case_folder,'case.foam')
        pathlib.Path(case_file).touch()


        return

    def copy_last_phi(self):
        """
        this function look for the last phi file in the cfd folder
        if in the cfd simulation the flow is in m3/s it just copies it to the
        fluned folders. If the flows is in kg/s it converts it to m3/s
        """

        last_time = self.cfd_simulation.last_time
        target_file = os.path.join(self.fluned_path,'0','phi')
        origin_file = os.path.join(self.cfd_path,str(last_time),'phi')

        if self.cfd_simulation.volumetric_flag:
            shutil.copyfile( origin_file,target_file)
        else:
            print ("converting phi values to volumetric flow ...")
            self.cfd_simulation.convert_phi_to_volumetric(target_file)


        return

    def copyLastU(self):
        """ this function look for the last U file in the cfd folder"""

        folderItems = os.listdir(self.cfd_path)

        folderTimes=[int(itm) for itm in folderItems if checkInt(itm) == True]

        lastTime = max(folderTimes)

        targetFile = os.path.join(self.fluned_path,'0','U')
        originFile = os.path.join(self.cfd_path,str(lastTime),'U')
        shutil.copyfile( originFile,targetFile)


        return

    def copyLastNut(self):
        """ this function look for the last nut file in the cfd folder,
        if it does not exists it means the simulation is laminar"""

        folderItems = os.listdir(self.cfd_path)

        folderTimes=[int(itm) for itm in folderItems if checkInt(itm) == True]

        lastTime = max(folderTimes)

        targetFile = os.path.join(self.fluned_path,'0','nut')
        originFile = os.path.join(self.cfd_path,str(lastTime),'nut')
        checkFile = os.path.isfile(originFile)
        if checkFile:
            shutil.copyfile( originFile,targetFile)

        return

    def copyPolyMesh(self):
        """this function copy the poly Mesh files in the FLUNED folder"""

        sourceFolder = os.path.join(self.cfd_path,'constant','polyMesh')
        targetFolder = os.path.join(self.fluned_path,'constant','polyMesh')

        # shutil do not allow to copy into an existing folder
        dirCheck = os.path.isdir(targetFolder)
        if dirCheck:
            shutil.rmtree(targetFolder)

        shutil.copytree(sourceFolder, targetFolder)

        return

    def reconstructFaces(self):
        """ this function examines the polyMesh data and reconstruct which
         faces are input/output/wall"""

        faces = {}

        polyMeshFolder = os.path.join(self.fluned_path,'constant','polyMesh')
        boundaryFile = os.path.join(polyMeshFolder,'boundary')

        try:
            inpFile = open(boundaryFile,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open boundary file")
            sys.exit()
        with inpFile:

            blockVector = []
            text = inpFile.read()
            faceDefPat = re.compile("\d+[\n\r\s]+?\(.*?[\n\r\s]+?\)",
                                    re.MULTILINE | re.DOTALL )
            faceNumberPat = re.compile("(\d+)[\n\r\s]+?\(.*?\)",
                                       re.MULTILINE | re.DOTALL )
            boundaryPat = re.compile("[^\s]+[\n\r\s]+?\{.*?\}",
                                     re.MULTILINE | re.DOTALL )
            boundaryNamePat = re.compile("([^\s]+)[\n\r\s]+?\{.*?\}",
                                         re.MULTILINE | re.DOTALL )
            faceNPat = re.compile("nFaces.*?(\d+)")
            firstFacePat = re.compile("startFace.*?(\d+)")
            defBlock = faceDefPat.findall(text)[0]
            faceNumber = int(faceNumberPat.findall(text)[0])
            #print (defBlock)
            #print (faceNumber)
            boundaryBlocks = boundaryPat.findall(defBlock)
            #print (boundaryBlocks)

            #here it get the description of the boundary surfaces
            boundaryVec = []
            for block in boundaryBlocks:
                boundaryDic = {}
                boundaryDic['faceID'] = boundaryNamePat.findall(block)[0]
                #print (boundaryDic['faceID'])
                boundaryDic['nFaces'] = int(faceNPat.findall(block)[0])
                boundaryDic['firstFace'] = int(firstFacePat.findall(block)[0])
                boundaryDic['faces']=list(range(boundaryDic['firstFace'],
                        boundaryDic['firstFace'] + boundaryDic['nFaces']))
                boundaryVec.append(boundaryDic)
            #print (boundaryVec)


        phifilename = os.path.join(self.fluned_path,'0','phi')
        try:
            inpFile = open(phifilename,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open phi file")
            sys.exit()
        with inpFile:

            facePhiPat = re.compile("\((.{1,}?)\)",
                    re.MULTILINE | re.DOTALL )
            text = inpFile.read()
            wallFacePat = re.compile("value\s+uniform\s+0")

            #print ("face phi")
            for face in boundaryVec:
                faceBlockPat2=re.compile(face['faceID'] + "[\n\r\s]+?\{.*?\}",
                                         re.MULTILINE | re.DOTALL )
                faceBloc = faceBlockPat2.findall(text)[0]
                #print (faceBloc)
                facePhis = facePhiPat.findall(faceBloc)
                #print (facePhis)
                wallConfirm = wallFacePat.findall(faceBloc)
                #print (wallConfirm)
                if len(facePhis) == 0 and len(wallConfirm)!=0 :
                    face['type'] = "wall"
                    face['phis'] = np.zeros(face['nFaces'])
                else:
                    phiList=facePhiPat.findall(faceBloc)[0].strip().split('\n')
                    phiList = np.array([float(val) for val in phiList])


                    if all(i >= 0 for i in phiList):
                        face['type'] = 'outlet'
                    elif all(i <= 0 for i in phiList):
                        face['type'] = 'inlet'
                    else:
                        print("Error, phis with mixed sign in boundary")
                        sys.exit()



        self.faces = boundaryVec


        return


    def generate_system_file(self):
        """
        this function creates the files needed in the system folder for an
        openFOAM simulation - in later development it will apply the case
        parameters
        """

        control_dict_text_transient = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}

application    FLUNED-solver;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         100;

deltaT          0.1;

writeControl    timeStep;

writeInterval   100;

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

adjustTimeStep  yes;

maxCo           0.5;




functions
{

"""

        control_dict_text = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}

application    FLUNED-solver;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         1000;

deltaT          1;

writeControl    timeStep;

writeInterval   100;

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;





functions
{

"""

        vol_flow_text = """
    {}
    {{

        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        patch   {};
        fields  (phi);

        operation sum;
        regionType  patch;
        name        $patch;

        writeFields     false;
        writeControl {};

    }}

 """
        vol_t_flow_text = """
    {}
    {{

        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        patch   {};
        fields  (T);
        weightField phi;

        operation sum;
        regionType  patch;
        name        $patch;

        writeFields     false;
        writeControl {};

    }}

 """
        vol_td_flow_text = """
    {}
    {{

        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        patch   {};
        fields  (Td);
        weightField phi;

        operation sum;
        regionType  patch;
        name        $patch;

        writeFields     false;
        writeControl {};

    }}
 """

        vol_ta_flow_text = """
    {}
    {{

        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        patch   {};
        fields  (Ta);
        weightField phi;

        operation sum;
        regionType  patch;
        name        $patch;

        writeFields     false;
        writeControl {};

    }}

 """

        vol_tr_flow_text = """
    {}
    {{

        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        patch   {};
        fields  (Tr);
        weightField phi;

        operation sum;
        regionType  patch;
        name        $patch;

        writeFields     false;
        writeControl {};

    }}

 """
        system_folder = os.path.join(self.fluned_path,'system')
        control_dict_path = os.path.join(system_folder,'controlDict')



        with open(control_dict_path,'w',encoding='utf-8') as fw:
            if self.time_treatment =='steadystate':
                fw.write(control_dict_text)
                for face in self.faces:
                    if face['type'] in ['inlet','outlet']:
                        fw.write(vol_flow_text.format(
                                     "volFlow-"+face['faceID'],
                                     face['faceID'],
                                     'outputTime'))
                        fw.write(vol_t_flow_text.format(
                                     "volTFlow-"+face['faceID'],
                                     face['faceID'],
                                     'outputTime'))
                        fw.write(vol_tr_flow_text.format(
                                     "volTrFlow-"+face['faceID'],
                                     face['faceID'],
                                     'outputTime'))
                        fw.write(vol_td_flow_text.format(
                                     "volTdFlow-"+face['faceID'],
                                     face['faceID'],
                                     'outputTime'))
                        fw.write(vol_ta_flow_text.format(
                                     "volTaFlow-"+face['faceID'],
                                     face['faceID'],
                                     'outputTime'))
            elif self.time_treatment =='transient':
                fw.write(control_dict_text_transient)
                for face in self.faces:
                    if face['type'] in ['inlet','outlet']:
                        fw.write(vol_flow_text.format(
                                     "volFlow-"+face['faceID'],
                                     face['faceID'],
                                     'timeStep'))
                        fw.write(vol_t_flow_text.format(
                                     "volTFlow-"+face['faceID'],
                                     face['faceID'],
                                     'timeStep'))
                        fw.write(vol_tr_flow_text.format(
                                     "volTrFlow-"+face['faceID'],
                                     face['faceID'],
                                     'timeStep'))
                        fw.write(vol_td_flow_text.format(
                                     "volTdFlow-"+face['faceID'],
                                     face['faceID'],
                                     'timeStep'))
                        fw.write(vol_ta_flow_text.format(
                                     "volTaFlow-"+face['faceID'],
                                     face['faceID'],
                                     'timeStep'))


            fw.write("}")


        fv_scheme_text = """
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}}

ddtSchemes
{{
    default         {};
}}

gradSchemes
{{
    default         cellLimited Gauss linear 0.5;
}}

divSchemes
{{
    default         none;
    div(phi,T)      Gauss linearUpwind default;
    div(phi,Ta)      Gauss linearUpwind default;
    div(phi,Td)      Gauss linearUpwind default;
    div(phi,Tr)      Gauss linearUpwind default;
}}

laplacianSchemes
{{
    default         none;
    laplacian(DT,T) Gauss linear corrected;
    laplacian(Dturbulent,T) Gauss linear corrected;
    laplacian(DT,Tr) Gauss linear corrected;
    laplacian(Dturbulent,Tr) Gauss linear corrected;
    laplacian(DT,Ta) Gauss linear corrected;
    laplacian(Dturbulent,Ta) Gauss linear corrected;
    laplacian(DT,Td) Gauss linear corrected;
    laplacian(Dturbulent,Td) Gauss linear corrected;
}}

interpolationSchemes
{{
    default         linear;
}}

snGradSchemes
{{
    default         limited 1;
}}

"""

        schemes_path = os.path.join(system_folder,'fvSchemes')


        with open(schemes_path,'w',encoding='utf-8') as fw:
            if self.time_treatment == 'steadystate':
                fw.write(fv_scheme_text.format('steadyState'))
            if self.time_treatment == 'transient':
                fw.write(fv_scheme_text.format('Euler'))

        fv_solution_text = """
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


solvers
{{
    T
    {{
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-07;
        relTol          {};
    }}
    Tr
    {{
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-07;
        relTol          {};
    }}
    Ta
    {{
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-07;
        relTol          {};
    }}
    Td
    {{
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-07;
        relTol          {};
    }}
}}

SIMPLE
{{
    nNonOrthogonalCorrectors 0;
    consistent      yes;

    residualControl
    {{

        T               1e-6;
        Tr              1e-6;
        Ta              1e-6;
        Td              1e-6;

    }}


}}

"""
        solution_path = os.path.join(system_folder,'fvSolution')


        with open(solution_path,'w',encoding='utf-8') as fw:
            if self.time_treatment == 'steadystate':
                fw.write(fv_solution_text.format(0.01,0.01,0.01,0.01))
            elif self.time_treatment == 'transient':
                fw.write(fv_solution_text.format(0,0,0,0))

        parallel_dict_text = """
FoamFile
{
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains 4;

method          scotch;

// *********************************************************************** //
"""
        parallel_dict_path = os.path.join(system_folder,'decomposeParDict')
        with open(parallel_dict_path,'w',encoding='utf-8') as fw:
                fw.write(parallel_dict_text)

        return


    def launchVolFuncObjects(self):
        """this function launch the utilities that calculates the volumes
        - this field is needed only for the activation rate
        interpolation"""

        print ("calculating volumes  ...")
        origFolder = os.getcwd()
        os.chdir(self.fluned_path)
        launchVolumes = "postProcess -func writeCellVolumes".split()
        with open('log', "a") as outfile:
            proc = subprocess.Popen(launchVolumes, stdout=outfile).wait()
        os.chdir(origFolder)

        return

    def launchCentroidFuncObjects(self):
        """this function launch the utilities that calculates the volume and
        the centroids - this fields are needed only for the activation rate
        interpolation"""

        print ("calculating centroids ...")
        origFolder = os.getcwd()
        os.chdir(self.fluned_path)
        launchCentroids = "postProcess -func writeCellCentres".split()
        with open('log', "a") as outfile:
            proc = subprocess.Popen(launchCentroids, stdout=outfile).wait()
        os.chdir(origFolder)

        return

    def launch_solver(self):
        """
        this function launches the scalar calculation
        """

        print ("launching FLUNED solver ...")
        current_folder = os.getcwd()
        os.chdir(self.fluned_path)
        launch_calc_string = "FLUNED-solver".split()
        with open('simulation_log', "a", encoding='utf-8') as outfile:
            subprocess.Popen(launch_calc_string, stdout=outfile).wait()
        os.chdir(current_folder)

        return


    def readVolumes(self):
        """ this function reads the volumes from the V file located in the
        zero folder """

        # common patterns
        internalBlockPat = re.compile("internalField.*?\((.{1,}?)\)",
                                          re.MULTILINE | re.DOTALL )

        vFile = os.path.join(self.fluned_path,'0','V')
        try:
            inpFile = open(vFile,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open Volume V file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            numInternalBlocks = internalBlockPat.findall(text)
            internalVolumes = numInternalBlocks[0].split('\n')[1:-1]
            self.Volumes = np.array([float(val) for val in internalVolumes])

        return

    def readCentroids(self):
        """ this function reads the centroids from the 0 folder"""

        # common patterns
        internalBlockPat = re.compile("internalField.*?\((.{1,}?)\n\\s*\)",
                                          re.MULTILINE | re.DOTALL )

        cFile = os.path.join(self.fluned_path,'0','C')
        try:
            inpFile = open(cFile,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open  C file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            numInternalBlocks = internalBlockPat.findall(text)
            internalCentroids = numInternalBlocks[0].split('\n')[1:]
            internalCentroids=[val.strip('()') for val in internalCentroids]
            self.Centroids = np.array([[float(val) for val in v.split()]
                                for v in internalCentroids])


        return


    def write_cell_zones(self):
        """
        function to write the cell zones file in the polyMesh folder
        """
        cell_zones_file_path = os.path.join(self.cfd_path,'constant',
                                           'polyMesh', 'cellZones')

        cell_zones_string = """
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       regIOobject;
    location    "constant/polyMesh";
    object      cellZones;

}

0()

"""
        with open(cell_zones_file_path,'w',encoding='utf-8') as fo2:
            fo2.write(cell_zones_string)

        return



    def writeNeighbour_multi_h5(self):
        """
        this function writes all the cell neighbours,
        it selects those that belongs to the fluid partition
        """

        print ("writing multiblock openFOAM neighbour ... ")

        filename = os.path.join(self.cfd_path,self.casH5file)

        neighPat = re.compile(".*/faces/c1/\d+\Z", re.IGNORECASE)

        neighbourFilePath = os.path.join(self.cfd_path,'constant','polyMesh',
                                      'neighbour')

        neighHeader = """
/*--------------------------------*- C++ -*-------------------------------*\\
| =========                 |                                               |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\    /   O peration     | Version:  2112                                |
|   \\  /    A nd           | Website:  www.openfoam.com                    |
|    \\/     M anipulation  |                                               |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       labelList;
    location    "constant/polyMesh";
    object      neighbour;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""
        neighData = get_h5_dataset_multi(filename,neighPat)


        for face in self.faceList:
            if (face['fluid']) and face['faceType'] == 'internal':
                nWaterNeigh = face['nFaces']
                neighTable = extract_multiblock(face['minID'],face['maxID'],
                            neighData)
                break

        with  open(neighbourFilePath,'w') as fo:

            fo.write(neighHeader)
            fo.write(str(nWaterNeigh) + '\n')
            fo.write("(\n")
            neStr = "{:d}\n"
            for val in neighTable:
                fo.write(neStr.format(int(val -self.fluid_cellID_min)))
            fo.write(")")


        return

    def writeFaces_multi_h5(self):
        """
        this function write the face definition of the mesh extracted from a
        multiblock fluent simulation
        """

        print ("writing multiblock openFOAM faces ... ")

        filename = os.path.join(self.cfd_path,self.casH5file)

        facePat = re.compile(".*/faces/nodes/\d+/nnodes\Z", re.IGNORECASE)
        face2Pat = re.compile(".*/faces/nodes/\d+/nodes\Z", re.IGNORECASE)

        facesFilePath = os.path.join(self.cfd_path,'constant','polyMesh',
                                      'faces')

        facesHeader = """
/*-------------------------------*- C++ -*--------------------------------*\\
| =========                |                                                |
| \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox          |
|  \\    /   O peration    | Version:  2112                                 |
|   \\  /    A nd          | Website:  www.openfoam.com                     |
|    \\/     M anipulation |                                                |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       faceList;
    location    "constant/polyMesh";
    object      faces;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

        self.faces1H5Path,nFacesVec = get_h5_path_dataset(filename,facePat)
        self.faces2H5Path,facesDef = get_h5_path_dataset(filename,face2Pat)




        pointIDTemp = 1
        uniquePoints = np.zeros((0,), dtype =np.uint64)
        for face in self.faceList:
            # store the number of nodes per face
            nPointsTemp = nFacesVec[(face['minID']-1):face['maxID']]
            # number of nodes needed to define one face
            face['nPoints'] = sum(nPointsTemp)
            face['minPointID'] = pointIDTemp
            face['maxPointID'] = face['minPointID'] + face['nPoints'] - 1
            pointIDTemp += face['nPoints']
            #print (face['faceName'], face['minPointID'],face['maxPointID'])
            if face['fluid'] == True:
                pointsTemp = facesDef[(face['minPointID']-1):
                        face['maxPointID']]
                #print ("fluid face", min(pointsTemp),max(pointsTemp))
                uniquePoints = np.append(uniquePoints,pointsTemp)


        self.uniquePoints = np.unique(uniquePoints)

        #print ("CHECK")
        #print (len(facesDef))
        #print (pointIDTemp)

        faceList = sorted(self.faceList, reverse=False,
                                      key=lambda x:x['newOrderID'])


        with  open(facesFilePath,'w') as fo:
            fo.write(facesHeader)
            fo.write(str(self.nWaterFaces) + '\n')
            fo.write("(\n")
            j=0
            for face in faceList:
                if face['fluid'] == False:
                    continue
                nPointsTemp=nFacesVec[(face['minID']-1):face['maxID']]
                pointsTemp = facesDef[(face['minPointID']-1):
                        face['maxPointID']]
                i = 0
                #print (len(nPointsTemp))
                #print (sum(nPointsTemp))
                #print (len(pointsTemp))
                for val in nPointsTemp:
                    faceStr1 = "{:d}("
                    newString = faceStr1.format(val)
                    pointList = reversed(pointsTemp[i:(i+val)])
                    for pt in pointList:
                        ptNew = np.searchsorted(self.uniquePoints,pt)
                        newString += " {:d} ".format(int(ptNew))
                    newString += " )\n"
                    fo.write(newString)
                    i += val
                    j += val

            fo.write(")")





        faceZonesFilePath = os.path.join(self.cfd_path,'constant',
                                           'polyMesh', 'faceZones')

        faceZones = """
/*-------------------------------*- C++ -*--------------------------------*\\
| =========                 |                                               |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\    /   O peration     | Version:  2112                                |
|   \\  /    A nd           | Website:  www.openfoam.com                    |
|    \\/     M anipulation  |                                               |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       regIOobject;
    location    "constant/polyMesh";
    object      faceZones;

}
// * * * * * * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * //

0()

// ********************************************************************** //
"""
        with open(faceZonesFilePath,'w') as fo2:
            fo2.write(faceZones)



        return

    def writeBoundary_multi_h5(self):
        """
        this function writes the boundary file for a simulation
        extracted from a multi block fluent simulation
        the self.faceList should be already arranged for the work
        """

        boundaryFilePath = os.path.join(self.cfd_path,'constant','polyMesh',
                                      'boundary')
        boundaryHeader = """
/*--------------------------------*- C++ -*------------------------------*\\
| =========                 |                                               |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\    /   O peration     | Version:  2112                                |
|   \\  /    A nd           | Website:  www.openfoam.com                    |
|    \\/     M anipulation  |                                               |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       polyBoundaryMesh;
    location    "constant/polyMesh";
    object      boundary;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""
        nWaterFaces = 0
        for face in self.faceList:
            if (face['fluid'] == True) and face['faceType'] != 'internal':
                nWaterFaces += 1

        faceList = sorted(self.faceList, reverse=False,
                                      key=lambda x:x['newOrderID'])


        with open(boundaryFilePath,'w') as fi:
            fi.write(boundaryHeader)
            fi.write("{:d}\n".format(nWaterFaces ))
            fi.write("(\n")
            for face in faceList:
                if face['fluid'] == False:
                    continue
                if face['faceType'] == 'internal':
                    continue

                fi.write("   {}\n".format(face['faceName']))
                fi.write("    {\n")

                if face["faceType"] == 'wall':
                    fi.write("        type            wall;\n")
                else:
                    fi.write("        type            patch;\n")

                fi.write("        nFaces          {:d};\n".format(
                          face['nFaces']))
                fi.write("        startFace       {:d};\n".format(
                         int(face['newMinID'])-1))
                fi.write("    }\n")
            fi.write(")\n")



        return


        return

    def writeOwner_multi_h5(self):
        """ This function writes all the cell owners,
        it selects those that belong to the fluid partition"""

        print ("writing multiblock openFOAM owner ... ")

        filename = os.path.join(self.cfd_path,self.casH5file)

        ownerPat = re.compile(".*/faces/c0/\d+\Z", re.IGNORECASE)

        ownerFilePath = os.path.join(self.cfd_path,'constant','polyMesh',
                                      'owner')


        ownerHeader = """
/*-------------------------------*- C++ -*--------------------------------*\\
| =========                 |                                               |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\    /   O peration     | Version:  2112                                |
|   \\  /    A nd           | Website:  www.openfoam.com                    |
|    \\/     M anipulation  |                                               |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       labelList;
    location    "constant/polyMesh";
    object      owner;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

        self.ownerH5Path,ownerList = get_h5_path_dataset(filename,ownerPat)


        ownerTable = np.zeros((self.nWaterFaces,), dtype=np.uint64)



        for face in self.faceList:
            if face['fluid'] == True:

                ownerTable[(face['newMinID']-1):face['newMaxID']] = (
                    ownerList[(face['minID']-1):face['maxID']])


        #check = np.count_nonzero(ownerTable == 0)
        #print (check)

        with  open(ownerFilePath,'w') as fo:
            fo.write(ownerHeader)
            fo.write(str(self.nWaterFaces) + '\n')
            fo.write("(\n")
            owStr = "{:d}\n"
            for val in ownerTable:
                fo.write(owStr.format(int(val - self.fluid_cellID_min)))
            fo.write(")")


        return

    def writeNodes_multi(self):
        """
        this function writes the node of a converted multiblock fluent
        simulation. This function limits the writing to the points present
        in the fluid region written in the self.uniquePoints attribute
        """

        print ("writing multiblock openFOAM points ... ")

        filename = os.path.join(self.cfd_path,self.casfile)

        pointsPat = re.compile("^\(3010\s+\(.*?\).*\)",
                                          re.MULTILINE | re.DOTALL )

        pointsFilePath = os.path.join(self.cfd_path,'constant','polyMesh',
                                      'points')


        pointsHeader = """
/*-------------------------------*- C++ -*--------------------------------*\\
| =========                 |                                               |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\    /   O peration     | Version:  2112                                |
|   \\  /    A nd           | Website:  www.openfoam.com                    |
|    \\/     M anipulation  |                                               |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       vectorField;
    location    "constant/polyMesh";
    object      points;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""
        #self.nodesH5Path,points = get_h5_path_dataset(filename,nodesPat)

        get_fluent_binarray_double(filename,'3010',3)



        #with open(pointsFilePath,'w') as fo:

        #    fo.write(pointsHeader)
        #    fo.write(str(len(self.uniquePoints)) + '\n')
        #    fo.write("(\n")
        #    ptStr = "({:24.18e} {:24.18e} {:24.18e})\n"
        #    for index in self.uniquePoints:
        #
        #        cord = points[int(index-1)]

        #        fo.write(ptStr.format(*cord))

        #    fo.write(")")
        #




        pointsZonesFilePath = os.path.join(self.cfd_path,'constant',
                                           'polyMesh', 'pointZones')

        pointZones = """
/*--------------------------------*- C++ -*-------------------------------*\\
| =========                 |                                               |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\    /   O peration     | Version:  2112                                |
|   \\  /    A nd           | Website:  www.openfoam.com                    |
|    \\/     M anipulation  |                                               |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       regIOobject;
    location    "constant/polyMesh";
    object      pointZones;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

0()

// *********************************************************************** //
"""
        with open(pointsZonesFilePath,'w') as fo2:
            fo2.write(pointZones)

        return

    def writeNodes_multi_h5(self):
        """
        this function writes the node of a converted multiblock fluent
        simulation. This function limits the writing to the points present
        in the fluid region written in the self.uniquePoints attribute
        """

        print ("writing multiblock openFOAM points ... ")

        filename = os.path.join(self.cfd_path,self.casH5file)

        nodesPat = re.compile(".*/nodes/coords/\d+\Z", re.IGNORECASE)

        pointsFilePath = os.path.join(self.cfd_path,'constant','polyMesh',
                                      'points')


        pointsHeader = """
/*-------------------------------*- C++ -*--------------------------------*\\
| =========                 |                                               |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\    /   O peration     | Version:  2112                                |
|   \\  /    A nd           | Website:  www.openfoam.com                    |
|    \\/     M anipulation  |                                               |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       vectorField;
    location    "constant/polyMesh";
    object      points;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""
        self.nodesH5Path,points = get_h5_path_dataset(filename,nodesPat)


        with open(pointsFilePath,'w') as fo:

            fo.write(pointsHeader)
            fo.write(str(len(self.uniquePoints)) + '\n')
            fo.write("(\n")
            ptStr = "({:24.18e} {:24.18e} {:24.18e})\n"
            for index in self.uniquePoints:

                cord = points[int(index-1)]

                fo.write(ptStr.format(*cord))

            fo.write(")")





        pointsZonesFilePath = os.path.join(self.cfd_path,'constant',
                                           'polyMesh', 'pointZones')

        pointZones = """
/*--------------------------------*- C++ -*-------------------------------*\\
| =========                 |                                               |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\    /   O peration     | Version:  2112                                |
|   \\  /    A nd           | Website:  www.openfoam.com                    |
|    \\/     M anipulation  |                                               |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       regIOobject;
    location    "constant/polyMesh";
    object      pointZones;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

0()

// *********************************************************************** //
"""
        with open(pointsZonesFilePath,'w') as fo2:
            fo2.write(pointZones)

        return


    def getH5files(self):


        casH5files = []
        datH5files = []
        casH5FilePat = re.compile("\.cas.h5\Z", re.IGNORECASE)
        datH5FilePat = re.compile("\.dat.h5\Z", re.IGNORECASE)

        folder =self.cfd_path

        for filename in os.listdir(folder):
            casH5file = casH5FilePat.findall(filename)
            datH5file = datH5FilePat.findall(filename)
            if len(casH5file) == 1:
                casH5files.append(filename)
            if len(datH5file) == 1:
                datH5files.append(filename)

        if len(casH5files) == 1:
            self.casH5file = casH5files[0]
        else:
            print ("ERROR zero or more than one cas.h5 files")
            sys.exit()

        if len(datH5files) == 1:
            self.datH5file = datH5files[0]
        else:
            print ("ERROR zero or more than one dat.h5 files")
            sys.exit()

        return

    def getCASDATfiles(self):


        casfiles = []
        datfiles = []
        casFilePat = re.compile("\.cas\Z", re.IGNORECASE)
        datFilePat = re.compile("\.dat\Z", re.IGNORECASE)

        folder =self.cfd_path

        for filename in os.listdir(folder):
            casfile = casFilePat.findall(filename)
            datfile = datFilePat.findall(filename)
            if len(casfile) == 1:
                casfiles.append(filename)
            if len(datfile) == 1:
                datfiles.append(filename)

        if len(casfiles) == 1:
            self.casfile = casfiles[0]
        else:
            print ("ERROR zero or more than one cas file")
            sys.exit()

        if len(datfiles) == 1:
            self.datfile = datfiles[0]
        else:
            print ("ERROR zero or more than one dat file")
            sys.exit()

        return

    def generateTrSourceFile(self):
        """
        this file create the source file for the time residency fictiotious
        scalar in the zero folder.
        There are two modes:
        1) steady state: the source is 1 in every cell
        2) transient the source is 0 in every cell
         """



        zeroFolder = os.path.join(self.fluned_path,'0')
        zeroSourcePath = os.path.join(zeroFolder,'TrSource')


        sHeaderText="""
/*------------------------------*- C++ -*----------------------------------*\
| =========               |                                                 |
| \\      /  F ield       | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration   | Version:  2.3.1                                 |
|   \\  /    A nd         | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulatio |                                                 |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      TrSource;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

"""
        intFieldText = """

internalField   uniform {};

        """
        boundaryText = """
boundaryField
{

"""
        closerText = """

}


// *********************************************************************** //
"""

        with open(zeroSourcePath,'w') as fw:
            fw.write(sHeaderText)

            if self.time_treatment == 'steadystate':
                fw.write(intFieldText.format(1))
            elif self.time_treatment == 'transient':
                fw.write(intFieldText.format(0))

            fw.write(boundaryText)

            for face in self.faces:
                fw.write("    " + face['faceID'] + '\n    {\n')
                fw.write("        type            fixedValue;\n")
                valString="        value          uniform 0;\n"
                fw.write(valString)
                fw.write("    }\n" )

            fw.write(closerText)



        return


    def generateSourceFile(self):
        """this file create the source file in the zero folder.
        There are three modes:
        1) no activation, then the file contains only zeros.
        2) constant activation (with input value)
        3) with a source file
         """


        if (self.activation_file == ''):
            if self.activation_const == 0:
                # case with no rad source
                activ_sources = [0 for i in range(self.num_internal_cells)]
            else:
                # case with constant value source
                self.readVolumes()
                activ_sources = ([self.activation_const*vol for vol
                    in self.Volumes])

        else:
            #case with source file
            self.launchCentroidFuncObjects()
            self.readVolumes()
            self.readCentroids()

            # 1.sample the activation file
            print ("Sampling Reaction Rate file ... ")
            sampledRates = sample_coordinates_VTK(self.activation_file,
                                                self.activation_dataset,
                                                self.Centroids)

            # 1.1 if present sample the vtk file to get the error array
            if self.activation_dataset_error != '':
                print ("Sampling Reaction Rate MCNP errors  ... ")
                sampledStatErr = sample_coordinates_VTK(self.activation_file,
                                            self.activation_dataset_error,
                                            self.Centroids)

            # 2.use the activation const as a factor
            if self.activation_const == 0:
                factor = 1
            else:
                factor = self.activation_const
            activ_sources = [factor*rate for rate in sampledRates ]

            # 3.apply a normalization factor if provided
            if self.activation_normalization != 0:
                vec = ([rate*vol for rate,vol in
                    zip(activ_sources,self.Volumes)])

                totalSampled = sum(vec)
                print ("total sampled atoms/s")
                print (totalSampled)
                normFactor = self.activation_normalization/totalSampled

                activ_sources = [rate*normFactor for rate in activ_sources ]


                nVec = ([rate*vol for rate,vol in
                    zip(activ_sources,self.Volumes)])

                print ("new total sampled atoms/s")
                print (sum(nVec))

        self.activation_sources  = activ_sources

        zeroFolder = os.path.join(self.fluned_path,'0')
        zeroSourcePath = os.path.join(zeroFolder,'Source')


        sHeaderText="""
/*------------------------------*- C++ -*----------------------------------*\
| =========               |                                                 |
| \\      /  F ield       | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration   | Version:  2.3.1                                 |
|   \\  /    A nd         | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulatio |                                                 |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Source;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

internalField   nonuniform List<scalar>
"""
        boundaryText = """
boundaryField
{

"""
        closerText = """

}


// *********************************************************************** //
"""

        with open(zeroSourcePath,'w') as fw:
            fw.write(sHeaderText)
            fw.write("{:d}\n".format(self.num_internal_cells))
            fw.write("(\n")
            for val in self.activation_sources:
                fw.write("{:e}\n".format(val))
            fw.write(")\n;\n\n")

            fw.write(boundaryText)

            for face in self.faces:
                fw.write("    " + face['faceID'] + '\n    {\n')
                fw.write("        type            fixedValue;\n")
                valString="        value          uniform 0;\n"
                fw.write(valString)
                fw.write("    }\n" )

            fw.write(closerText)


        if self.activation_dataset_error != '':
            zeroSourceErrorPath = os.path.join(zeroFolder,'SourceError')


            eHeaderText="""
/*------------------------------*- C++ -*----------------------------------*\
| =========               |                                                 |
| \\      /  F ield       | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration   | Version:  2.3.1                                 |
|   \\  /    A nd         | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulatio |                                                 |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      SourceError;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

internalField   nonuniform List<scalar>
"""
            with open(zeroSourceErrorPath,'w') as fw:
                fw.write(eHeaderText)
                fw.write("{:d}\n".format(self.num_internal_cells))
                fw.write("(\n")
                for val in sampledStatErr:
                    fw.write("{:e}\n".format(val))
                fw.write(")\n;\n\n")

                fw.write(boundaryText)

                for face in self.faces:
                    fw.write("    " + face['faceID'] + '\n    {\n')
                    fw.write("        type            fixedValue;\n")
                    valString="        value          uniform 0;\n"
                    fw.write(valString)
                    fw.write("    }\n" )

                fw.write(closerText)

        return


    def writeSpeed_multi_h5(self):
        """
        this function write the U velocity files for t=1 and t=0
        """

        print ("writing multiblock openFOAM U files ... ")

        filename = os.path.join(self.cfd_path,self.datH5file)

        uCellXPat = re.compile(".*/cells/SV_U/.*",
                                     re.IGNORECASE)
        uCellYPat = re.compile(".*/cells/SV_V/.*",
                                     re.IGNORECASE)
        uCellZPat = re.compile(".*/cells/SV_W/.*",
                                     re.IGNORECASE)

        uFaceXPat = re.compile(".*/faces/SV_U/.*",
                                     re.IGNORECASE)
        uFaceYPat = re.compile(".*/faces/SV_V/.*",
                                     re.IGNORECASE)

        uFaceZPat = re.compile(".*/faces/SV_W/.*",
                                     re.IGNORECASE)


        uOneFilePath = os.path.join(self.cfd_path,'1','U')


        xCellVelocity = get_h5_dataset_multi(filename,uCellXPat)
        yCellVelocity = get_h5_dataset_multi(filename,uCellYPat)
        zCellVelocity = get_h5_dataset_multi(filename,uCellZPat)

        xFaceVelocity = get_h5_dataset_multi(filename,uFaceXPat)
        yFaceVelocity = get_h5_dataset_multi(filename,uFaceYPat)
        zFaceVelocity = get_h5_dataset_multi(filename,uFaceZPat)


        u1Header = """
/*--------------------------------*- C++ -*-------------------------------*\\
| =========                 |                                               |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\    /   O peration     | Version:  2112                                |
|   \\  /    A nd           | Website:  www.openfoam.com                    |
|    \\/     M anipulation  |                                               |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       volVectorField;
    location    "1";
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];
"""
        faceList = sorted(self.faceList, reverse=False,
                                      key=lambda x:x['newOrderID'])

        with  open(uOneFilePath,'w') as fo:


            uStr = "({:24.18e} {:24.18e} {:24.18e})\n"
            fo.write(u1Header)
            fo.write("internalField   nonuniform List<vector>\n")
            fo.write("{:d}\n".format(self.fluid_cellN))
            fo.write("(\n")
            uXValues = extract_multiblock(self.fluid_cellID_min,
                                          self.fluid_cellID_max,
                                          xCellVelocity)
            uYValues = extract_multiblock(self.fluid_cellID_min,
                                          self.fluid_cellID_max,
                                          yCellVelocity)
            uZValues = extract_multiblock(self.fluid_cellID_min,
                                          self.fluid_cellID_max,
                                          zCellVelocity)
            for uX,uY,uZ in zip(uXValues,
                                uYValues,
                                uZValues):
                fo.write(uStr.format(uX,uY,uZ))
            fo.write(")\n")
            fo.write(";\n")
            fo.write("\n")
            fo.write("boundaryField\n")
            fo.write("{\n")

            for face in faceList:

                if face['faceType'] == 'internal':
                    continue
                if face['fluid'] :

                    fo.write("    {}\n".format(face["faceName"]))
                    fo.write("    {\n")
                    if face['faceType'] == 'wall':
                        fo.write("        type        movingWallVelocity;\n")
                        fo.write("        value       uniform (0 0 0);\n")
                    elif face['faceType'] == 'outlet':
                        fo.write("        type            zeroGradient;\n")
                    elif face['faceType'] == 'inlet':


                        uXValues = extract_multiblock(face['minID'],
                                                      face['maxID'],
                                                      xFaceVelocity)
                        uYValues = extract_multiblock(face['minID'],
                                                      face['maxID'],
                                                      yFaceVelocity)
                        uZValues = extract_multiblock(face['minID'],
                                                      face['maxID'],
                                                      zFaceVelocity)
                        fo.write("        type            fixedValue;\n")
                        fo.write("        value   nonuniform List<vector>\n")
                        fo.write("{:d}\n".format(len(uXValues)))
                        fo.write("(\n")
                        for v1,v2,v3 in zip(uXValues,uYValues,uZValues):
                            fo.write(uStr.format(v1,v2,v3))
                        fo.write(")\n")
                        fo.write(";\n")
                    fo.write("    }\n")


            fo.write("}\n")

        return


    def defineWalls_multi_h5(self):
        """
        this function just check the phi values to distinguish between
        wall, inlet and outlet
        """

        filename = os.path.join(self.cfd_path,self.datH5file)

        phiPat = re.compile(".*/faces/SV_FLUX/.*",
                                     re.IGNORECASE)

        phiValues = get_h5_dataset_multi(filename,phiPat)

        for face in self.faceList:
            if not face['fluid']  :
                continue
            if (face['fluid'])  and face['faceType'] == 'internal':
                continue


            phiTemp=extract_multiblock(face['minID'],face['maxID'],phiValues)

            if all(val == 0 for val in phiTemp):
                face["faceType"] = "wall"
            elif all(val >= 0 for val in phiTemp) :
                face["faceType"] = "outlet"
                if not all(val > 0 for val in phiTemp) :
                    print("WARNING in outlet ", face['faceName'])
                    print("some zero phis are present")
            elif all(val <= 0 for val in phiTemp) :
                face["faceType"] = "inlet"
                if not all(val < 0 for val in phiTemp) :
                    print("WARNING in inlet ", face['faceName'])
                    print("some zero phis are present")
            else:
                print("ERROR in face ", face['faceName'])
                print("nonzero phis of different signs are present")
                sys.exit()

        #for dic in self.faceList:
        #    print (dic)
        #df = pd.DataFrame(self.faceList)
        #df.to_csv('faces.csv', index=True)

        return

    def writePhi_multi_h5(self):
        """
        this function works for a multi-block fluent simulation and does
        three things:
        - write the phi values at the end of the simulation
        - scan the phi values to understand which are the wall and which the
          inlet/outlet. It updates the fluentFacesVector vector with the
          results
        """


        print ("writing multiblock openFOAM phi files ... ")

        filename = os.path.join(self.cfd_path,self.datH5file)

        phiPat = re.compile(".*/faces/SV_FLUX/.*",
                                     re.IGNORECASE)


        phiOneFilePath = os.path.join(self.cfd_path,'1','phi')

        phiValues = get_h5_dataset_multi(filename,phiPat)





        phi1Header = """
/*--------------------------------*- C++ -*-------------------------------*\\
| =========                 |                                               |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\    /   O peration     | Version:  2112                                |
|   \\  /    A nd           | Website:  www.openfoam.com                    |
|    \\/     M anipulation  |                                               |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       surfaceScalarField;
    location    "1";
    object      phi;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 3 -1 0 0 0 0];


        """


        faceList = sorted(self.faceList, reverse=False,
                                      key=lambda x:x['newOrderID'])

        with  open(phiOneFilePath,'w') as fo:

            fo.write(phi1Header)
            fo.write("internalField   nonuniform List<scalar>\n")
            # write the internal mesh
            for face in faceList:
                if (face['fluid'])  and face['faceType'] == 'internal':
                    phiTemp=extract_multiblock(face['minID'],face['maxID'],
                            phiValues)
                    fo.write("{:d}\n".format(len(phiTemp)))
                    fo.write("(\n")
                    phiStr = "{:24.18e}\n"
                    for val in phiTemp:
                        fo.write(phiStr.format(val/self.fluentDensity))
            fo.write(")\n")
            fo.write(";\n")
            fo.write("\n")
            fo.write("boundaryField\n")
            fo.write("{\n")

            for face in faceList:

                if face['faceType'] == 'internal':
                    continue
                if face['fluid'] :
                    phiTemp=extract_multiblock(face['minID'],face['maxID'],
                            phiValues)
                    fo.write("    {}\n".format(face["faceName"]))
                    fo.write("    {\n")
                    fo.write("        type            calculated;\n")
                    if face["faceType"] == "wall":
                        fo.write("        value           uniform 0;\n")
                    if face["faceType"] != "wall":
                        fo.write("        value   nonuniform List<scalar>\n")
                        phiTemp=extract_multiblock(face['minID'],
                                                   face['maxID'],
                                                   phiValues)
                        fo.write("{:d}\n".format(len(phiTemp)))
                        fo.write("(\n")
                        for val in phiTemp:
                            fo.write(phiStr.format(val/self.fluentDensity))
                        fo.write(")\n")
                        fo.write(";\n")
                    fo.write("    }\n")


            fo.write("}\n")


        return



    def getFluidFaces(self):
        """
        this function reads the zone topology information relative to
        the faces of the multi block model. Using this info it redefine the
        ranges of the faces to have the internal mesh at the beginning and
        the remaining faces later. Ideally these ranges will help the rest
        of the face definition.
        """

        filename = os.path.join(self.cfd_path,self.casfile)

        regionPat1 = re.compile('\([^()]*?\)')

        nameMatches = regionPat1.findall(self.cfdPostMeshInfoString)

        faceNames = [name.strip('()') for name in nameMatches[1:]]


        faceList = []

        for i,val in enumerate(faceNames):
            newDict = {}
            name = val.split()[0]
            zoneTypeString = val.split()[1]
            ownerName = val.split()[2]
            if len(val.split()) > 3:
                neighbourName = val.split()[3]
            else:
                neighbourName = ''

            newDict['faceName'] = name
            newDict['origOrderID'] = i

            # look for the face ID
            filtVect = [v for v in self.headerDictCas if v['index']==39]
            for entry in filtVect:
                if entry['header'].split()[2] == name:
                    newDict['regionID'] = int(entry['header'].split()[0])

            # look for the min and max ID
            filtVect = [v for v in self.headerDictCas if v['index']==3013]
            for entry in filtVect:
                headerVec = []
                for  val in entry['header'].split():
                    intVal = int(val,16)
                    headerVec.append(intVal)


                if headerVec[0] == newDict['regionID']:
                    newDict['minID'] = int(headerVec[1])
                    newDict['maxID'] = int(headerVec[2])
                    newDict['nFaces']=newDict['maxID']-newDict['minID'] + 1
                    newDict['zoneType'] = int(headerVec[3])
                    newDict['zoneTypeString'] = zoneTypeString
                    for region in self.regionList:
                        if region['name'] == ownerName:
                            newDict['faceOwner'] = int(region['regionID'])
                    if neighbourName == '':
                        newDict['faceNeighbour'] = 0
                    else:
                        for region in self.regionList:
                            if region['name'] == neighbourName:
                                newDict['faceNeighbour'] = int(
                                        region['regionID'])

                    if newDict['faceOwner'] == self.fluid_region_id:
                        newDict['fluid'] = True
                    else:
                        newDict['fluid'] = False

                    #up to this point i read only boundary walls
                    newDict['faceType'] = 'wall'


            faceList.append(newDict)

        storedFaces = [v['faceName'] for v in faceList]
        storedRegions = [v['name'] for v in self.regionList]
        stored = storedFaces + storedRegions
        maxI = max([v['origOrderID'] for v in faceList])
        i+=1

        print (storedFaces)

        filtVect = [v for v in self.headerDictCas if v['index']==39]

        for entry in filtVect:
            name = entry['header'].split()[2]

            if name  not in stored:

                ownerName = name.split('-')[-1]
                print (ownerName)

                newDict = {}
                newDict['faceName'] = name
                newDict['origOrderID'] = i
                i+=1
                newDict['regionID'] = int(entry['header'].split()[0])

                # look for the min and max ID
                filtVct1=[v for v in self.headerDictCas if v['index']==3013]
                for entry2 in filtVct1:
                    headerVec = []
                    for  val in entry2['header'].split():
                        intVal = int(val,16)
                        headerVec.append(intVal)
                    if headerVec[0] == newDict['regionID']:
                        newDict['faceType'] = 'internal'
                        newDict['minID'] = int(headerVec[1])
                        newDict['maxID'] = int(headerVec[2])
                        newDict['nFaces']= (newDict['maxID']-
                                newDict['minID'] + 1 )
                        newDict['zoneType'] = int(headerVec[3])
                        newDict['zoneTypeString'] = 'internal'

                        for region in self.regionList:
                            if region['name'] == ownerName:
                                newDict['faceOwner'] = int(
                                        region['regionID'])
                                newDict['faceNeighbour']=newDict['faceOwner']

                        if newDict['faceOwner'] == self.fluid_region_id:
                            newDict['fluid'] = True
                        else:
                            newDict['fluid'] = False

                faceList.append(newDict)


        faceList = sorted(faceList, reverse=False,
                                      key=lambda x:x['origOrderID'])

        # reorder the groups putting the internal fluidmesh at the beginning
        for face in faceList:
            if face['fluid'] == True and face['faceType'] == 'internal':
                face['newOrderID'] = 0
                face['newMinID'] =   1
                face['newMaxID'] = face['nFaces']
                idTemp = face['nFaces']
                break

        orderID = 1
        for face in faceList:
            if face['fluid'] == False :
                continue
            if face['fluid'] == True :
                if face['faceType'] == 'internal':
                    continue
                else:
                    face['newOrderID'] = orderID
                    face['newMinID'] =   idTemp + 1
                    face['newMaxID'] = idTemp + face['nFaces']
                    idTemp = face['newMaxID']
                    orderID += 1

        for face in faceList:
            if face['fluid'] == True :
                continue
            else:
                face['newOrderID'] = orderID
                face['newMinID'] =   idTemp + 1
                face['newMaxID'] = idTemp + face['nFaces']
                idTemp = face['newMaxID']
                orderID += 1

        faceList = sorted(faceList, reverse=False,
                                      key=lambda x:x['newOrderID'])



        nWaterOwner = 0

        for face in faceList:
            if face['fluid'] == True:
                nWaterOwner += face['nFaces']


        self.faceList = faceList
        self.nWaterFaces = nWaterOwner

        for dic in faceList:
            print (dic)


        return

    def getFluidFaces_h5(self):
        """
        this function reads the zone topology information relative to
        the faces of the multi block model. Using this info it redefine the
        ranges of the faces to have the internal mesh at the beginning and
        the remaining faces later. Ideally these ranges will help the rest
        of the face definition.
        """

        filename = os.path.join(self.cfd_path,self.casH5file)

        facesNamePat = re.compile(".*/faces/zoneTopology/name\Z",
                                     re.IGNORECASE)
        facesPatMin = re.compile(".*/faces/zoneTopology/minId\Z",
                                     re.IGNORECASE)

        facesPatMax = re.compile(".*/faces/zoneTopology/maxId\Z",
                                     re.IGNORECASE)

        zoneTypePat = re.compile(".*/faces/zoneTopology/zoneType\Z",
                                     re.IGNORECASE)

        faceOwnerPat = re.compile(".*/faces/zoneTopology/c0\Z",
                                     re.IGNORECASE)

        faceNeighPat = re.compile(".*/faces/zoneTopology/c1\Z",
                                     re.IGNORECASE)

        faceNames = get_h5_dataset(filename,facesNamePat)
        faceNames = list(faceNames[0].decode('UTF-8').split(';'))
        minIDs = list(get_h5_dataset(filename, facesPatMin))
        maxIDs = list(get_h5_dataset(filename, facesPatMax))
        zoneTypes = list(get_h5_dataset(filename, zoneTypePat))
        faceOwners = list(get_h5_dataset(filename, faceOwnerPat))
        faceNeighbours = list(get_h5_dataset(filename, faceNeighPat))

        sortedLists = sorted(zip(faceNames,minIDs,maxIDs,zoneTypes,
            faceOwners,faceNeighbours), key=lambda x:x[1])

        faceNames, minIDs,maxIDs,zoneTypes, faceOwners, faceNeighbours = (
                zip(*sortedLists))

        faceList = []

        for i,val in enumerate(faceNames):
            newDict = {}
            newDict['faceName'] = val
            newDict['origOrderID'] = i
            newDict['minID'] = int(minIDs[i])
            newDict['maxID'] = int(maxIDs[i])
            newDict['nFaces']=newDict['maxID'] - newDict['minID'] + 1
            newDict['origNFaces']=newDict['maxID'] - newDict['minID'] + 1
            newDict['zoneType'] = int(zoneTypes[i])
            newDict['faceOwner'] = int(faceOwners[i])
            newDict['faceNeighbour'] = int(faceNeighbours[i])
            if newDict['faceOwner'] == self.fluid_region_id:
                newDict['fluid'] = True
            else:
                newDict['fluid'] = False

            if newDict['zoneType'] == 2:
                newDict['faceType'] = 'internal'
            else:
                newDict['faceType'] = 'wall'
            faceList.append(newDict)

        faceList = sorted(faceList, reverse=False,
                                      key=lambda x:x['origOrderID'])

        # reorder the groups putting the internal fluidmesh at the beginning
        for face in faceList:
            if face['fluid'] == True and face['faceType'] == 'internal':
                face['newOrderID'] = 0
                face['newMinID'] =   1
                face['newMaxID'] = face['nFaces']
                idTemp = face['nFaces']
                break

        orderID = 1
        for face in faceList:
            if face['fluid'] == False :
                continue
            if face['fluid'] == True :
                if face['faceType'] == 'internal':
                    continue
                else:
                    face['newOrderID'] = orderID
                    face['newMinID'] =   idTemp + 1
                    face['newMaxID'] = idTemp + face['nFaces']
                    idTemp = face['newMaxID']
                    orderID += 1

        for face in faceList:
            if face['fluid'] == True :
                continue
            else:
                face['newOrderID'] = orderID
                face['newMinID'] =   idTemp + 1
                face['newMaxID'] = idTemp + face['nFaces']
                idTemp = face['newMaxID']
                orderID += 1

        #faceList = sorted(faceList, reverse=False,
        #                              key=lambda x:x['newOrderID'])

        #for dic in faceList:
        #    print (dic)


        nWaterOwner = 0

        for face in faceList:
            if face['fluid'] == True:
                nWaterOwner += face['nFaces']


        self.faceList = faceList
        self.nWaterFaces = nWaterOwner

        return

    def createHeaderDictionary(self):
        """
        this function parse all the headers in the cas/dat file and generate
        a dictionary that can be used in the other functions to get the data
        """

        filename = os.path.join(self.cfd_path,self.casfile)

        self.headerDictCas = get_fluent_parse_headers(filename)

        return

    def getFluidCells(self):
        """
        this function reads all the cells in the multi region file and
        stores those that are of the fluid type.
        """

        # get region names

        regionPat = re.compile('\(.*?\)')
        regionPat1 = re.compile('\([^()]*?\)')
        filename = os.path.join(self.cfd_path,self.casfile)
        regionNames = get_fluent_parse_regions(filename)
        self.cfdPostMeshInfoString = regionNames

        regionNames = regionPat1.findall(regionNames)

        regionNames = regionNames[0].strip('()').split()


        ## define a list of dictionary with the region info
        regionList = []
        for name in regionNames:
            newDict = {}
            newDict['name'] = name

            # look for the region ID
            filtVect = [v for v in self.headerDictCas if v['index']==39]
            for entry in filtVect:
                if entry['header'].split()[2] == name:
                    newDict['regionID'] = int(entry['header'].split()[0])

            # look for the min and max ID
            filtVect = [v for v in self.headerDictCas if v['index']==3012]
            for entry in filtVect:
                headerVec = []
                for  val in entry['header'].split():
                    intVal = int(val,16)
                    headerVec.append(intVal)


                if headerVec[0] == newDict['regionID']:
                    newDict['minID'] = int(headerVec[1])
                    newDict['maxID'] = int(headerVec[2])

            if newDict['name'] == self.fluent_fluid_region_name:
                newDict['fluid'] = True
                self.fluid_region_id = int(newDict['regionID'])
                self.fluid_cellID_min = int(newDict['minID'])
                self.fluid_cellID_max = int(newDict['maxID'])
            else:
                newDict['fluid'] = False

            regionList.append(newDict)

        self.regionList = regionList

        #for v in self.regionList:
        #    print(v)


        return

    def getFluidCells_h5(self):
        """
        this function reads all the cells in the multi region file and
        stores those that are of the fluid type.
        """

        filename = os.path.join(self.cfd_path,self.casH5file)

        regionNamePat = re.compile(".*/cells/zoneTopology/name\Z",
                                     re.IGNORECASE)

        regionIDPat = re.compile(".*/cells/zoneTopology/id\Z",
                                     re.IGNORECASE)

        partitionIdPat = re.compile(".*/cells/partition/.*/partition-ids\Z",
                                     re.IGNORECASE)

        regionPatMin = re.compile(".*/cells/zoneTopology/minId\Z",
                                     re.IGNORECASE)

        regionPatMax = re.compile(".*/cells/zoneTopology/maxId\Z",
                                     re.IGNORECASE)

        regionMinIDPaths = []
        regionMaxIDPaths = []

        regionNames = get_h5_dataset(filename,regionNamePat)
        regionNames = regionNames[0].decode('UTF-8').split(';')
        regionIDs = get_h5_dataset(filename,regionIDPat)
        minIDs = get_h5_dataset(filename, regionPatMin)
        maxIDs = get_h5_dataset(filename, regionPatMax)

        # define a list of dictionary with the region info
        regionList = []
        for i,val in enumerate(regionNames):
            newDict = {}
            newDict['name'] = val
            newDict['minID'] = int(minIDs[i])
            newDict['maxID'] = int(maxIDs[i])
            newDict['regionID'] = int(regionIDs[i])
            if val == self.fluent_fluid_region_name:
                newDict['fluid'] = True
                self.fluid_region_id = int(regionIDs[i])
                self.fluid_cellID_min = int(minIDs[i])
                self.fluid_cellID_max = int(maxIDs[i])
                self.fluid_cellN=self.fluid_cellID_max-self.fluid_cellID_min+1
            else:
                newDict['fluid'] = False

            regionList.append(newDict)

        self.regionList = regionList

        return

    def writeNut_multi_h5(self):
        """
        this file write the turbulent viscosity in the 1 folder, it is
        assumed this values is always calculated: meaning that we always
        convert from a turbulent simulation.
        """
        print ("writing multiblock openFOAM nut files ... ")

        filename = os.path.join(self.cfd_path,self.datH5file)

        mutCellPat = re.compile(".*/cells/SV_MU_T/.*",
                                     re.IGNORECASE)

        #mutCellValues = get_h5_dataset(filename,mutCellPat)
        mutCellValues = get_h5_dataset_multi(filename,mutCellPat)


        nutOneFilePath = os.path.join(self.cfd_path,'1','nut')

        filenameCAS = os.path.join(self.cfd_path,self.casH5file)

        ownerPat = re.compile(".*/faces/c0/\d+\Z", re.IGNORECASE)

        ownerList = get_h5_dataset(filenameCAS,ownerPat)

        nut1Header = """
/*------------------------------*- C++ -*---------------------------------*\\
| =========                 |                                               |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \\    /   O peration     | Version:  2112                                |
|   \\  /    A nd           | Website:  www.openfoam.com                    |
|    \\/     M anipulation  |                                               |
\*-------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       volScalarField;
    location    "1";
    object      nut;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];
"""

        faceList = sorted(self.faceList, reverse=False,
                                      key=lambda x:x['newOrderID'])

        mutValues = extract_multiblock(self.fluid_cellID_min,
                                       self.fluid_cellID_max,
                                       mutCellValues)

        with  open(nutOneFilePath,'w') as fo:

            nutStr = "{:24.18e}\n"
            fo.write(nut1Header)
            fo.write("internalField   nonuniform List<scalar>\n")
            fo.write("{:d}\n".format(self.fluid_cellN))
            fo.write("(\n")
            for val in mutValues:
                fo.write(nutStr.format(val/self.fluentDensity))

            fo.write(")\n")
            fo.write(";\n")

            fo.write("\n")
            fo.write("boundaryField\n")
            fo.write("{\n")

            for face in faceList:

                if face['faceType'] == 'internal':
                    continue
                if face['fluid'] :

                    fo.write("    {}\n".format(face["faceName"]))
                    fo.write("    {\n")
                    if face['faceType'] == 'wall':
                        fo.write("        type         nutkWallFunction;\n")
                        fo.write("        blending            stepwise;\n")
                        fo.write("        Cmu            0.09;\n")
                        fo.write("        kappa            0.41;\n")
                        fo.write("        E            9.8;\n")
                        fo.write("        value     uniform  0;\n")
                    else :
                        fo.write("        type            calculated;\n")
                        fo.write("        value   nonuniform List<scalar>\n")
                        faceOwn = ownerList[(face['minID']-1):face['maxID']]
                        fo.write("{:d}\n".format(len(faceOwn)))
                        fo.write("(\n")
                        for ow in faceOwn:
                            fo.write(nutStr.format(
                                mutValues[int(ow-self.fluid_cellID_min)]/
                                            self.fluentDensity))
                        fo.write(")\n")
                        fo.write(";\n")
                    fo.write("    }\n")

            fo.write("}\n")

        return


    def getNumCells(self):
        """this function create an attribute to the class that specifies the
        number of internal cells. It does so by reading the U file"""

        folderItems = os.listdir(self.cfd_path)

        folderTimes=[int(itm) for itm in folderItems if checkInt(itm) == True]

        lastTime = max(folderTimes)

        velocityFile = os.path.join(self.cfd_path,str(lastTime),'U')

        internalBlockPat = re.compile("internalField.*?(\d+).*?\(",
                                          re.MULTILINE | re.DOTALL )

        try:
            inpFile = open(velocityFile,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open velocity file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            cellNumber = internalBlockPat.findall(text)[0]

        self.num_internal_cells = int(cellNumber)

        return


    def readDensity_h5(self):
        """
        this function reads the density
        """

        filename = os.path.join(self.cfd_path,self.datH5file)

        denCellPat = re.compile(".*/cells/SV_DENSITY/.*",
                                     re.IGNORECASE)


        cell_density_path = []


        with h5py.File(filename, "r") as fi:

            paths = get_dataset_keys(fi)

            for path in paths:
                denMatches = denCellPat.findall(path)
                if len(denMatches) == 1:
                    cell_density_path.append(path)

            if len(cell_density_path) == 1:
                self.denCellPath = cell_density_path[0]
                if hasattr(self,'minIDregion'):
                    self.fluentDensity = (fi[self.denCellPath]
                                            [self.minIDregion - 1])
                else:
                    self.fluentDensity = fi[self.denCellPath][0]
            else:
                print ("ERROR zero or more than one density datasets found")
                sys.exit()

        return



def main():



    parser = argparse.ArgumentParser(description = "FLUNED case generator")
    parser.add_argument('-i','--input', type=str, help="input")
    parser.add_argument('-l',"--launch_simulation", action ='store_true' ,
            help="launch simulation", default = False)
    parser.add_argument('-t',"--input_template", action ='store_true' ,
            help="create an input template", default = False)
    args = parser.parse_args()

    if not args.input and not args.input_template:
        print ("WARNING no input provided")
        print ("printing template and exiting")
        create_input_template()
        sys.exit()

    if args.input_template:
        print ("printing template and exiting")
        create_input_template()
        sys.exit()

    inputCases = read_input_file(args.input)

    # define a vector of FLUNED cases

    for case in inputCases:

        print ("creating FLUNED case...")

        fCase = FlunedCase(case)

        if fCase.cfd_type == 'fluent-h5-multi':
            print ("parsing h5 files ... ")
            fCase.getH5files()
            fCase.create_case_folders_fluent()
            fCase.getFluidCells_h5()
            fCase.getFluidFaces_h5()
            fCase.readDensity_h5()
            fCase.defineWalls_multi_h5()
            fCase.writeOwner_multi_h5()
            fCase.writeNeighbour_multi_h5()
            fCase.write_cell_zones()
            fCase.writeBoundary_multi_h5()
            fCase.writeFaces_multi_h5()
            fCase.writeNodes_multi_h5()
            fCase.writePhi_multi_h5()
            fCase.writeSpeed_multi_h5()
            fCase.writeNut_multi_h5()

        if fCase.cfd_type == 'fluent-multi':
            print ("parsing binary cas/dat files not implemented yet ... ")
            sys.exit()
            #fCase.getCASDATfiles()
            #fCase.createCaseFoldersFluent()
            #fCase.createHeaderDictionary()
            #fCase.getFluidCells()
            #fCase.getFluidFaces()
            #fCase.writeNodes_multi()


        print ("copying last CFD iteration files ... ")
        fCase.create_case_folder()
        fCase.copy_last_phi()
        fCase.copyLastU()
        fCase.getNumCells()
        fCase.copyLastNut()
        fCase.copyPolyMesh()
        fCase.reconstructFaces()
        fCase.generate_system_file()
        fCase.launchVolFuncObjects()
        fCase.generateConstantFiles()
        fCase.generate_zero_t()
        fCase.generate_zero_ta()
        fCase.generate_zero_td()
        fCase.generate_zero_tr()
        fCase.generateSourceFile()
        fCase.generateTrSourceFile()

        if args.launch_simulation:

            fCase.launch_solver()

        print ("FINISHED!")



    return

if __name__ == '__main__':
    main()
