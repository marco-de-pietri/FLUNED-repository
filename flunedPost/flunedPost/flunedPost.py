# -*- coding: utf-8 -*-
import re
import sys
import os
import copy
#import matplotlib.pyplot as plt
import math
import numpy as np
import subprocess
import argparse
import linecache
import tracemalloc
import vtk
import pyvista as pv
from vtk.util import numpy_support as VN

def mergeContinueRuns(time_lists,data_lists):
    """
    this function takes time series segments and join them in a single one -
    the overlapping sections are removed
    """

    if len(time_lists) == 0 or len(data_lists) == 0 :
        print ("ERROR could not read post processing data")
        sys.exit()

    if len(time_lists) != len(data_lists) :
        print ("ERROR mismatch in length of data series")
        print ("number of time lists: ", len(time_lists))
        print ("number of data lists: ", len(data_lists))
        sys.exit()

    tot_len_time_lists = sum([len(x) for x in time_lists])
    tot_len_data_lists = sum([len(x) for x in data_lists])

    if tot_len_time_lists != tot_len_data_lists:
        print ("ERROR mismatch in length of data series")
        print ("number of time points: ", tot_len_time_lists)
        print ("number of data points: ", tot_len_data_lists)

    temp_list = list(zip(time_lists,data_lists))

    #time_lists_sorted = sorted(time_lists,
    #                           key=lambda x:x[0])

    sorted_temp_list = sorted(temp_list, key= lambda x:x[0][0])

    time_lists_sorted,data_lists_sorted = zip(*sorted_temp_list)
    #data_lists_sorted = sorted(data_lists,
    #                           key=lambda x:x[0])

    time_series = time_lists_sorted[0]
    data_series = data_lists_sorted[0]


    if len(time_lists) > 1:
        for i,l in enumerate(time_lists_sorted):
            time_series_temp = copy.deepcopy(time_series)
            if i == 0:
                continue
            appended = False
            for j, val in enumerate(time_series_temp):
                if val >= time_lists_sorted[i][0]:
                    time_series = time_series[0:j]
                    time_series.extend(l)
                    data_series = data_series[0:j]
                    data_series.extend(data_lists_sorted[i])
                    appended = True
                    break

            if not appended:
                time_series.extend(l)
                data_series.extend(data_lists_sorted[i])

    #for time,data in zip(time_series,data_series):
    #    print (f"{time:.0f}")
    #    print (f"{data:.0e}")

    return time_series, data_series

def getPostFiles(fPath):
    """
    this function crawl the folder to reach the data in the post process file
    """

    filePaths = []
    folderItms = os.listdir(fPath)
    folderItms = sorted(folderItms, reverse=False,key=lambda x:float(x))

    for fld in folderItms:

        fPath1 = os.path.join(fPath,fld)
        if os.path.isdir(fPath1):
            fileItms = os.listdir(fPath1)
            completePath = os.path.join(fPath1,fileItms[0])
            filePaths.append(completePath)

    return filePaths

def postFileArea(fPath):
    """
    this function extracts the face area (in m2) contained in the post
    processing file
    """


    try:
        postFile = open(fPath,'r',encoding="utf8", errors='ignore')
    except IOError:
        print("couldn't open postprocess file")
        sys.exit()
    with postFile:
        lines = postFile.readlines()
        for line in lines:
            line = line.replace('#','')
            wrds = line.split()
            if 'area' in wrds[0].lower():
                area = float(wrds[-1])
                break

    return area

def postFileArray(fPathList,name):
    """
    this function extracts a generic array contained in the post
    processing file
    """



    arrayList = []

    for fPath in fPathList:

        array = []
        index = -1

        try:
            postFile = open(fPath,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open postprocess file")
            sys.exit()
        with postFile:
            lines = postFile.readlines()
            for line in lines:

                line = line.replace('#','')

                wrds = line.split()

                if len(wrds) == 0:
                    continue

                if 'time' in wrds[0].lower():
                    index = wrds.index(name)
                    continue

                if index != -1:
                    array.append(float(wrds[index]))

            arrayList.append(array)


    return arrayList

def getVTKsize(vtkFile):
    """this function reads the unstructured mesh and get the cartesian
    boundaries"""

    checkFile = os.path.isfile(vtkFile)
    if not checkFile:
        print ("ERROR vtk file not found")
        sys.exit()


    #read the vtk file with an unstructured grid
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtkFile)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    # return is in meter (for the standard openfoam vtk file)
    bounds =list(data.GetBounds())
    return bounds

def getTotalAtomsVtk(vtkFile, datasetName):
    """this function reads the vtk and calulates the total conc"""

    checkFile = os.path.isfile(vtkFile)
    if not checkFile:
        print ("ERROR vtk file not found")
        sys.exit()


    #read the vtk file with an unstructured grid
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtkFile)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()

    elementConc = data.GetCellData().GetArray(datasetName)
    elementConcArray = VN.vtk_to_numpy(elementConc)



    #pyvista part
    mesh = pv.read(vtkFile)
    mesh = mesh.compute_cell_sizes()
    volume_data =  mesh.cell_data['Volume']

    totAtom=sum([conc*vol for conc,vol in zip(elementConcArray,volume_data)])


    return totAtom

def sampleCoordinatesVTK(vtkFile, datasetName, coordinates):
    """this function reads the vtk and sample the reaction rates"""

    checkFile = os.path.isfile(vtkFile)
    if not checkFile:
        print ("ERROR vtk file not found")
        sys.exit()

    print ("Sampling vtk ... ")

    #read the vtk file with an unstructured grid
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtkFile)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()



    #define probe

    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(coordinates))

    for i,val in enumerate(coordinates):
        points.SetPoint(i,val[0],val[1],val[2])



    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)


    #Perform the interpolation
    probeFilter =vtk.vtkProbeFilter()
    probeFilter.SetSourceData(data)
    probeFilter.SetInputData(polydata)
    probeFilter.Update()

    vtkArray = probeFilter.GetOutput().GetPointData().GetArray(datasetName)

    concentrations = VN.vtk_to_numpy(vtkArray)

    for i in range(len(concentrations)): # remove eventual tiny negative conc
        if concentrations[i] < 0 :
            concentrations[i] = 0




    return  concentrations

def formatValues(vector):
    maxLen = 70
    returnString = ''
    newLine = ''

    for item in vector:
        newNumber = '{:.7e}'.format(item)
        if returnString == '' and newLine == '':
            newLine = newNumber
            continue

        if len(newLine + ' ' + newNumber) > maxLen :
            returnString +=newLine + '\n'
            newLine = newNumber

        else:
            newLine = newLine + ' ' + newNumber

    returnString += newLine + '\n'

    return returnString

def display_top(snapshot, key_type='lineno', limit=10):
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        print("#%s: %s:%s: %.1f KiB"
              % (index, frame.filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))

    return


def heronFormula(points):

    dist1 = pointDistance(points[0],points[1])
    dist2 = pointDistance(points[1],points[2])
    dist3 = pointDistance(points[2],points[0])

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

def pointDistance(point1,point2):

    dist = pow((sum([pow((point1[i] - point2[i]),2) for i in [0,1,2]])),0.5)

    return dist

def mod(vec1):

    sum = 0
    for el in vec1:
        sum += el**2

    return pow(sum,0.5)

def calculateArea(points):
    """ this function calculate the area of the element surface - these can
    be tetra with a triangular face or wedge with a quadrilateral face """

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

def check_float(str):
    """
    this function checks if a string is a float
    """
    try:
        float(str)
        return True
    except ValueError:
        return False

class flunedCase:
    """
    fluned case class
    """
    def __init__(self, path):
        """
        initialize case and create FLUNED case folder
        """

        self.fluned_path = os.path.join(path)
        self.case = os.path.split(path)[1]
        self.n_elements = 0


    def readPostProcess_flows(self):
        """
        this function reads the post process folders and extract from all
        the folders the flows, area and type of the faces
        """


        print ("reconstructing interface properties...")

        postFolder = os.path.join(self.fluned_path,'postProcessing')

        # read the files in the postprocess folder
        postFlowVec = []

        folderItms = os.listdir(postFolder)
        flowFolders = [itm for itm in folderItms if itm[0:8]=='volFlow-']

        for itm in flowFolders:
            postDic = {}
            postDic['faceID'] = itm[8:]
            postDic['folder'] = itm
            postDic['flowFiles'] = getPostFiles(os.path.join(postFolder,itm))
            postDic['areaFile'] = postFileArea(postDic['flowFiles'][0])

            timeLists = postFileArray(postDic['flowFiles'],
                                         'Time')

            flowLists =  postFileArray(postDic['flowFiles'],
                                       'sum(phi)')

            timeListSorted,flowListSorted = mergeContinueRuns (
                                            timeLists,
                                            flowLists,
                                            )
            postDic['timeFile']    = timeListSorted
            postDic['phiFlowFile'] = flowListSorted
            postDic['phiFlowFileLast'] = postDic['phiFlowFile'][-1]

            if postDic['phiFlowFileLast'] > 0:
                postDic['typeFile'] = 'outlet'
            elif postDic['phiFlowFileLast'] == 0:
                postDic['typeFile'] = 'wall'
            elif postDic['phiFlowFileLast'] < 0:
                postDic['typeFile'] = 'inlet'
            #print (postDic)
            postFlowVec.append(postDic)


        self.facesPost = postFlowVec



        return

    def readPostProcess_Trflows(self):
        """
        this function reads the post process folders and extract from all
        the folders the flow of the Tr scalar
        """



        postFolder = os.path.join(self.fluned_path,'postProcessing')

        # read the files in the postprocess folder
        postTrFlowVec = []

        folderItms = os.listdir(postFolder)
        flowFolders = [itm for itm in folderItms if itm[0:10]=='volTrFlow-']

        if len(flowFolders) == 0:
            # if the output relative to time is are not present
            # initialize to zero the following variables
            for face in self.facesPost:
                face['avTrFile']   =  [0]
                face['avTrGrad']   =  [0]
                face['avTrFrac']   =  0
                face['rtdResTime'] =  0
                face['rtdDecRate'] =  0
            return

        for itm in flowFolders:
            postDic = {}
            postDic['faceID'] = itm[10:]
            postDic['folder'] = itm
            postDic['flowFiles'] = getPostFiles(os.path.join(postFolder,itm))

            timeLists = postFileArray(postDic['flowFiles'],
                                         'Time')

            flowLists =  postFileArray(postDic['flowFiles'],
                                       'sum(Tr)')

            timeListSorted,flowListSorted = mergeContinueRuns (
                                            timeLists,
                                            flowLists,
                                            )

            postDic['TrFlowFile'] = flowListSorted
            #postDic['TrFlowFile']=postFileArray(
            #        postDic['flowFiles'],'sum(Tr)')
            postDic['TrFlowFileLast'] = postDic['TrFlowFile'][-1]
            #print (postDic)
            postTrFlowVec.append(postDic)


        for face in self.facesPost:
            face['TrFlowFileLast'] = ([v['TrFlowFileLast'] for
                                       v in postTrFlowVec
                                       if v['faceID']==face['faceID']][0])
            face['TrFlowFile'] = ([v['TrFlowFile'] for
                                       v in postTrFlowVec
                                       if v['faceID']==face['faceID']][0])
            face['avTrFileLast'] = (face['TrFlowFileLast']
                                   / face['phiFlowFileLast'])
            face['avTrFile'] = ([t/f for t,f in
                                zip(face['TrFlowFile'],face['phiFlowFile'])])
            face['avTrGrad'] = (list((np.gradient(face['avTrFile'],
                                     face['timeFile'],
                                     edge_order=1))))


            dt = face['timeFile'][1] - face['timeFile'][0]

            face['avTrFrac'] = [dt*g for g in face['avTrGrad']]

            face['rtdResTime'] = sum([t*g for t,g in zip(face['timeFile'],
                                    face['avTrFrac'])])

            face['rtdDecRate'] = sum([g*math.exp(-self.decay_constant*t)
                                     for t,g in zip(face['timeFile'],
                                     face['avTrFrac'])])





        return

    def readPostProcess_Tflows(self):
        """
        this function reads the post process folders and extract from all
        the folders the flow of the T scalar
        """



        postFolder = os.path.join(self.fluned_path,'postProcessing')

        # read the files in the postprocess folder
        postTFlowVec = []

        folderItms = os.listdir(postFolder)
        flowFolders = [itm for itm in folderItms if itm[0:9]=='volTFlow-']

        for itm in flowFolders:
            postDic = {}
            postDic['faceID'] = itm[9:]
            postDic['folder'] = itm
            postDic['flowFiles'] = getPostFiles(os.path.join(postFolder,itm))
            timeLists = postFileArray(postDic['flowFiles'],
                                         'Time')

            flowLists =  postFileArray(postDic['flowFiles'],
                                       'sum(T)')

            timeListSorted,flowListSorted = mergeContinueRuns (
                                            timeLists,
                                            flowLists,
                                            )
            postDic['timeFile']    = timeListSorted
            postDic['TFlowFile'] = flowListSorted
            #postDic['TFlowFile']=postFileArray(postDic['flowFiles'],'sum(T)')
            postDic['TFlowFileLast'] = postDic['TFlowFile'][-1]
            #print (postDic)
            postTFlowVec.append(postDic)


        for face in self.facesPost:
            face['TFlowFileLast'] = ([v['TFlowFileLast'] for
                                      v in postTFlowVec if
                                      v['faceID']==face['faceID']][0])
            face['TFlowFile'] = ([v['TFlowFile'] for v in postTFlowVec if
                                      v['faceID']==face['faceID']][0])
            face['avTFileLast']=face['TFlowFileLast']/face['phiFlowFileLast']
            face['avTFile'] = ([t/f for t,f in
                                zip(face['TFlowFile'],face['phiFlowFile'])])

        return



    def readDefVector(self,elID):
        """ this function looks for the studied element - and output the
        definition - no storing"""


        polyMeshFolder = os.path.join(self.fluned_path,'constant','polyMesh')
        facefilename = os.path.join(polyMeshFolder,'faces')
        try:
            inpFile = open(facefilename,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open face file")
            sys.exit()
        with inpFile:

            lineStartPat = re.compile("\A\s*\(")
            defPat = re.compile("(?<=\().*(?=\))")

            start = False
            startLine = 0

            targetLine = elID

            for i,line in enumerate(inpFile):
                if start == False:
                    if len(lineStartPat.findall(line)) != 0:
                        startLine = i
                        targetLine = startLine + elID + 1
                        start = True

                else:
                    if i == targetLine:
                        defLine = line
                        break

            definition = defPat.findall(defLine)
            defVector = [int(num) for num in definition[0].split()]

        return defVector

    def readPointDefinitions(self,pointIDVector):
        """ this function looks for the studied points - and output the
        definition - no storing"""

        #add test for empty lines - later


        returnDic = {}

        polyMeshFolder = os.path.join(self.fluned_path,'constant','polyMesh')
        pointfilename = os.path.join(polyMeshFolder,'points')
        try:
            inpFile=open(pointfilename,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open point file")
            sys.exit()
        with inpFile:

            lineStartPat = re.compile("\A\s*\(")
            defPat = re.compile("(?<=\().*(?=\))")

            start = False
            startLine = 0

            targetLines = sorted(pointIDVector)

            for i,line in enumerate(inpFile):
                if start == False:
                    if len(lineStartPat.findall(line)) != 0:
                        startLine = i
                        targetLines=([(val+startLine+1) for
                            val in targetLines])
                        start = True

                else:
                    if i in targetLines:
                        defLine = line
                        definition = defPat.findall(defLine)
                        defVector=([float(num) for num in
                            definition[0].split()])
                        returnDic[i-startLine-1] = defVector


        return returnDic

    def readElDefinitions(self,elIDVector):
        """ this function looks for the studied element - and output the
        definition - no storing"""

        #add test for empty lines - later


        returnDic = {}

        polyMeshFolder = os.path.join(self.fluned_path,'constant','polyMesh')
        facefilename = os.path.join(polyMeshFolder,'faces')
        try:
            inpFile = open(facefilename,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open face file")
            sys.exit()
        with inpFile:

            lineStartPat = re.compile("\A\s*\(")
            defPat = re.compile("(?<=\().*(?=\))")

            start = False
            startLine = 0

            targetLines = sorted(elIDVector)

            for i,line in enumerate(inpFile):
                if start == False:
                    if len(lineStartPat.findall(line)) != 0:
                        startLine = i
                        targetLines=([(val+startLine+1) for
                            val in targetLines])
                        start = True

                else:
                    if i in targetLines:
                        defLine = line
                        definition = defPat.findall(defLine)
                        defVector=[int(num) for num in definition[0].split()]
                        returnDic[i - startLine -1] = defVector


        return returnDic

    def readPointCoordinate(self,pointID):
        """ this function looks for the studied point - and output the
        coordinates - no storing"""


        polyMeshFolder = os.path.join(self.fluned_path,'constant','polyMesh')
        pointfilename = os.path.join(polyMeshFolder,'points')
        try:
            inpFile = open(pointfilename,'r',
                    encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open points file")
            sys.exit()
        with inpFile:

            lineStartPat = re.compile("\A\s*\(")
            defPat = re.compile("(?<=\().*(?=\))")

            start = False
            startLine = 0

            targetLine = pointID

            for i,line in enumerate(inpFile):
                if start == False:
                    if len(lineStartPat.findall(line)) != 0:
                        startLine = i
                        targetLine = startLine + pointID + 1
                        start = True

                else:
                    if i == targetLine:
                        defLine = line
                        break

            definition = defPat.findall(defLine)
            defVector = [float(num) for num in definition[0].split()]

        return defVector






    def readVolumes(self):
        """ this function reads the volumes from the V file located in the
        zero folder """


        print ("reading volume values...")

        # common patterns
        internalBlockPat = re.compile("internalField.*?\((.{1,}?)\)",
                                          re.MULTILINE | re.DOTALL )

        nElPat = re.compile("internalField.*?(\d+).*?\(",
                                          re.MULTILINE | re.DOTALL )
        vFile = os.path.join(self.fluned_path,'0','V')
        try:
            inpFile = open(vFile,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open Volume V file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            cellNumberText = nElPat.findall(text)
            self.n_elements = (int(cellNumberText[0]))
            numInternalBlocks = internalBlockPat.findall(text)
            internalVolumes = numInternalBlocks[0].split('\n')[1:-1]
            self.Volumes = np.zeros(self.n_elements)
            for i in range(self.n_elements):
                self.Volumes[i] = float(internalVolumes[i])

        return


    def readVelocities(self):
        """ this function reads the velocity from the U file located in the
        zero folder """

        print ("reading velocity values...")

        # common patterns
        internalBlockPat = re.compile("internalField.*?\((.{1,}?)\n\\s*\)",
                                          re.MULTILINE | re.DOTALL )

        uFile = os.path.join(self.fluned_path,'0','U')
        try:
            inpFile = open(uFile,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open  U file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            numInternalBlocks = internalBlockPat.findall(text)
            internalVels = numInternalBlocks[0].split('\n')[1:]
            internalVels=[val.strip('()') for val in internalVels]

            self.Velocities = np.array(
                              [np.array(
                              [float(val) for val in v.split()])
                                for v in internalVels])


        return

    def readGradT(self):
        """ this function reads the T gradient"""

        print ("reading gradient values...")

        # find the last folder
        #folderItms = os.listdir(self.fluned_path)
        #folderTimes=[int(itm) for itm in folderItms if checkInt(itm) == True]
        #lastTime = max(folderTimes)

        # common patterns
        internalBlockPat = re.compile("internalField.*?\((.{1,}?)\n\\s*\)",
                                          re.MULTILINE | re.DOTALL )

        gradFile = os.path.join(self.last_time_step,'grad(T)')
        try:
            inpFile = open(gradFile,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open grad file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            numInternalBlocks = internalBlockPat.findall(text)
            internalVels = numInternalBlocks[0].split('\n')[1:]
            internalVels=[val.strip('()') for val in internalVels]
            self.Gradients = np.array([
                             np.array([float(val) for val in v.split()])
                                for v in internalVels])


        return
    def readT(self):
        """
        this function reads the T scalar
        """


        print ("reading scalar values...")


        # find the last folder
        #folderItms = os.listdir(self.fluned_path)
        #folderTimes=[int(itm) for itm in folderItms if checkInt(itm) == True]
        #lastTime = max(folderTimes)

        # common patterns
        internalBlockPat = re.compile("internalField.*?\((.{1,}?)\)",
                                          re.MULTILINE | re.DOTALL )

        tFile = os.path.join(self.last_time_step,'T')

        try:
            inpFile = open(tFile,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open Volume T file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            numInternalBlocks = internalBlockPat.findall(text)
            internalScalar = numInternalBlocks[0].split('\n')[1:-1]
            self.TScalar = np.zeros(self.n_elements)
            for i in range(self.n_elements):
                self.TScalar[i] = float(internalScalar[i])

        return

    def getTimeTreatment(self):
        """
        this function parse the fvSchemes file to check if the simulation
        is steady state or transient
        """

        # common patterns
        ddtPat = re.compile("ddtSchemes.*?\{.{1,}?\}",
                                          re.MULTILINE | re.DOTALL )

        cFile = os.path.join(
                self.fluned_path,'system','fvSchemes')
        try:
            inpFile = open(cFile,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open fvSchemes file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            ddtText = ddtPat.findall(text)[0]
            if 'euler'  in ddtText.lower():
                self.time_treatment = 'transient'
            elif 'steadystate'  in ddtText.lower():
                self.time_treatment = 'steadystate'
            else:
                print("couldn't open recognise the time scheme")
                sys.exit()


        return

    def parseConstants(self):
        """ this function parses the constant properties to get the decay
        variable and the others"""

        # common patterns
        dtPat = re.compile("DT\s*DT.*")
        lambdaPat = re.compile("lambda\s*lambda.*")
        schPat = re.compile("Sct\s*Sct.*")

        cFile = os.path.join(
                self.fluned_path,'constant','transportProperties')
        try:
            inpFile = open(cFile,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open transportProperties file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            dtLines = dtPat.findall(text)
            lambdaLines = lambdaPat.findall(text)
            schLines = schPat.findall(text)
            if len(dtLines) != 0:
                vals = dtLines[0].strip(' ;').split()
                val = vals[-1]
                self.molecular_dif = float(val)
            else:
                self.molecular_dif = 0

            if len(lambdaLines) != 0:
                vals = lambdaLines[0].strip(' ;').split()
                val = vals[-1]
                self.decay_constant = float(val)
            else:
                self.decay_constant = 0

            if len(schLines) != 0:
                vals = schLines[0].strip(' ;').split()
                val = vals[-1]
                self.schmidt_number = float(val)
            else:
                self.schmidt_number = 0

        return

    def write_summary(self,arguments):
        """
        write a bunch of results of the simulation
        """


        print ("printing final Summary ... ")


        if self.time_treatment == 'steadystate':
            self.write_summary_steady(arguments)
        elif self.time_treatment == 'transient':
            self.write_summary_transient(arguments)

        return

    def write_summary_transient(self,arguments):

        resFolder = os.path.join(self.fluned_path,'RESULTS')

        dirCheck = os.path.isdir(resFolder)
        if not dirCheck:
            os.mkdir(resFolder)

        sumFile = os.path.join(resFolder,'SUMMARY.csv')

        inletAtoms = 0
        inletFlow = 0
        outletAtoms = 0
        outletFlow = 0

        for face in self.facesPost:

            if face['typeFile'] == 'inlet':
                inletAtoms += abs(face['TFlowFileLast'])
                inletFlow += abs(face['phiFlowFileLast'])

            if face['typeFile'] == 'outlet':
                outletAtoms += abs(face['TFlowFileLast'])
                outletFlow += abs(face['phiFlowFileLast'])

        totalInletConc = inletAtoms/inletFlow
        totalInActivity = totalInletConc*self.decay_constant

        totalOutletConc = outletAtoms/outletFlow
        totalOutActivity = totalOutletConc*self.decay_constant

        cellActivity = np.zeros(len(self.TScalar))

        for i in range(len(self.TScalar)):
            if self.TScalar[i] > 0:
                cellActivity[i] = (self.TScalar[i]*
                                   self.Volumes[i]*
                                   self.decay_constant)
            else:
                cellActivity[i] = 0

        totActivity = sum(cellActivity)
        avgActivity = totActivity/sum(self.Volumes)

        averageVolume = sum(self.Volumes)/len(self.Volumes)

        negCheck = any(n < 0 for n in self.TScalar)

        with open(sumFile,'w') as fw:
            fw.write("FLUNED SIMULATION SUMMARY\n")
            fw.write("CASE,{},\n".format(self.case))
            fw.write("TRANSIENT SIMULATION\n")
            fw.write("N ELEMENTS,{},\n".format(len(self.TScalar)))
            fw.write("ISOTOPE,{},\n".format(self.isotope))
            fw.write("DECAY CONSTANT,{:e},\n".format(self.decay_constant))
            fw.write("MOL DIFFUSION,{:e},\n".format(self.molecular_dif))
            fw.write("TURB SCHMIDT N,{:f},\n".format(self.schmidt_number))
            fw.write("\n")
            fw.write("\n")
            fw.write("QUALITY\n")
            fw.write("AVG VOL [m3],{:e},\n".format(averageVolume))

            if negCheck:
                fw.write("WARNING some elements are negative")

            if arguments.source_CDGS:
                fw.write("\n")
                fw.write("\n")
                fw.write("SOURCE SAMPLING\n")
                fw.write("SAMPLING RESOLUTION [m],{:f},\n".format(
                    self.precision))
                fw.write("SAMPLED VOXELS [#],{:d},\n".format(self.xInts*
                                                             self.yInts*
                                                             self.zInts))
                fw.write("VTK EMISSION RATE [#/s],{:e},\n".format(
                    self.originalEmissionRate))
                sampString = "SAMPLED EMISSION RATE (UNSCALED) [#/s],{:e},\n"
                fw.write(sampString.format(self.unscaledEmissionRate))

            fw.write("\n")
            fw.write("\n")
            fw.write("ACTIVATION\n")
            fw.write("INLET ATOMS FINAL [#/s],{:.5e},\n".format(inletAtoms))
            fw.write("OUTLET ATOMS FINAL [#/s],{:.5e},\n".
                    format(outletAtoms))
            fw.write("TOT IN ACTIVITY FINAL [Bq/m3],{:.5e},\n".format(
                totalInActivity))
            fw.write("TOT OUT ACTIVITY FINAL [Bq/m3],{:.5e},\n".format(
                totalOutActivity))
            fw.write("TOT CASE ACTIVITY FINAL [Bq],{:.5e},\n".format(
                totActivity))
            fw.write("TOT AVG ACTIVITY FINAL [Bq/m3],{:.5e},\n".format(
                avgActivity))
            if totalInActivity != 0:
                normAvgActivity = avgActivity/totalInActivity
                stri="INLET-NORMALIZED AVG ACTIVITY FINAL [Bq/m3],{:.5e},\n"
                fw.write(stri.format(normAvgActivity))





            if inletAtoms != 0:
                reductionRate = outletAtoms/inletAtoms
                fw.write("OUT/IN RATIO FINAL,{:.5e},\n".
                        format(reductionRate))


            fw.write("\n")
            fw.write("\n")
            fw.write("FACES\n")

            for face in self.facesPost:
                if face['typeFile'] != 'wall':
                    fw.write(face['faceID'])
                    fw.write("\n")
                    fw.write("TYPE,{},\n".format(face['typeFile']))
                    fw.write("AREA [m2],{:.5e},\n".format(face['areaFile']))
                    fw.write("FLUID FLOW FINAL [m3/s],{:.5e},\n".format(
                        face['phiFlowFileLast']))
                    fw.write("ATOM FLOW FINAL [#/s],{:.5e},\n".format(
                        face['TFlowFileLast']))
                    fw.write("ATOM CONC FINAL [#/m3],{:.5e},\n".format(
                        face['avTFileLast']))
                    fw.write("SPECIFIC ACTIVITY FINAL [Bq/m3],{:.5e},\n".
                            format(abs(face['avTFileLast'])*
                            self.decay_constant))
                    if face['typeFile'] == 'outlet':
                        fluxFrac = abs(face['phiFlowFileLast']/outletFlow)
                        fw.write("AVG RTD RES T [s],{:.5f},\n".format(
                                     abs(face['rtdResTime'])))
                        fw.write("RTD FACE RED RATE,{:.5f},\n".format(
                                     fluxFrac*(face['rtdDecRate'])))
                    fw.write("\n")

        for face in self.facesPost:
            if face['typeFile'] == 'wall':
                continue

            # write the flowing atoms at inlet-outlet
            faceSummary=('face-atom-flow-' + face['faceID'] +
                         '-' + face['typeFile'] + '.csv')
            sumFile1 = os.path.join(resFolder,faceSummary)
            with open(sumFile1,'w') as fw:
                for t,c in zip(face['timeFile'],face['TFlowFile']):
                    fw.write('{:.3f},{:.5e},\n'.format(t,c))

            # write the concentration at inlet-outlet
            faceSummary=('face-conc-' + face['faceID'] +
                         '-' + face['typeFile'] + '.csv')
            sumFile1 = os.path.join(resFolder,faceSummary)
            with open(sumFile1,'w') as fw:
                for t,c in zip(face['timeFile'],face['avTFile']):
                    fw.write('{:.3f},{:.5e},\n'.format(t,c))

            # write the specific activity at inlet-outlet
            faceSummary=('face-specific-activity-' + face['faceID'] +
                         '-' + face['typeFile'] + '.csv')
            sumFile1 = os.path.join(resFolder,faceSummary)
            with open(sumFile1,'w') as fw:
                for t,c in zip(face['timeFile'],face['avTFile']):
                    fw.write('{:.3f},{:.5e},\n'.format(
                             t,c*self.decay_constant))

            # write the specific time-conc at inlet-outlet
            faceSummary=('face-fictitious-time-' + face['faceID'] +
                         '-' + face['typeFile'] + '.csv')
            sumFile1 = os.path.join(resFolder,faceSummary)
            with open(sumFile1,'w') as fw:
                for t,c in zip(face['timeFile'],face['avTrFile']):
                    fw.write('{:.3f},{:.5e},\n'.format(t,c))


            # write the RTD at inlet-outlet
            faceSummary=('face-RTD-raw-' + face['faceID'] +
                         '-' + face['typeFile'] + '.csv')
            sumFile1 = os.path.join(resFolder,faceSummary)
            with open(sumFile1,'w') as fw:
                for t,g in zip(face['timeFile'],face['avTrGrad']):
                    fw.write('{:.3f},{:.5e},\n'.format(t,g))


        return


    def write_summary_steady(self,arguments):
        """
        This function writes the summary file in the RESULTS/ folder
        """

        result_dir = os.path.join(self.fluned_path,'RESULTS')

        dir_check = os.path.isdir(result_dir)
        if not dir_check:
            os.mkdir(result_dir)

        summary_file = os.path.join(result_dir,'SUMMARY.csv')

        inletAtoms = 0
        inletFlow = 0
        outletAtoms = 0
        outletFlow = 0

        for face in self.facesPost:

            if face['typeFile'] == 'inlet':
                inletAtoms += abs(face['TFlowFileLast'])
                inletFlow += abs(face['phiFlowFileLast'])

            if face['typeFile'] == 'outlet':
                outletAtoms += abs(face['TFlowFileLast'])
                outletFlow += abs(face['phiFlowFileLast'])

        totalInletConc = inletAtoms/inletFlow
        totalInActivity = totalInletConc*self.decay_constant

        totalOutletConc = outletAtoms/outletFlow
        totalOutActivity = totalOutletConc*self.decay_constant

        cellActivity = np.zeros(len(self.TScalar))

        for i,_ in enumerate(self.TScalar):
            if self.TScalar[i] > 0:
                cellActivity[i] = (self.TScalar[i]*
                                   self.Volumes[i]*
                                   self.decay_constant)
            else:
                cellActivity[i] = 0

        totActivity = sum(cellActivity)
        avgActivity = totActivity/sum(self.Volumes)

        averageVolume = sum(self.Volumes)/len(self.Volumes)



        #parameter 1 based on residence times
        #averageCellReTime = sum(self.resTimes)/len(self.resTimes)

        #qualityPar = [t*c*lbda/(ln2*totS) for t,c
        #              in zip(self.resTimes,self.TScalar)]

        #avgMeshQualParameter = sum(qualityPar)/len(qualityPar)




        if arguments.check:

            lbda = self.decay_constant
            ln2 = np.log(2)

            totS = sum(self.TScalar)

            #parameter 2 based on transit times

            transitTimes = [(vol**(1/3))/mod(v) for vol,v in
                            zip(self.Volumes,self.Velocities)]

            averageTransitTime = sum(transitTimes)/len(transitTimes)

            qualityPar2 = [t*c*lbda/(ln2*totS) for t,c
                          in zip(transitTimes,self.TScalar)]

            avgMeshQualParameter2 = sum(qualityPar2)/len(qualityPar2)

            #parameter 3 based on scalar gradient

            gradients = [mod(g) for g in self.Gradients]

            avgGrad = sum(gradients)/len(gradients)

            gradientDist = [g*(vol**(1/3)) for vol,g in
                            zip(self.Volumes,gradients)]

            avgGradDist = sum(gradientDist)/len(gradientDist)

            qualityPar3 = [gd*c/totS for gd,c in
                           zip(gradientDist,self.TScalar)]

            avgMeshQualParameter3 = sum(qualityPar3)/len(qualityPar3)

        negCheck = any(n < 0 for n in self.TScalar)

        with open(summary_file,'w') as fw:
            fw.write("FLUNED SIMULATION SUMMARY\n")
            fw.write("CASE,{},\n".format(self.case))
            fw.write("STEADY STATE SIMULATION\n")
            fw.write("N ELEMENTS,{},\n".format(len(self.TScalar)))
            fw.write("ISOTOPE,{},\n".format(self.isotope))
            fw.write("DECAY CONSTANT,{:e},\n".format(self.decay_constant))
            fw.write("MOL DIFFUSION,{:e},\n".format(self.molecular_dif))
            fw.write("TURB SCHMIDT N,{:f},\n".format(self.schmidt_number))
            fw.write("\n")
            fw.write("\n")
            fw.write("QUALITY\n")
            fw.write("AVG VOL [m3],{:e},\n".format(averageVolume))
            #fw.write("AVG CELL RES T [s],{:f},\n".format(averageCellReTime))
            #fw.write("AVG CELL Q1,{:e},\n".format(avgMeshQualParameter))
            if  arguments.check:
                fw.write("AVG CELL TRNST T [s],{:f},\n".format(
                             averageTransitTime))
                fw.write("AVG CELL Q2,{:e},\n".format(avgMeshQualParameter2))
                fw.write("AVG CELL GRADxLEN,{:e},\n".format(avgGradDist))
                fw.write("AVG CELL Q3,{:e},\n".format(avgMeshQualParameter3))
            if negCheck:
                fw.write("WARNING some elements are negative")
            if arguments.source_CDGS:
                fw.write("\n")
                fw.write("\n")
                fw.write("SOURCE SAMPLING\n")
                fw.write("SAMPLING RESOLUTION [m],{:f},\n".format(
                    self.precision))
                fw.write("SAMPLED VOXELS [#],{:d},\n".format(self.xInts*
                                                             self.yInts*
                                                             self.zInts))
                fw.write("VTK EMISSION RATE [#/s],{:e},\n".format(
                    self.originalEmissionRate))
                sampString = "SAMPLED EMISSION RATE (UNSCALED) [#/s],{:e},\n"
                fw.write(sampString.format(self.unscaledEmissionRate))
            fw.write("\n")
            fw.write("\n")
            fw.write("ACTIVATION\n")
            fw.write("INLET ATOMS [#/s],{:.5e},\n".format(inletAtoms))
            fw.write("OUTLET ATOMS [#/s],{:.5e},\n".format(outletAtoms))
            fw.write("TOT IN ACTIVITY [Bq/m3],{:.5e},\n".format(
                totalInActivity))
            fw.write("TOT OUT ACTIVITY[Bq/m3],{:.5e},\n".format(
                totalOutActivity))
            fw.write("TOT CASE ACTIVITY[Bq],{:.5e},\n".format(totActivity))
            fw.write("TOT AVG ACTIVITY[Bq/m3],{:.5e},\n".format(avgActivity))
            if totalInActivity != 0:
                normAvgActivity = avgActivity/totalInActivity
                string = "INLET-NORMALIZED AVG ACTIVITY[Bq/m3],{:.5e},\n"
                fw.write(string.format(normAvgActivity))





            if inletAtoms != 0:
                reductionRate = outletAtoms/inletAtoms
                fw.write("OUT/IN RATIO,{:.5e},\n".format(reductionRate))


            fw.write("\n")
            fw.write("\n")
            fw.write("FACES\n")

            for face in self.facesPost:
                if face['typeFile'] != 'wall':
                    fw.write(face['faceID'])
                    fw.write("\n")
                    fw.write("TYPE,{},\n".format(face['typeFile']))
                    fw.write("AREA [m2],{:.5e},\n".format(face['areaFile']))
                    fw.write("FLUID FLOW [m3/s],{:.5e},\n".format(
                        face['phiFlowFileLast']))
                    fw.write("ATOM FLOW [#/s],{:.5e},\n".format(
                        face['TFlowFileLast']))
                    fw.write("ATOM CONC [#/m3],{:.5e},\n".format(
                        face['avTFileLast']))
                    fw.write("SPECIFIC ACTIVITY [Bq/m3],{:.5e},\n".format(
                            abs(face['avTFileLast'])*self.decay_constant))
                    fw.write("AVG RES T [s],{:.5f},\n".format(
                                 abs(face['avTrFileLast'])))
                    fw.write("\n")


        return

    def launch_func_object(self):
        """
        launches the utility to calculate T gradient
        """
        print ("launching gradient function object")
        orig_folder = os.getcwd()
        os.chdir(self.fluned_path)
        launch_grad = "postProcess -func 'grad(T)'".split()
        with open('log', "a",encoding='utf-8') as out_file:
            subprocess.Popen(launch_grad, stdout=out_file).wait()
        os.chdir(orig_folder)

        return

    def generate_vtk(self):
        """
        launches the utility to create a vtk file
        """

        print ("launching FoamToVTK utility")
        orig_folder = os.getcwd()
        os.chdir(self.fluned_path)
        cmd_str_1="foamToVTK -latestTime -noFaceZones -noFunctionObjects "
        cmd_str_3=" -excludePatches (\".*\")"
        launch_f2vtk = cmd_str_1.split()
        launch_f2vtk.append("-fields")
        launch_f2vtk.append("(T Ta Td)")
        launch_f2vtk.extend(cmd_str_3.split())
        #print (launch_f2vtk)
        with open('log', "a",encoding='utf-8') as outfile:
            subprocess.Popen(launch_f2vtk, stdout=outfile).wait()
        os.chdir(orig_folder)

        return

    def getVTKPath(self):
        """
        this function looks for the vtk file in the VTK folder and
        store the complete path
        """

        vtkFilePat = re.compile("\.vtk\Z", re.IGNORECASE)

        vtkFiles = []

        folder = os.path.join(self.fluned_path, 'VTK')

        self.vtk_path = ''

        for filename in os.listdir(folder):
            matchVTKfiles = vtkFilePat.findall(filename)
            if len(matchVTKfiles) == 1:
                vtkFiles.append(filename)

        if len(vtkFiles) == 1:
            self.vtk_path = os.path.join(self.fluned_path,'VTK',vtkFiles[0])
        else:
            print ("ERROR zero or more than one vtk files")
            sys.exit()



        return

    def getIsotope(self):
        """
        using the decay constant this function understand if we are
        considering N-16, N-17 or O-19 and assign the spectrum
        accordingly. If it is not possible a dummy spectrum is
        assigned
        """

        N16_decay_constant = 0.09721559
        N17_decay_constant = 0.1661825
        O19_decay_constant = 0.02578672546

        N16_branching_ratio = 0.749577982
        N17_branching_ratio = 0.951951
        O19_branching_ratio = 1.5299368

        n16Spectrum = [
        [0.0000000e+00,0          ],
        [9.8649013e-01,3.49995E-05],
        [9.8650986e-01,0          ],
        [1.7549824e+00,0.00140000],
        [1.7550175e+00,0          ],
        [1.9547805e+00,0.00040000],
        [1.9548195e+00,0          ],
        [2.7414726e+00,0.00840000],
        [2.7415274e+00,0          ],
        [2.8224718e+00,1.3000E-05],
        [2.8225282e+00,0          ],
        [6.0481395e+00,0.00013000],
        [6.0482605e+00,0          ],
        [6.1291087e+00,0.688000],
        [6.1292313e+00,0          ],
        [6.9154308e+00,0.00040000],
        [6.9155692e+00,0          ],
        [7.1150788e+00,0.05000],
        [7.1152212e+00,0          ],
        [8.8691113e+00,0.00080000],
        [8.8692887e+00,0]
        ]

        # JEFF 3.3
        n17Spectrum = [
        [0.0000000    ,0.0],
        [0.3859999    ,0.377],
        [0.3860001    ,0.0],
        [0.8859999    ,0.006],
        [0.8860001    ,0.0],
        [1.1629999    ,0.498],
        [1.1630001    ,0.0],
        [1.6889999    ,0.069],
        [1.6890001    ,0.0],
        [3.6139999    ,0.00024],
        [3.6140001    ,0.0],
        [3.8199999    ,0.00012],
        [3.8200001    ,0.0],
        ]

        o19Spectrum = [
        [0.0000000e+00,0],
        [1.0989890e-01,0.019048696],
        [1.0990110e-01,0],
        [1.9719803e-01,0.608294369],
        [1.9720197e-01,0],
        [1.3569864e+00,0.34351231],
        [1.3570136e+00,0],
        [1.4439856e+00,0.019048696],
        [1.4440144e+00,0],
        [1.5539845e+00,0.008889519],
        [1.5540155e+00,0],
        [1.5969840e+00,0.000190487],
        [1.5970160e+00,0],
        [2.5819742e+00,0.000190487],
        [2.5820258e+00,0],
        [4.1789582e+00,0.000825437],
        [4.1790418e+00,0],
        ]
        dummySpectrum = [
        [ 0 ,0 ],
        [1  , 1],
        ]
        if math.isclose(self.decay_constant,N16_decay_constant,rel_tol=1e-3):
            self.isotope = 'N16'
            spectrum =  n16Spectrum
            self.branching_ratio = N16_branching_ratio
        elif math.isclose(self.decay_constant,N17_decay_constant,
                rel_tol=1e-3):
            self.isotope = 'N17'
            spectrum =  n17Spectrum
            self.branching_ratio = N17_branching_ratio
        elif math.isclose(self.decay_constant,O19_decay_constant,
                rel_tol=1e-3):
            self.isotope = 'O19'
            spectrum =  o19Spectrum
            self.branching_ratio = O19_branching_ratio
        else:
            self.isotope = 'dummy'
            spectrum =  dummySpectrum
            print ("WARNING decay constant isotope not recognized")
            print (self.decay_constant)
            print ("a dummy spectrum is assigned")

        # normalize spectrum


        sumOriginal = sum([val[1] for val in spectrum])

        normalizingFactor = 1/sumOriginal

        for val in spectrum:
            val.append(val[1]*normalizingFactor)

        self.spectrum = spectrum



        return

    def getOriginalEmission(self):
        """
        this function reads the initial vtk and calculates the total
        emission
        """

        totalAtoms = getTotalAtomsVtk(self.vtk_path, self.dataset)

        totalEmissionVtk = (totalAtoms*
                            self.decay_constant*
                            self.branching_ratio*
                            self.scaling)

        self.originalEmissionRate = totalEmissionVtk

        return

    def calculateSamplingCoordinates(self):
        """
        self explainatory
        """


        vtkBoundaries = getVTKsize(self.vtk_path)
        xBounds = ([math.floor(vtkBoundaries[0]*100),
                    math.ceil(vtkBoundaries[1]*100)])
        yBounds = ([math.floor(vtkBoundaries[2]*100),
                    math.ceil(vtkBoundaries[3]*100)])
        zBounds = ([math.floor(vtkBoundaries[4]*100),
                    math.ceil(vtkBoundaries[5]*100)])

        samplingResolution = self.precision

        xInts = math.ceil((xBounds[1] - xBounds[0])/(samplingResolution*100))
        yInts = math.ceil((yBounds[1] - yBounds[0])/(samplingResolution*100))
        zInts = math.ceil((zBounds[1] - zBounds[0])/(samplingResolution*100))

        self.xInts = xInts
        self.yInts = yInts
        self.zInts = zInts

        xNodes = ([(xBounds[0] + samplingResolution*100*i )
            for i in range(xInts+1)])
        yNodes = ([(yBounds[0] + samplingResolution*100*i )
            for i in range(yInts+1)])
        zNodes = ([(zBounds[0] + samplingResolution*100*i )
            for i in range(zInts+1)])

        self.xNodes = xNodes
        self.yNodes = yNodes
        self.zNodes = zNodes

        voxelVector = []
        sampleCoordinates = []
        id = 1

        for i in range(xInts) :
            xVoxelNodes = [xNodes[i], xNodes[i+1]]
            xVoxelCenter = (xNodes[i+1] + xNodes[i])/2
            for j in range(yInts) :
                yVoxelNodes = [yNodes[j], yNodes[j+1]]
                yVoxelCenter = (yNodes[j+1] + yNodes[j])/2
                for k in range(zInts) :
                    zVoxelNodes = [zNodes[k], zNodes[k+1]]
                    zVoxelCenter = (zNodes[k+1] + zNodes[k])/2
                    newDic = {}
                    newDic['ID'] = id
                    newDic['centCoords'] = ([xVoxelCenter,
                                             yVoxelCenter,
                                             zVoxelCenter])
                    sampleCoordinates.append([cord/100 for cord
                        in newDic['centCoords']])
                    voxelVector.append(newDic)
                    id += 1

        self.voxelVector = voxelVector
        self.sampleCoordinates = sampleCoordinates
        self.voxelVolume = (self.precision*100)**3 # cubic cc

        return

    def sampleCoordinatesValues(self):
        """
        add the sampled value to the sample coordinates vector
        """

        vtkFile = self.vtk_path
        vtkDataSet = self.dataset
        spectrumVector = self.spectrum



        sampledRates = sampleCoordinatesVTK(vtkFile, vtkDataSet,
                                            self.sampleCoordinates)

        totalEmission = 0

        voxelVolume = self.voxelVolume

        for voxel,concentration in zip(self.voxelVector,sampledRates):
            #voxel['emittingDensity'] = concentration*1e-06 #from m3 to cm3
            voxel['emission'] = (concentration*
                                 voxelVolume*
                                 self.decay_constant*
                                 self.branching_ratio*
                                 self.scaling*
                                 1e-06)  #atoms per m3 to cm3
            #voxel['emittingSpectrum'] = [val[2]*voxel['emission']
            #                            for val in spectrumVector]

            totalEmission += voxel['emission']

        self.unscaledEmissionRate = totalEmission

        ratioVtkSampling = self.originalEmissionRate/totalEmission

        totalEmissionScaled = 0

        for voxel,concentration in zip(self.voxelVector,sampledRates):
            #voxel['emittingDensity'] = concentration*1e-06 #from m3 to cm3
            voxel['emission'] = voxel['emission']*ratioVtkSampling
            #voxel['emittingSpectrum'] = [val[2]*voxel['emission']
            #                            for val in spectrumVector]

            totalEmissionScaled += voxel['emission']

        self.scaledEmissionRate = totalEmissionScaled

        return

    def writeCDGS(self):
        """
        this function write the sample CDGS file
        """

        cdgsFile = self.vtk_path[:-3] + 'CDGS'

        spectrumVector = self.spectrum

        with open(cdgsFile,'w') as fw:
            fw.write("num_meshes 1\n")
            fw.write("global_source {:e}\n".format(self.scaledEmissionRate))
            fw.write("mesh_id 1\n")
            fw.write("cdgs from vtk\n")
            fw.write("Cooling_time 0.0\n")
            fw.write("total_source {:e}\n".format(self.scaledEmissionRate))
            fw.write("energy_type {}\n".format('bins'))
            fw.write("energy_boundaries {:d}\n".format(len(spectrumVector)))
            # WRITE SPECTRUM BINS
            specValues = [val[0] for val in spectrumVector]
            specString = formatValues(specValues)
            fw.write(specString)

            fw.write("mesh_type rec\n")
            fw.write("mesh_boundaries {:d} {:d} {:d}\n".format(self.xInts+1,
                                                               self.yInts+1,
                                                               self.zInts+1))
            fw.write("0.000000e+00  0.000000e+00  0.000000e+00\n")
            fw.write("1.000000e+00  0.000000e+00  0.000000e+00\n")
            fw.write("0.000000e+00  1.000000e+00  0.000000e+00\n")
            xString = formatValues(self.xNodes)
            fw.write(xString)
            yString = formatValues(self.yNodes)
            fw.write(yString)
            zString = formatValues(self.zNodes)
            fw.write(zString)
            fw.write("source_data\n")

            voxelString1 = "{:d} {:.5e} {:.5e} 1\n"
            voxelString2 = "0 1.0 {:.5e}\n"
            specErrorString = formatValues([0]*(len(spectrumVector)-1))

            for vox in self.voxelVector:
                if vox['emission'] > 0:
                    fw.write(voxelString1.format(vox['ID'],
                                                 vox['emission'],
                                                 self.voxelVolume))
                    fw.write(voxelString2.format(vox['emission']))
                    emittingSpectrum =  ([val[2]*vox['emission']
                                         for val in self.spectrum])
                    spectrumString=formatValues(emittingSpectrum[:-1])
                    fw.write(spectrumString)
                    fw.write(specErrorString)

            fw.write("end_source_data\n")



        return

    def get_last_time_step(self):
        """
        this function save as a class attribute the path of the last time step
        """

        folder_itms = os.listdir(self.fluned_path)
        folder_times= ([[float(itm),itm] for itm in folder_itms
                     if check_float(itm) == True])

        folder_times.sort(key=lambda x: x[0])

        self.last_time_step = os.path.join(self.fluned_path,folder_times[-1][1])

        return



def main():


    parser = argparse.ArgumentParser(description="flunedPost v0.0.1")
    parser.add_argument('-d',"--debug", action ='store_true' ,
            help="debug mode", default = False)
    parser.add_argument('-c',"--check", action ='store_true' ,
            help="check mode", default = False)
    parser.add_argument('-s',"--source_CDGS", action ='store_true' ,
            help="print cdgs", default = False)
    parser.add_argument("--precision", type=float, default=0.01,
            help="sampling precision in m")
    parser.add_argument("--dataset", type=str, default='T',
            help="dataset in the vtk file to sample")
    parser.add_argument("--scaling", type=float, default=1,
            help="scaling factor for the total emission")
    args = parser.parse_args()


    if args.debug:
        tracemalloc.start()


    simulationFolder = os.getcwd()
    simCase = flunedCase(simulationFolder)
    simCase.parseConstants()
    simCase.getIsotope()
    simCase.getTimeTreatment()
    simCase.get_last_time_step()

    if  args.check:
        simCase.launch_func_object()
        simCase.readVelocities()
        simCase.readGradT()



    simCase.readVolumes()
    simCase.readT()
    simCase.readPostProcess_flows()
    simCase.readPostProcess_Tflows()
    simCase.readPostProcess_Trflows()
    simCase.generate_vtk()





    if args.source_CDGS:
        simCase.precision = args.precision
        simCase.dataset = args.dataset
        simCase.scaling = args.scaling
        simCase.getVTKPath()
        simCase.getOriginalEmission()
        simCase.calculateSamplingCoordinates()
        simCase.sampleCoordinatesValues()
        simCase.writeCDGS()

    simCase.write_summary(args)

    print ("Finished!")

    if args.debug:
        print(tracemalloc.get_traced_memory())
        snapshot = tracemalloc.take_snapshot()
        display_top(snapshot)



    return

if __name__ == '__main__':
    main()
