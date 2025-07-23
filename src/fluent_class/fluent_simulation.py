import re
import os
import h5py
import sys
import numpy as np
from .fluned_h5_utils import get_dataset_keys
from .fluned_h5_utils import get_h5_dataset
from .fluned_h5_utils import get_h5_dataset_multi
from .fluned_h5_utils import get_h5_path_dataset
from .fluned_h5_utils import extract_multiblock


class fluentSimulation:
    """
    A class to represent an ANSYS Fluent simulation
    some functions are present to parse and generate file

    Attributes:
    -----------
    path : str
        The complete path to the simulation folder.
    """

    def __init__(self, h5_path: str, region_name):
        if os.path.isdir(h5_path):
            self.casH5file, self.datH5file = self.get_h5_files(h5_path)
        else:
            raise ValueError("folder containing the h5 file expected")

        self.fluent_path = h5_path
        self.fluid_region_name = region_name

        self.cell_density_path = ""
        self.fluid_region_id = 0
        self.min_id_region = 0
        self.data_h5_path = os.path.join(self.fluent_path, self.datH5file)
        self.case_h5_path = os.path.join(self.fluent_path, self.casH5file)

        return

    def parse_multidataset_h5_file(self):
        """
        self explanatory
        """

        self.parse_fluid_regions_data_h5()
        self.parse_fluid_faces_data_h5()
        self.fluent_density = self.read_density_h5()

        return

    def convert_to_openfoam(self, target_openfoam_path):
        self.target_openfoam_path = target_openfoam_path

        self.define_walls_multidataset_h5()
        self.write_nut_multidataset_h5()
        self.write_owner_multidataset_h5()
        self.write_neighbour_multidataset_h5()
        self.write_boundary_multidataset_h5()
        self.write_faces_multidataset_h5()
        self.write_phi_multidataset_h5()
        self.write_velocity_multidataset_h5()
        self.write_nodes_multidataset_h5()
        self.write_cell_zones()

        return

    def get_h5_files(self, folder):
        casH5files = []
        datH5files = []
        casH5FilePat = re.compile(r"\.cas.h5\Z", re.IGNORECASE)
        datH5FilePat = re.compile(r"\.dat.h5\Z", re.IGNORECASE)

        for filename in os.listdir(folder):
            casH5file = casH5FilePat.findall(filename)
            datH5file = datH5FilePat.findall(filename)
            if len(casH5file) == 1:
                casH5files.append(filename)
            if len(datH5file) == 1:
                datH5files.append(filename)

        if len(casH5files) == 1:
            fluent_case_file = casH5files[0]
        else:
            raise ValueError("ERROR zero or more than one dat.h5 files")

        if len(datH5files) == 1:
            fluent_data_file = datH5files[0]
        else:
            raise ValueError("ERROR zero or more than one dat.h5 files")

        return fluent_case_file, fluent_data_file

    def parse_fluid_regions_data_h5(self):
        """
        this function reads all the cells in the multi region file and
        stores those that are of the fluid type.
        """

        filename = self.case_h5_path

        regionNamePat = re.compile(r".*/cells/zoneTopology/name\Z", re.IGNORECASE)

        regionIDPat = re.compile(r".*/cells/zoneTopology/id\Z", re.IGNORECASE)

        # )

        regionPatMin = re.compile(r".*/cells/zoneTopology/minId\Z", re.IGNORECASE)

        regionPatMax = re.compile(r".*/cells/zoneTopology/maxId\Z", re.IGNORECASE)

        # regionMinIDPaths = []
        # regionMaxIDPaths = []

        regionNames = get_h5_dataset(filename, regionNamePat)
        regionNames = regionNames[0].decode("UTF-8").split(";")
        regionIDs = get_h5_dataset(filename, regionIDPat)
        minIDs = get_h5_dataset(filename, regionPatMin)
        maxIDs = get_h5_dataset(filename, regionPatMax)

        # define a list of dictionary with the region info
        regionList = []
        for i, val in enumerate(regionNames):
            newDict = {}
            newDict["name"] = val
            newDict["minID"] = int(minIDs[i])
            newDict["maxID"] = int(maxIDs[i])
            newDict["regionID"] = int(regionIDs[i])
            if val == self.fluid_region_name:
                newDict["fluid"] = True
                self.fluid_region_id = int(regionIDs[i])
                self.fluid_cellID_min = int(minIDs[i])
                self.fluid_cellID_max = int(maxIDs[i])
                self.fluid_cellN = self.fluid_cellID_max - self.fluid_cellID_min + 1
            else:
                newDict["fluid"] = False

            regionList.append(newDict)

        self.regionList = regionList

        return

    def read_density_h5(self):
        """
        this function reads the density stored in the h5 fluent file
        """

        filename = self.data_h5_path

        cell_density_pat = re.compile(".*/cells/SV_DENSITY/.*", re.IGNORECASE)

        cell_density_path = []

        found_density = 0

        with h5py.File(filename, "r") as fi:
            paths = get_dataset_keys(fi)

            for path in paths:
                denMatches = cell_density_pat.findall(path)
                if len(denMatches) == 1:
                    cell_density_path.append(path)

            if len(cell_density_path) == 1:
                self.cell_density_path = cell_density_path[0]
                if hasattr(self, "min_id_region"):
                    found_density = fi[self.cell_density_path][self.min_id_region - 1]
                else:
                    found_density = fi[self.cell_density_path][0]
            else:
                raise ValueError("ERROR zero or more than one density datasets found")

        return found_density

    def parse_fluid_faces_data_h5(self):
        """
        this function reads the zone topology information relative to
        the faces of the multi block model. Using this info it redefine the
        ranges of the faces to have the internal mesh at the beginning and
        the remaining faces later. Ideally these ranges will help the rest
        of the face definition.
        """

        filename = self.case_h5_path

        facesNamePat = re.compile(r".*/faces/zoneTopology/name\Z", re.IGNORECASE)
        facesPatMin = re.compile(r".*/faces/zoneTopology/minId\Z", re.IGNORECASE)

        facesPatMax = re.compile(r".*/faces/zoneTopology/maxId\Z", re.IGNORECASE)

        zoneTypePat = re.compile(r".*/faces/zoneTopology/zoneType\Z", re.IGNORECASE)

        faceOwnerPat = re.compile(r".*/faces/zoneTopology/c0\Z", re.IGNORECASE)

        faceNeighPat = re.compile(r".*/faces/zoneTopology/c1\Z", re.IGNORECASE)

        faceNames = get_h5_dataset(filename, facesNamePat)
        faceNames = list(faceNames[0].decode("UTF-8").split(";"))
        minIDs = list(get_h5_dataset(filename, facesPatMin))
        maxIDs = list(get_h5_dataset(filename, facesPatMax))
        zoneTypes = list(get_h5_dataset(filename, zoneTypePat))
        faceOwners = list(get_h5_dataset(filename, faceOwnerPat))
        faceNeighbours = list(get_h5_dataset(filename, faceNeighPat))

        sortedLists = sorted(
            zip(faceNames, minIDs, maxIDs, zoneTypes, faceOwners, faceNeighbours),
            key=lambda x: x[1],
        )

        faceNames, minIDs, maxIDs, zoneTypes, faceOwners, faceNeighbours = zip(
            *sortedLists
        )

        faceList = []

        for i, val in enumerate(faceNames):
            newDict = {}
            newDict["faceName"] = val
            newDict["origOrderID"] = i
            newDict["minID"] = int(minIDs[i])
            newDict["maxID"] = int(maxIDs[i])
            newDict["nFaces"] = newDict["maxID"] - newDict["minID"] + 1
            newDict["origNFaces"] = newDict["maxID"] - newDict["minID"] + 1
            newDict["zoneType"] = int(zoneTypes[i])
            newDict["faceOwner"] = int(faceOwners[i])
            newDict["faceNeighbour"] = int(faceNeighbours[i])
            if newDict["faceOwner"] == self.fluid_region_id:
                newDict["fluid"] = True
            else:
                newDict["fluid"] = False

            if newDict["zoneType"] == 2:
                newDict["faceType"] = "internal"
            else:
                newDict["faceType"] = "wall"
            faceList.append(newDict)

        faceList = sorted(faceList, reverse=False, key=lambda x: x["origOrderID"])

        # reorder the groups putting the internal fluidmesh at the beginning
        for face in faceList:
            if face["fluid"] and face["faceType"] == "internal":
                face["newOrderID"] = 0
                face["newMinID"] = 1
                face["newMaxID"] = face["nFaces"]
                idTemp = face["nFaces"]
                break

        orderID = 1
        for face in faceList:
            if not face["fluid"]:
                continue
            if face["fluid"]:
                if face["faceType"] == "internal":
                    continue
                else:
                    face["newOrderID"] = orderID
                    face["newMinID"] = idTemp + 1
                    face["newMaxID"] = idTemp + face["nFaces"]
                    idTemp = face["newMaxID"]
                    orderID += 1

        for face in faceList:
            if face["fluid"]:
                continue
            else:
                face["newOrderID"] = orderID
                face["newMinID"] = idTemp + 1
                face["newMaxID"] = idTemp + face["nFaces"]
                idTemp = face["newMaxID"]
                orderID += 1

        # faceList = sorted(faceList, reverse=False,
        #                              key=lambda x:x['newOrderID'])

        # for dic in faceList:
        #    print (dic)

        nWaterOwner = 0

        for face in faceList:
            if face["fluid"]:
                nWaterOwner += face["nFaces"]

        self.faceList = faceList
        self.nWaterFaces = nWaterOwner

        return

    def define_walls_multidataset_h5(self):
        """
        this function just check the phi values to distinguish between
        wall, inlet and outlet
        """

        filename = self.data_h5_path

        phiPat = re.compile(".*/faces/SV_FLUX/.*", re.IGNORECASE)

        phiValues = get_h5_dataset_multi(filename, phiPat)

        for face in self.faceList:
            if not face["fluid"]:
                continue
            if (face["fluid"]) and face["faceType"] == "internal":
                continue

            phiTemp = extract_multiblock(face["minID"], face["maxID"], phiValues)

            if all(val == 0 for val in phiTemp):
                face["faceType"] = "wall"
            elif all(val >= 0 for val in phiTemp):
                face["faceType"] = "outlet"
                if not all(val > 0 for val in phiTemp):
                    print("WARNING in outlet ", face["faceName"])
                    print("some zero phis are present")
            elif all(val <= 0 for val in phiTemp):
                face["faceType"] = "inlet"
                if not all(val < 0 for val in phiTemp):
                    print("WARNING in inlet ", face["faceName"])
                    print("some zero phis are present")
            else:
                print("ERROR in face ", face["faceName"])
                print("nonzero phis of different signs are present")
                sys.exit()

        # for dic in self.faceList:
        #    print (dic)
        # df = pd.DataFrame(self.faceList)
        # df.to_csv('faces.csv', index=True)

        return

    def write_owner_multidataset_h5(self):
        """
        This function writes all the cell owners,
        it selects those that belong to the fluid partition
        """

        print("writing multiblock openFOAM owner ... ")

        filename = self.case_h5_path

        ownerPat = re.compile(r".*/faces/c0/\d+\Z", re.IGNORECASE)

        ownerFilePath = os.path.join(
            self.target_openfoam_path, "constant", "polyMesh", "owner"
        )

        ownerHeader = """
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

        self.ownerH5Path, ownerList = get_h5_path_dataset(filename, ownerPat)

        ownerTable = np.zeros((self.nWaterFaces,), dtype=np.uint64)

        for face in self.faceList:
            if face["fluid"]:
                ownerTable[(face["newMinID"] - 1) : face["newMaxID"]] = ownerList[
                    (face["minID"] - 1) : face["maxID"]
                ]

        # check = np.count_nonzero(ownerTable == 0)
        # print (check)

        with open(ownerFilePath, "w") as fo:
            fo.write(ownerHeader)
            fo.write(str(self.nWaterFaces) + "\n")
            fo.write("(\n")
            owStr = "{:d}\n"
            for val in ownerTable:
                fo.write(owStr.format(int(val - self.fluid_cellID_min)))
            fo.write(")")

        return

    def write_neighbour_multidataset_h5(self):
        """
        this function writes all the cell neighbours,
        it selects those that belongs to the fluid partition
        """

        print("writing multiblock openFOAM neighbour ... ")

        filename = self.case_h5_path

        neighPat = re.compile(r".*/faces/c1/\d+\Z", re.IGNORECASE)

        neighbourFilePath = os.path.join(
            self.target_openfoam_path, "constant", "polyMesh", "neighbour"
        )

        neighHeader = """
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
        neighData = get_h5_dataset_multi(filename, neighPat)

        for face in self.faceList:
            if (face["fluid"]) and face["faceType"] == "internal":
                nWaterNeigh = face["nFaces"]
                neighTable = extract_multiblock(face["minID"], face["maxID"], neighData)
                break

        with open(neighbourFilePath, "w") as fo:
            fo.write(neighHeader)
            fo.write(str(nWaterNeigh) + "\n")
            fo.write("(\n")
            neStr = "{:d}\n"
            for val in neighTable:
                fo.write(neStr.format(int(val - self.fluid_cellID_min)))
            fo.write(")")

        return

    def write_boundary_multidataset_h5(self):
        """
        this function writes the boundary file for a simulation
        extracted from a multi block fluent simulation
        the self.faceList should be already arranged for the work
        """

        boundaryFilePath = os.path.join(
            self.target_openfoam_path, "constant", "polyMesh", "boundary"
        )
        boundaryHeader = """
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
            if face["fluid"] and (face["faceType"] != "internal"):
                nWaterFaces += 1

        faceList = sorted(self.faceList, reverse=False, key=lambda x: x["newOrderID"])

        with open(boundaryFilePath, "w") as fi:
            fi.write(boundaryHeader)
            fi.write("{:d}\n".format(nWaterFaces))
            fi.write("(\n")
            for face in faceList:
                if not face["fluid"]:
                    continue
                if face["faceType"] == "internal":
                    continue

                fi.write("   {}\n".format(face["faceName"]))
                fi.write("    {\n")

                if face["faceType"] == "wall":
                    fi.write("        type            wall;\n")
                else:
                    fi.write("        type            patch;\n")

                fi.write("        nFaces          {:d};\n".format(face["nFaces"]))
                fi.write(
                    "        startFace       {:d};\n".format(int(face["newMinID"]) - 1)
                )
                fi.write("    }\n")
            fi.write(")\n")

        return

    def write_faces_multidataset_h5(self):
        """
        this function write the face definition of the mesh extracted from a
        multiblock fluent simulation
        """

        print("writing multiblock openFOAM faces ... ")

        filename = self.case_h5_path

        facePat = re.compile(r".*/faces/nodes/\d+/nnodes\Z", re.IGNORECASE)
        face2Pat = re.compile(r".*/faces/nodes/\d+/nodes\Z", re.IGNORECASE)

        facesFilePath = os.path.join(
            self.target_openfoam_path, "constant", "polyMesh", "faces"
        )

        facesHeader = """
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

        self.faces1H5Path, nFacesVec = get_h5_path_dataset(filename, facePat)
        self.faces2H5Path, facesDef = get_h5_path_dataset(filename, face2Pat)

        pointIDTemp = 1
        uniquePoints = np.zeros((0,), dtype=np.uint64)
        for face in self.faceList:
            # store the number of nodes per face
            nPointsTemp = nFacesVec[(face["minID"] - 1) : face["maxID"]]
            # number of nodes needed to define one face
            face["nPoints"] = sum(nPointsTemp)
            face["minPointID"] = pointIDTemp
            face["maxPointID"] = face["minPointID"] + face["nPoints"] - 1
            pointIDTemp += face["nPoints"]
            # print (face['faceName'], face['minPointID'],face['maxPointID'])
            if face["fluid"]:
                pointsTemp = facesDef[(face["minPointID"] - 1) : face["maxPointID"]]
                # print ("fluid face", min(pointsTemp),max(pointsTemp))
                uniquePoints = np.append(uniquePoints, pointsTemp)

        self.uniquePoints = np.unique(uniquePoints)

        # print ("CHECK")
        # print (len(facesDef))
        # print (pointIDTemp)

        faceList = sorted(self.faceList, reverse=False, key=lambda x: x["newOrderID"])

        with open(facesFilePath, "w") as fo:
            fo.write(facesHeader)
            fo.write(str(self.nWaterFaces) + "\n")
            fo.write("(\n")
            j = 0
            for face in faceList:
                if not face["fluid"]:
                    continue
                nPointsTemp = nFacesVec[(face["minID"] - 1) : face["maxID"]]
                pointsTemp = facesDef[(face["minPointID"] - 1) : face["maxPointID"]]
                i = 0
                # print (len(nPointsTemp))
                # print (sum(nPointsTemp))
                # print (len(pointsTemp))
                for val in nPointsTemp:
                    faceStr1 = "{:d}("
                    newString = faceStr1.format(val)
                    pointList = reversed(pointsTemp[i : (i + val)])
                    for pt in pointList:
                        ptNew = np.searchsorted(self.uniquePoints, pt)
                        newString += " {:d} ".format(int(ptNew))
                    newString += " )\n"
                    fo.write(newString)
                    i += val
                    j += val

            fo.write(")")

        faceZonesFilePath = os.path.join(
            self.target_openfoam_path, "constant", "polyMesh", "faceZones"
        )

        faceZones = """
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
        with open(faceZonesFilePath, "w") as fo2:
            fo2.write(faceZones)

        return

    def write_nodes_multidataset_h5(self):
        """
        this function writes the node of a converted multiblock fluent
        simulation. This function limits the writing to the points present
        in the fluid region written in the self.uniquePoints attribute
        """

        print("writing multiblock openFOAM points ... ")

        filename = self.case_h5_path

        nodesPat = re.compile(r".*/nodes/coords/\d+\Z", re.IGNORECASE)

        pointsFilePath = os.path.join(
            self.target_openfoam_path, "constant", "polyMesh", "points"
        )

        pointsHeader = """
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
        self.nodesH5Path, points = get_h5_path_dataset(filename, nodesPat)

        with open(pointsFilePath, "w") as fo:
            fo.write(pointsHeader)
            fo.write(str(len(self.uniquePoints)) + "\n")
            fo.write("(\n")
            ptStr = "({:24.18e} {:24.18e} {:24.18e})\n"
            for index in self.uniquePoints:
                cord = points[int(index - 1)]

                fo.write(ptStr.format(*cord))

            fo.write(")")

        pointsZonesFilePath = os.path.join(
            self.target_openfoam_path, "constant", "polyMesh", "pointZones"
        )

        pointZones = """
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
        with open(pointsZonesFilePath, "w") as fo2:
            fo2.write(pointZones)

        return

    def write_phi_multidataset_h5(self):
        """
        this function works for a multi-block fluent simulation and does
        three things:
        - write the phi values at the end of the simulation
        - scan the phi values to understand which are the wall and which the
          inlet/outlet. It updates the fluentFacesVector vector with the
          results
        """

        print("writing multiblock openFOAM phi files ... ")

        filename = self.data_h5_path

        phiPat = re.compile(".*/faces/SV_FLUX/.*", re.IGNORECASE)

        phiOneFilePath = os.path.join(self.target_openfoam_path, "1", "phi")

        phiValues = get_h5_dataset_multi(filename, phiPat)

        phi1Header = """
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

        faceList = sorted(self.faceList, reverse=False, key=lambda x: x["newOrderID"])

        with open(phiOneFilePath, "w") as fo:
            fo.write(phi1Header)
            fo.write("internalField   nonuniform List<scalar>\n")
            # write the internal mesh
            for face in faceList:
                if (face["fluid"]) and face["faceType"] == "internal":
                    phiTemp = extract_multiblock(
                        face["minID"], face["maxID"], phiValues
                    )
                    fo.write("{:d}\n".format(len(phiTemp)))
                    fo.write("(\n")
                    phiStr = "{:24.18e}\n"
                    for val in phiTemp:
                        fo.write(phiStr.format(val / self.fluent_density))
            fo.write(")\n")
            fo.write(";\n")
            fo.write("\n")
            fo.write("boundaryField\n")
            fo.write("{\n")

            for face in faceList:
                if face["faceType"] == "internal":
                    continue
                if face["fluid"]:
                    phiTemp = extract_multiblock(
                        face["minID"], face["maxID"], phiValues
                    )
                    fo.write("    {}\n".format(face["faceName"]))
                    fo.write("    {\n")
                    fo.write("        type            calculated;\n")
                    if face["faceType"] == "wall":
                        fo.write("        value           uniform 0;\n")
                    if face["faceType"] != "wall":
                        fo.write("        value   nonuniform List<scalar>\n")
                        phiTemp = extract_multiblock(
                            face["minID"], face["maxID"], phiValues
                        )
                        fo.write("{:d}\n".format(len(phiTemp)))
                        fo.write("(\n")
                        for val in phiTemp:
                            fo.write(phiStr.format(val / self.fluent_density))
                        fo.write(")\n")
                        fo.write(";\n")
                    fo.write("    }\n")

            fo.write("}\n")

        return

    def write_velocity_multidataset_h5(self):
        """
        this function write the U velocity files for t=1 and t=0
        """

        print("writing multiblock openFOAM U files ... ")

        filename = self.data_h5_path

        uCellXPat = re.compile(".*/cells/SV_U/.*", re.IGNORECASE)
        uCellYPat = re.compile(".*/cells/SV_V/.*", re.IGNORECASE)
        uCellZPat = re.compile(".*/cells/SV_W/.*", re.IGNORECASE)

        uFaceXPat = re.compile(".*/faces/SV_U/.*", re.IGNORECASE)
        uFaceYPat = re.compile(".*/faces/SV_V/.*", re.IGNORECASE)

        uFaceZPat = re.compile(".*/faces/SV_W/.*", re.IGNORECASE)

        uOneFilePath = os.path.join(self.target_openfoam_path, "1", "U")

        xCellVelocity = get_h5_dataset_multi(filename, uCellXPat)
        yCellVelocity = get_h5_dataset_multi(filename, uCellYPat)
        zCellVelocity = get_h5_dataset_multi(filename, uCellZPat)

        xFaceVelocity = get_h5_dataset_multi(filename, uFaceXPat)
        yFaceVelocity = get_h5_dataset_multi(filename, uFaceYPat)
        zFaceVelocity = get_h5_dataset_multi(filename, uFaceZPat)

        u1Header = """
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
        faceList = sorted(self.faceList, reverse=False, key=lambda x: x["newOrderID"])

        with open(uOneFilePath, "w") as fo:
            uStr = "({:24.18e} {:24.18e} {:24.18e})\n"
            fo.write(u1Header)
            fo.write("internalField   nonuniform List<vector>\n")
            fo.write("{:d}\n".format(self.fluid_cellN))
            fo.write("(\n")
            uXValues = extract_multiblock(
                self.fluid_cellID_min, self.fluid_cellID_max, xCellVelocity
            )
            uYValues = extract_multiblock(
                self.fluid_cellID_min, self.fluid_cellID_max, yCellVelocity
            )
            uZValues = extract_multiblock(
                self.fluid_cellID_min, self.fluid_cellID_max, zCellVelocity
            )
            for uX, uY, uZ in zip(uXValues, uYValues, uZValues):
                fo.write(uStr.format(uX, uY, uZ))
            fo.write(")\n")
            fo.write(";\n")
            fo.write("\n")
            fo.write("boundaryField\n")
            fo.write("{\n")

            for face in faceList:
                if face["faceType"] == "internal":
                    continue
                if face["fluid"]:
                    fo.write("    {}\n".format(face["faceName"]))
                    fo.write("    {\n")
                    if face["faceType"] == "wall":
                        fo.write("        type        movingWallVelocity;\n")
                        fo.write("        value       uniform (0 0 0);\n")
                    elif face["faceType"] == "outlet":
                        fo.write("        type            zeroGradient;\n")
                    elif face["faceType"] == "inlet":
                        uXValues = extract_multiblock(
                            face["minID"], face["maxID"], xFaceVelocity
                        )
                        uYValues = extract_multiblock(
                            face["minID"], face["maxID"], yFaceVelocity
                        )
                        uZValues = extract_multiblock(
                            face["minID"], face["maxID"], zFaceVelocity
                        )
                        fo.write("        type            fixedValue;\n")
                        fo.write("        value   nonuniform List<vector>\n")
                        fo.write("{:d}\n".format(len(uXValues)))
                        fo.write("(\n")
                        for v1, v2, v3 in zip(uXValues, uYValues, uZValues):
                            fo.write(uStr.format(v1, v2, v3))
                        fo.write(")\n")
                        fo.write(";\n")
                    fo.write("    }\n")

            fo.write("}\n")

        return

    def write_nut_multidataset_h5(self):
        """
        this file write the turbulent viscosity in the 1 folder, it is
        assumed this values is always calculated: meaning that we always
        convert from a turbulent simulation.
        """
        print("writing multiblock openFOAM nut files ... ")

        filename = self.data_h5_path
        filename_case = self.case_h5_path

        mutCellPat = re.compile(".*/cells/SV_MU_T/.*", re.IGNORECASE)

        # mutCellValues = get_h5_dataset(filename,mutCellPat)
        mutCellValues = get_h5_dataset_multi(filename, mutCellPat)

        nutOneFilePath = os.path.join(self.target_openfoam_path, "1", "nut")

        ownerPat = re.compile(r".*/faces/c0/\d+\Z", re.IGNORECASE)

        ownerList = get_h5_dataset(filename_case, ownerPat)

        nut1Header = """
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

        faceList = sorted(self.faceList, reverse=False, key=lambda x: x["newOrderID"])

        mutValues = extract_multiblock(
            self.fluid_cellID_min, self.fluid_cellID_max, mutCellValues
        )

        with open(nutOneFilePath, "w") as fo:
            nutStr = "{:24.18e}\n"
            fo.write(nut1Header)
            fo.write("internalField   nonuniform List<scalar>\n")
            fo.write("{:d}\n".format(self.fluid_cellN))
            fo.write("(\n")
            for val in mutValues:
                fo.write(nutStr.format(val / self.fluent_density))

            fo.write(")\n")
            fo.write(";\n")

            fo.write("\n")
            fo.write("boundaryField\n")
            fo.write("{\n")

            for face in faceList:
                if face["faceType"] == "internal":
                    continue
                if face["fluid"]:
                    fo.write("    {}\n".format(face["faceName"]))
                    fo.write("    {\n")
                    if face["faceType"] == "wall":
                        fo.write("        type         nutkWallFunction;\n")
                        fo.write("        blending            stepwise;\n")
                        fo.write("        Cmu            0.09;\n")
                        fo.write("        kappa            0.41;\n")
                        fo.write("        E            9.8;\n")
                        fo.write("        value     uniform  0;\n")
                    else:
                        fo.write("        type            calculated;\n")
                        fo.write("        value   nonuniform List<scalar>\n")
                        faceOwn = ownerList[(face["minID"] - 1) : face["maxID"]]
                        fo.write("{:d}\n".format(len(faceOwn)))
                        fo.write("(\n")
                        for ow in faceOwn:
                            fo.write(
                                nutStr.format(
                                    mutValues[int(ow - self.fluid_cellID_min)]
                                    / self.fluent_density
                                )
                            )
                        fo.write(")\n")
                        fo.write(";\n")
                    fo.write("    }\n")

            fo.write("}\n")

        return

    def write_cell_zones(self):
        """
        function to write the cell zones file in the polyMesh folder
        """
        cell_zones_file_path = os.path.join(
            self.target_openfoam_path, "constant", "polyMesh", "cellZones"
        )

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
        with open(cell_zones_file_path, "w", encoding="utf-8") as fo2:
            fo2.write(cell_zones_string)

        return
