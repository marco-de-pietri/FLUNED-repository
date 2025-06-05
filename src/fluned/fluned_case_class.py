import os
import re
import sys
import pathlib
import shutil
import math
import numpy as np
import subprocess
import h5py
import pyvista as pv
from vtk.util import numpy_support as VN
import vtk
from .util import check_int
from .util import mod
from .util import check_float
from .util import formatValues
from .util import sample_coordinates_vtk
from .util import open_utf8_or_gzip
from isotopes.isotopes import load_isotopes
from fluent_parser.fluned_h5_utils import get_dataset_keys
from fluent_parser.fluned_h5_utils import get_h5_dataset
from fluent_parser.fluned_h5_utils import get_h5_dataset_multi
from fluent_parser.fluned_h5_utils import get_h5_path_dataset
from fluent_parser.fluned_h5_utils import extract_multiblock
from fluent_parser.fluned_bin_utils import get_fluent_binarray_double
from fluent_parser.fluned_bin_utils import get_fluent_parse_headers
from fluent_parser.fluned_bin_utils import get_fluent_parse_regions
from ofClass.ofClass import SimulationOF
from ofClass.ofClass import merge_continue_runs


class flunedCase:
    """
    FLUNED simulation class
    """

    def __init__(self, fluned_path=""):
        """
        initialize flunedCase
        """
        self.decay_constant = 0
        self.activation_normalization = 0
        self.activation_dataset = ""
        self.activation_dataset_error = ""
        self.activation_const = 0
        self.activation_file = ""
        self.fv_scheme = "stable"
        self.fluned_path = fluned_path
        self.case = os.path.split(fluned_path)[1]
        self.num_internal_cells = 0

    def generate_case(self, arg_dict):
        """
        fill simulation parameters
        """

        self.case = arg_dict["case"]
        # set decay constant
        if "decay_constant" in arg_dict:
            self.decay_constant = float(arg_dict["decay_constant"])

        if "isotope" in arg_dict:
            isotope = arg_dict["isotope"].lower().replace("-", "")
            if isotope not in ["n16", "o19", "n17", "f20", "custom"]:
                raise ValueError("isotope not recognized")
            self.isotope = isotope
        else:
            self.isotope = "custom"

        if "activation_const" in arg_dict:
            self.activation_const = float(arg_dict["activation_const"])

        if "activation_normalization" in arg_dict:
            self.activation_normalization = float(arg_dict["activation_normalization"])

        if "activation_dataset" in arg_dict:
            self.activation_dataset = arg_dict["activation_dataset"]

        if "activation_dataset_error" in arg_dict:
            self.activation_dataset_error = arg_dict["activation_dataset_error"]

        if "activation_file" in arg_dict:
            self.activation_file = os.path.normcase(arg_dict["activation_file"])

        if "fv_scheme" in arg_dict:
            if arg_dict["fv_scheme"].lower() in ["stable", "accurate"]:
                self.fv_scheme = arg_dict["fv_scheme"]
            else:
                raise ValueError("fv_scheme must be 'stable' or 'accurate'")

        if "cfd_type" not in arg_dict:
            print("ERROR: cfd type not specified")
            sys.exit()
        else:
            self.cfd_type = arg_dict["cfd_type"].lower()

        if "time_treatment" not in arg_dict:
            raise ValueError("ERROR: type of time treatment not specified")
        elif arg_dict["time_treatment"].lower() not in ["steadystate", "transient"]:
            err_string = "ERROR: time treatment must be 'steadyState' or 'transient'"
            raise ValueError(err_string)
        else:
            self.time_treatment = arg_dict["time_treatment"].lower()

        if self.cfd_type in ["fluent-h5-multi", "fluent-multi"]:
            if "fluent_fluid_region_name" not in arg_dict:
                print("ERROR: name of the fluid region to extract not")
                print("specified! use parameter FLUENT_FLUID_REGION_NAME")
                sys.exit()
            else:
                self.fluent_fluid_region_name = arg_dict["fluent_fluid_region_name"]

        self.molecular_diffusion = float(arg_dict["molecular_diffusion"])
        self.schmidt_number = float(arg_dict["schmidt_number"])
        self.inlet_conc = float(arg_dict["inlet_conc"])

        self.cfd_path = os.path.normcase(arg_dict["cfd_path"])

        if not os.path.isdir(self.cfd_path):
            raise OSError(f"Folder not found: {self.cfd_path}")

        self.fluned_path = os.path.join(self.cfd_path, self.case)

        if not os.path.exists(self.fluned_path):
            os.mkdir(self.fluned_path)

        return

    def initialize_cfd_class(self):
        """
        this function initialize the cfd class, after the initial openfoam or
        fluent simulation has been processed
        """

        self.cfd_simulation = SimulationOF(self.cfd_path)

        return

    def generate_zero_ta(self):
        """
        this function generate the Ta file at t=0
        """

        zero_folder = os.path.join(self.fluned_path, "0")
        zero_ta_path = os.path.join(zero_folder, "Ta")

        ta_header_text = """
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

        with open(zero_ta_path, "w", encoding="utf-8") as fw:
            fw.write(ta_header_text)

            for face in self.faces:
                fw.write("    " + face["faceID"] + "\n    {\n")
                if face["type"] == "wall":
                    fw.write("        type            zeroGradient;\n")
                elif face["type"] == "outlet":
                    fw.write("        type            zeroGradient;\n")
                elif face["type"] == "inlet":
                    fw.write("        type            fixedValue;\n")
                    vals_string = "        value          uniform 0;\n"
                    fw.write(vals_string)
                else:
                    print("ERROR face type not recognized")
                    print(face["type"])
                    sys.exit()

                fw.write("    }\n")

            fw.write(end_text)

        return

    def generate_zero_tr(self):
        """
        this function generate the Tr file at t=0
        """

        zero_folder = os.path.join(self.fluned_path, "0")
        zero_tr_path = os.path.join(zero_folder, "Tr")

        tr_header_text = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Tr;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

"""
        end_text = """

}
"""

        with open(zero_tr_path, "w", encoding="utf-8") as fw:
            fw.write(tr_header_text)

            for face in self.faces:
                fw.write("    " + face["faceID"] + "\n    {\n")
                if face["type"] == "wall":
                    fw.write("        type            zeroGradient;\n")
                elif face["type"] == "outlet":
                    fw.write("        type            zeroGradient;\n")
                elif face["type"] == "inlet":
                    fw.write("        type            fixedValue;\n")
                    if self.time_treatment == "steadystate":
                        vals_string = "        value          uniform 0;\n"
                    elif self.time_treatment == "transient":
                        vals_string = "        value          uniform 1;\n"
                    fw.write(vals_string)
                else:
                    print("ERROR face type not recognized")
                    print(face["type"])
                    sys.exit()

                fw.write("    }\n")

            fw.write(end_text)

        return

    def generate_zero_td(self):
        """
        this function generate the Td file at t=0
        """

        zero_folder = os.path.join(self.fluned_path, "0")
        zero_td_path = os.path.join(zero_folder, "Td")

        td_header_text = """
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

        with open(zero_td_path, "w", encoding="utf-8") as fw:
            fw.write(td_header_text)

            for face in self.faces:
                fw.write("    " + face["faceID"] + "\n    {\n")
                if face["type"] == "wall":
                    fw.write("        type            zeroGradient;\n")
                elif face["type"] == "outlet":
                    fw.write("        type            zeroGradient;\n")
                elif face["type"] == "inlet":
                    fw.write("        type            fixedValue;\n")
                    vals_string = "        value         uniform {};\n"
                    fw.write(vals_string.format(self.inlet_conc))
                else:
                    print("ERROR face type not recognized")
                    print(face["type"])
                    sys.exit()

                fw.write("    }\n")

            fw.write(end_text)

        return

    def generate_zero_t(self):
        """
        this function generate the T file at t=0
        """

        zero_folder = os.path.join(self.fluned_path, "0")
        zero_t_path = os.path.join(zero_folder, "T")

        t_header_text = """
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

        with open(zero_t_path, "w", encoding="utf-8") as fw:
            fw.write(t_header_text)

            for face in self.faces:
                fw.write("    " + face["faceID"] + "\n    {\n")
                if face["type"] == "wall":
                    fw.write("        type            zeroGradient;\n")
                elif face["type"] == "outlet":
                    fw.write("        type            zeroGradient;\n")
                elif face["type"] == "inlet":
                    fw.write("        type            fixedValue;\n")
                    vals_string = "        value         uniform {};\n"
                    fw.write(vals_string.format(self.inlet_conc))
                else:
                    print("ERROR face type not recognized")
                    print(face["type"])
                    sys.exit()

                fw.write("    }\n")

            fw.write(end_text)

        return

    def generate_constant_file(self):
        """
        This function creates the files in the constant folder
        """

        transport_prop_text = """
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

isotope          {};

DT             DT [0 2 -1 0 0 0 0] {};

lambda         lambda [0 0 -1 0 0 0 0] {};

Sct            Sct [ 0 0 0 0 0 0 0 ] {};


// ************************************************************************ //
"""

        constant_folder = os.path.join(self.fluned_path, "constant")
        transport_prop_path = os.path.join(constant_folder, "transportProperties")

        with open(transport_prop_path, "w", encoding="utf-8") as fw:
            fw.write(
                transport_prop_text.format(
                    self.isotope,
                    self.molecular_diffusion,
                    self.decay_constant,
                    self.schmidt_number,
                )
            )

        return

    def create_case_folders_fluent(self):
        """
        This function creates the folders for the FLUENT simulation converted
        into the OpenFOAM format
        """

        case_folder = self.cfd_path

        # T=0 folder
        zero_folder = os.path.join(case_folder, "0")
        if not os.path.exists(zero_folder):
            os.mkdir(zero_folder)

        # T=1 folder
        one_folder = os.path.join(case_folder, "1")
        if not os.path.exists(one_folder):
            os.mkdir(one_folder)

        # constant folder
        const_folder = os.path.join(case_folder, "constant")
        if not os.path.exists(const_folder):
            os.mkdir(const_folder)

        # polyMesh folder
        poly_folder = os.path.join(case_folder, "constant", "polyMesh")
        if not os.path.exists(poly_folder):
            os.mkdir(poly_folder)

        # system folder
        sys_folder = os.path.join(case_folder, "system")
        if not os.path.exists(sys_folder):
            os.mkdir(sys_folder)

        # case foam file
        case_file = os.path.join(case_folder, "case.foam")
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
        zero_folder = os.path.join(case_folder, "0")
        if not os.path.exists(zero_folder):
            os.mkdir(zero_folder)

        # constant folder
        const_folder = os.path.join(case_folder, "constant")
        if not os.path.exists(const_folder):
            os.mkdir(const_folder)

        # polyMesh folder
        poly_folder = os.path.join(case_folder, "constant", "polyMesh")
        if not os.path.exists(poly_folder):
            os.mkdir(poly_folder)

        # system folder
        sys_folder = os.path.join(case_folder, "system")
        if not os.path.exists(sys_folder):
            os.mkdir(sys_folder)

        # case foam file
        case_file = os.path.join(case_folder, "case.foam")
        pathlib.Path(case_file).touch()

        return

    def copy_last_phi(self):
        """
        this function look for the last phi file in the cfd folder
        if in the cfd simulation the flow is in m3/s it just copies it to the
        fluned folders. If the flows is in kg/s it converts it to m3/s
        """

        last_time = self.cfd_simulation.last_time
        target_file = os.path.join(self.fluned_path, "0", "phi")
        origin_file = os.path.join(self.cfd_path, str(last_time), "phi")

        if self.cfd_simulation.volumetric_flag:
            shutil.copyfile(origin_file, target_file)
        else:
            print("converting phi values to volumetric flow ...")
            self.cfd_simulation.convert_phi_to_volumetric(target_file)

        return

    def copy_last_u(self):
        """
        this function look for the last U file in the cfd folder
        """

        last_time = self.cfd_simulation.last_time
        target_file = os.path.join(self.fluned_path, "0", "U")
        origin_file = os.path.join(self.cfd_path, str(last_time), "U")
        shutil.copyfile(origin_file, target_file)

        return None

    def copyLastNut(self):
        """this function look for the last nut file in the cfd folder,
        if it does not exists it means the simulation is laminar"""

        folderItems = os.listdir(self.cfd_path)

        folderTimes = [int(itm) for itm in folderItems if check_int(itm)]

        lastTime = max(folderTimes)

        targetFile = os.path.join(self.fluned_path, "0", "nut")
        originFile = os.path.join(self.cfd_path, str(lastTime), "nut")
        checkFile = os.path.isfile(originFile)
        if checkFile:
            shutil.copyfile(originFile, targetFile)

        return

    def copyPolyMesh(self):
        """this function copy the poly Mesh files in the FLUNED folder"""

        sourceFolder = os.path.join(self.cfd_path, "constant", "polyMesh")
        targetFolder = os.path.join(self.fluned_path, "constant", "polyMesh")

        # shutil do not allow to copy into an existing folder
        dirCheck = os.path.isdir(targetFolder)
        if dirCheck:
            shutil.rmtree(targetFolder)

        shutil.copytree(sourceFolder, targetFolder)

        return

    def reconstructFaces(self):
        """this function examines the polyMesh data and reconstruct which
        faces are input/output/wall"""

        polyMeshFolder = os.path.join(self.fluned_path, "constant", "polyMesh")
        boundaryFile = os.path.join(polyMeshFolder, "boundary")

        try:
            inpFile = open(boundaryFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open boundary file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            faceDefPat = re.compile(
                r"\d+[\n\r\s]+?\(.*?[\n\r\s]+?\)", re.MULTILINE | re.DOTALL
            )
            boundaryPat = re.compile(
                r"[^\s]+[\n\r\s]+?\{.*?\}", re.MULTILINE | re.DOTALL
            )
            boundaryNamePat = re.compile(
                r"([^\s]+)[\n\r\s]+?\{.*?\}", re.MULTILINE | re.DOTALL
            )
            faceNPat = re.compile("nFaces.*?(\d+)")
            firstFacePat = re.compile("startFace.*?(\d+)")
            defBlock = faceDefPat.findall(text)[0]
            # print (defBlock)
            # print (faceNumber)
            boundaryBlocks = boundaryPat.findall(defBlock)
            # print (boundaryBlocks)

            # here it get the description of the boundary surfaces
            boundaryVec = []
            for block in boundaryBlocks:
                boundaryDic = {}
                boundaryDic["faceID"] = boundaryNamePat.findall(block)[0]
                # print (boundaryDic['faceID'])
                boundaryDic["nFaces"] = int(faceNPat.findall(block)[0])
                boundaryDic["firstFace"] = int(firstFacePat.findall(block)[0])
                boundaryDic["faces"] = list(
                    range(
                        boundaryDic["firstFace"],
                        boundaryDic["firstFace"] + boundaryDic["nFaces"],
                    )
                )
                boundaryVec.append(boundaryDic)
            # print (boundaryVec)

        phifilename = os.path.join(self.fluned_path, "0", "phi")
        try:
            inpFile = open(phifilename, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open phi file")
            sys.exit()
        with inpFile:
            facePhiPat = re.compile(r"\((.{1,}?)\)", re.MULTILINE | re.DOTALL)
            text = inpFile.read()
            wallFacePat = re.compile("value\s+uniform\s+0")

            # print ("face phi")
            for face in boundaryVec:
                faceBlockPat2 = re.compile(
                    face["faceID"] + r"[\n\r\s]+?\{.*?\}", re.MULTILINE | re.DOTALL
                )
                faceBloc = faceBlockPat2.findall(text)[0]
                # print (faceBloc)
                facePhis = facePhiPat.findall(faceBloc)
                # print (facePhis)
                wallConfirm = wallFacePat.findall(faceBloc)
                # print (wallConfirm)
                if len(facePhis) == 0 and len(wallConfirm) != 0:
                    face["type"] = "wall"
                    face["phis"] = np.zeros(face["nFaces"])
                else:
                    phiList = facePhiPat.findall(faceBloc)[0].strip().split("\n")
                    phiList = np.array([float(val) for val in phiList])

                    if all(i >= 0 for i in phiList):
                        face["type"] = "outlet"
                    elif all(i <= 0 for i in phiList):
                        face["type"] = "inlet"
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

includeDecayScalar true;




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

includeDecayScalar true;




functions
{

"""
        vol_calc_text = """
    volumeCalc
    {{
        type            writeCellVolumes;
        libs            ("libfieldFunctionObjects.so");
	    select      all;

        writeFields     false;
        writeControl {};

    }}

 """

        vol_flow_text = """
    {}
    {{

        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        patch   {};
        fields  (phi);

        operation sum;
        select  patch;
        name        $patch;

        writeFields     false;
        writeControl {};

    }}

 """
        vol_tx_flow_text = """
    {}
    {{

        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        patch   {};
        fields  ({});
        weightField phi;

        operation sum;
        select  patch;
        name        $patch;

        writeFields     false;
        writeControl {};

    }}

 """

        vol_tx_sum_text = """
    {}
    {{

        type            volFieldValue;
        libs            ("libfieldFunctionObjects.so");
        fields  ({});

        operation       volAverage;
	    select      all;

        writeFields     false;
        writeControl {};


    }}

 """
        system_folder = os.path.join(self.fluned_path, "system")
        control_dict_path = os.path.join(system_folder, "controlDict")

        with open(control_dict_path, "w", encoding="utf-8") as fw:
            if self.time_treatment == "steadystate":
                fw.write(control_dict_text)
                write_control = "outputTime"
            elif self.time_treatment == "transient":
                fw.write(control_dict_text_transient)
                write_control = "timeStep"

            fw.write(vol_calc_text.format(write_control))

            fw.write(vol_tx_sum_text.format("volTSum", "T", write_control))
            fw.write(vol_tx_sum_text.format("volTaSum", "Ta", write_control))
            fw.write(vol_tx_sum_text.format("volTdSum", "Td", write_control))

            for face in self.faces:
                if face["type"] in ["inlet", "outlet"]:
                    fw.write(
                        vol_flow_text.format(
                            "volFlow-" + face["faceID"], face["faceID"], write_control
                        )
                    )

                    fw.write(
                        vol_tx_flow_text.format(
                            "volTFlow-" + face["faceID"],
                            face["faceID"],
                            "T",
                            write_control,
                        )
                    )
                    fw.write(
                        vol_tx_flow_text.format(
                            "volTrFlow-" + face["faceID"],
                            face["faceID"],
                            "Tr",
                            write_control,
                        )
                    )
                    fw.write(
                        vol_tx_flow_text.format(
                            "volTdFlow-" + face["faceID"],
                            face["faceID"],
                            "Td",
                            write_control,
                        )
                    )
                    fw.write(
                        vol_tx_flow_text.format(
                            "volTaFlow-" + face["faceID"],
                            face["faceID"],
                            "Ta",
                            write_control,
                        )
                    )

            fw.write("}")

        fv_scheme_stable = """
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
    default         cellLimited Gauss linear 1;
}}

divSchemes
{{
    default         none;
    div(phi,T)      Gauss  upwind;
    div(phi,Ta)      Gauss upwind;
    div(phi,Td)      Gauss upwind;
    div(phi,Tr)      Gauss upwind;
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
    default         limited 0.5;
}}

"""

        fv_scheme_accurate = """
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

        schemes_path = os.path.join(system_folder, "fvSchemes")

        with open(schemes_path, "w", encoding="utf-8") as fw:
            # select the fv scheme
            if self.fv_scheme == "accurate":
                fv_scheme_text = fv_scheme_accurate
            elif self.fv_scheme == "stable":
                fv_scheme_text = fv_scheme_stable

            # select the time treatment
            elif self.time_treatment == "transient":
                fw.write(fv_scheme_text.format("Euler"))
            if self.time_treatment == "steadystate":
                fw.write(fv_scheme_text.format("steadyState"))
            elif self.time_treatment == "transient":
                fw.write(fv_scheme_text.format("Euler"))

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
        solution_path = os.path.join(system_folder, "fvSolution")

        with open(solution_path, "w", encoding="utf-8") as fw:
            if self.time_treatment == "steadystate":
                fw.write(fv_solution_text.format(0.01, 0.01, 0.01, 0.01))
            elif self.time_treatment == "transient":
                fw.write(fv_solution_text.format(0, 0, 0, 0))

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
        parallel_dict_path = os.path.join(system_folder, "decomposeParDict")
        with open(parallel_dict_path, "w", encoding="utf-8") as fw:
            fw.write(parallel_dict_text)

        return

    def launch_vol_func_object(self):
        """
        this function launch the utilities that calculates the volumes
        - this field is needed only for the activation rate
        interpolation
        """

        print("calculating volumes  ...")
        orig_folder = os.getcwd()
        os.chdir(self.fluned_path)
        launch_volumes = "postProcess -func writeCellVolumes".split()
        with open("log", "a", encoding="utf-8") as outfile:
            subprocess.Popen(launch_volumes, stdout=outfile).wait()
        os.chdir(orig_folder)

        return None

    def launchCentroidFuncObjects(self):
        """this function launch the utilities that calculates the volume and
        the centroids - this fields are needed only for the activation rate
        interpolation"""

        print("calculating centroids ...")
        origFolder = os.getcwd()
        os.chdir(self.fluned_path)
        launchCentroids = "postProcess -func writeCellCentres".split()
        with open("log", "a") as outfile:
            _ = subprocess.Popen(launchCentroids, stdout=outfile).wait()
        os.chdir(origFolder)

        return

    def launch_solver(self):
        """
        this function launches the scalar calculation
        """

        print("launching FLUNED solver ...")
        current_folder = os.getcwd()
        os.chdir(self.fluned_path)
        launch_calc_string = "FLUNED-solver".split()
        with open("simulation_log", "a", encoding="utf-8") as outfile:
            subprocess.Popen(launch_calc_string, stdout=outfile).wait()
        os.chdir(current_folder)

        return

    def read_volumes(self):
        """
        this function reads the volumes from the V file located in the
        zero folder
        """

        # common patterns
        internalBlockPat = re.compile(
            r"internalField.*?\((.{1,}?)\)", re.MULTILINE | re.DOTALL
        )

        v_file_path = os.path.join(self.fluned_path, "0")
        v_files = [f for f in os.listdir(v_file_path) if re.match(r"V(c)?(\..*)?$", f)]
        if not v_files:
            print("No V or Vc files found")
            sys.exit()

        v_file = os.path.join(v_file_path, v_files[0])

        text = open_utf8_or_gzip(v_file)
        numInternalBlocks = internalBlockPat.findall(text)
        internalVolumes = numInternalBlocks[0].split("\n")[1:-1]
        self.Volumes = np.array([float(val) for val in internalVolumes])
        self.num_internal_cells = len(internalVolumes)

        return

    # def read_volumes(self):
    #     """
    #     this function reads the volumes from the V file located in the
    #     zero folder
    #
    #     fluned post
    #     """
    #
    #     print("reading volume values...")
    #
    #     # common patterns
    #     internalBlockPat = re.compile(
    #         r"internalField.*?\((.{1,}?)\)", re.MULTILINE | re.DOTALL
    #     )
    #
    #     nElPat = re.compile("internalField.*?(\d+).*?\(", re.MULTILINE | re.DOTALL)
    #     v_file_path = os.path.join(self.fluned_path, "0")
    #     v_files = [f for f in os.listdir(v_file_path) if re.match(r"V(c)?(\..*)?$", f)]
    #     if not v_files:
    #         print("No V or Vc files found")
    #         sys.exit()
    #
    #     v_file = os.path.join(v_file_path, v_files[0])
    #
    #     text = open_utf8_or_gzip(v_file)
    #
    #     cellNumberText = nElPat.findall(text)
    #     self.num_internal_cells = int(cellNumberText[0])
    #     numInternalBlocks = internalBlockPat.findall(text)
    #     internalVolumes = numInternalBlocks[0].split("\n")[1:-1]
    #     self.Volumes = np.zeros(self.num_internal_cells)
    #     for i in range(self.num_internal_cells):
    #         self.Volumes[i] = float(internalVolumes[i])
    #
    #     return

    def readCentroids(self):
        """this function reads the centroids from the 0 folder"""

        # common patterns
        internalBlockPat = re.compile(
            "internalField.*?\((.{1,}?)\n\\s*\)", re.MULTILINE | re.DOTALL
        )

        cFile = os.path.join(self.fluned_path, "0", "C")
        try:
            inpFile = open(cFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open  C file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            numInternalBlocks = internalBlockPat.findall(text)
            internalCentroids = numInternalBlocks[0].split("\n")[1:]
            internalCentroids = [val.strip("()") for val in internalCentroids]
            self.Centroids = np.array(
                [[float(val) for val in v.split()] for v in internalCentroids]
            )

        return

    def write_cell_zones(self):
        """
        function to write the cell zones file in the polyMesh folder
        """
        cell_zones_file_path = os.path.join(
            self.cfd_path, "constant", "polyMesh", "cellZones"
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

    def writeNeighbour_multi_h5(self):
        """
        this function writes all the cell neighbours,
        it selects those that belongs to the fluid partition
        """

        print("writing multiblock openFOAM neighbour ... ")

        filename = os.path.join(self.cfd_path, self.casH5file)

        neighPat = re.compile(r".*/faces/c1/\d+\Z", re.IGNORECASE)

        neighbourFilePath = os.path.join(
            self.cfd_path, "constant", "polyMesh", "neighbour"
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

    def writeFaces_multi_h5(self):
        """
        this function write the face definition of the mesh extracted from a
        multiblock fluent simulation
        """

        print("writing multiblock openFOAM faces ... ")

        filename = os.path.join(self.cfd_path, self.casH5file)

        facePat = re.compile(r".*/faces/nodes/\d+/nnodes\Z", re.IGNORECASE)
        face2Pat = re.compile(r".*/faces/nodes/\d+/nodes\Z", re.IGNORECASE)

        facesFilePath = os.path.join(self.cfd_path, "constant", "polyMesh", "faces")

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
            self.cfd_path, "constant", "polyMesh", "faceZones"
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

    def writeBoundary_multi_h5(self):
        """
        this function writes the boundary file for a simulation
        extracted from a multi block fluent simulation
        the self.faceList should be already arranged for the work
        """

        boundaryFilePath = os.path.join(
            self.cfd_path, "constant", "polyMesh", "boundary"
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

        return

    def writeOwner_multi_h5(self):
        """This function writes all the cell owners,
        it selects those that belong to the fluid partition"""

        print("writing multiblock openFOAM owner ... ")

        filename = os.path.join(self.cfd_path, self.casH5file)

        ownerPat = re.compile(r".*/faces/c0/\d+\Z", re.IGNORECASE)

        ownerFilePath = os.path.join(self.cfd_path, "constant", "polyMesh", "owner")

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

    def writeNodes_multi(self):
        """
        this function writes the node of a converted multiblock fluent
        simulation. This function limits the writing to the points present
        in the fluid region written in the self.uniquePoints attribute
        """

        print("writing multiblock openFOAM points ... ")

        filename = os.path.join(self.cfd_path, self.casfile)

        get_fluent_binarray_double(filename, "3010", 3)

        pointsZonesFilePath = os.path.join(
            self.cfd_path, "constant", "polyMesh", "pointZones"
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

    def writeNodes_multi_h5(self):
        """
        this function writes the node of a converted multiblock fluent
        simulation. This function limits the writing to the points present
        in the fluid region written in the self.uniquePoints attribute
        """

        print("writing multiblock openFOAM points ... ")

        filename = os.path.join(self.cfd_path, self.casH5file)

        nodesPat = re.compile(r".*/nodes/coords/\d+\Z", re.IGNORECASE)

        pointsFilePath = os.path.join(self.cfd_path, "constant", "polyMesh", "points")

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
            self.cfd_path, "constant", "polyMesh", "pointZones"
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

    def getH5files(self):
        casH5files = []
        datH5files = []
        casH5FilePat = re.compile(r"\.cas.h5\Z", re.IGNORECASE)
        datH5FilePat = re.compile(r"\.dat.h5\Z", re.IGNORECASE)

        folder = self.cfd_path

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
            print("ERROR zero or more than one cas.h5 files")
            sys.exit()

        if len(datH5files) == 1:
            self.datH5file = datH5files[0]
        else:
            print("ERROR zero or more than one dat.h5 files")
            sys.exit()

        return

    def getCASDATfiles(self):
        casfiles = []
        datfiles = []
        casFilePat = re.compile(r"\.cas\Z", re.IGNORECASE)
        datFilePat = re.compile(r"\.dat\Z", re.IGNORECASE)

        folder = self.cfd_path

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
            print("ERROR zero or more than one cas file")
            sys.exit()

        if len(datfiles) == 1:
            self.datfile = datfiles[0]
        else:
            print("ERROR zero or more than one dat file")
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

        zeroFolder = os.path.join(self.fluned_path, "0")
        zeroSourcePath = os.path.join(zeroFolder, "TrSource")

        sHeaderText = """
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

        with open(zeroSourcePath, "w") as fw:
            fw.write(sHeaderText)

            if self.time_treatment == "steadystate":
                fw.write(intFieldText.format(1))
            elif self.time_treatment == "transient":
                fw.write(intFieldText.format(0))

            fw.write(boundaryText)

            for face in self.faces:
                fw.write("    " + face["faceID"] + "\n    {\n")
                fw.write("        type            fixedValue;\n")
                valString = "        value          uniform 0;\n"
                fw.write(valString)
                fw.write("    }\n")

            fw.write(closerText)

        return

    def generateSourceFile(self):
        """this file create the source file in the zero folder.
        There are three modes:
        1) no activation, then the file contains only zeros.
        2) constant activation (with input value)
        3) with a source file
        """

        if self.activation_file == "":
            if self.activation_const == 0:
                # case with no rad source
                activ_sources = [0 for i in range(self.num_internal_cells)]
            else:
                # case with constant value source
                self.read_volumes()
                activ_sources = [self.activation_const * vol for vol in self.Volumes]

        else:
            # case with source file
            self.launchCentroidFuncObjects()
            self.read_volumes()
            self.readCentroids()

            # 1.sample the activation file
            print("Sampling Reaction Rate file ... ")
            sampledRates = sample_coordinates_vtk(
                self.activation_file, self.activation_dataset, self.Centroids
            )

            # 1.1 if present sample the vtk file to get the error array
            if self.activation_dataset_error != "":
                print("Sampling Reaction Rate MCNP errors  ... ")
                sampledStatErr = sample_coordinates_vtk(
                    self.activation_file, self.activation_dataset_error, self.Centroids
                )

            # 2.use the activation const as a factor
            if self.activation_const == 0:
                factor = 1
            else:
                factor = self.activation_const
            activ_sources = [factor * rate for rate in sampledRates]

            # 3.apply a normalization factor if provided
            if self.activation_normalization != 0:
                vec = [rate * vol for rate, vol in zip(activ_sources, self.Volumes)]

                totalSampled = sum(vec)
                print("total sampled atoms/s")
                print(totalSampled)
                normFactor = self.activation_normalization / totalSampled

                activ_sources = [rate * normFactor for rate in activ_sources]

                nVec = [rate * vol for rate, vol in zip(activ_sources, self.Volumes)]

                print("new total sampled atoms/s")
                print(sum(nVec))

        self.activation_sources = activ_sources

        zeroFolder = os.path.join(self.fluned_path, "0")
        zeroSourcePath = os.path.join(zeroFolder, "Source")

        sHeaderText = """
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

        with open(zeroSourcePath, "w") as fw:
            fw.write(sHeaderText)
            fw.write("{:d}\n".format(self.num_internal_cells))
            fw.write("(\n")
            for val in self.activation_sources:
                fw.write("{:e}\n".format(val))
            fw.write(")\n;\n\n")

            fw.write(boundaryText)

            for face in self.faces:
                fw.write("    " + face["faceID"] + "\n    {\n")
                fw.write("        type            fixedValue;\n")
                valString = "        value          uniform 0;\n"
                fw.write(valString)
                fw.write("    }\n")

            fw.write(closerText)

        if self.activation_dataset_error != "":
            zeroSourceErrorPath = os.path.join(zeroFolder, "SourceError")

            eHeaderText = """
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
            with open(zeroSourceErrorPath, "w") as fw:
                fw.write(eHeaderText)
                fw.write("{:d}\n".format(self.num_internal_cells))
                fw.write("(\n")
                for val in sampledStatErr:
                    fw.write("{:e}\n".format(val))
                fw.write(")\n;\n\n")

                fw.write(boundaryText)

                for face in self.faces:
                    fw.write("    " + face["faceID"] + "\n    {\n")
                    fw.write("        type            fixedValue;\n")
                    valString = "        value          uniform 0;\n"
                    fw.write(valString)
                    fw.write("    }\n")

                fw.write(closerText)

        return

    def writeSpeed_multi_h5(self):
        """
        this function write the U velocity files for t=1 and t=0
        """

        print("writing multiblock openFOAM U files ... ")

        filename = os.path.join(self.cfd_path, self.datH5file)

        uCellXPat = re.compile(".*/cells/SV_U/.*", re.IGNORECASE)
        uCellYPat = re.compile(".*/cells/SV_V/.*", re.IGNORECASE)
        uCellZPat = re.compile(".*/cells/SV_W/.*", re.IGNORECASE)

        uFaceXPat = re.compile(".*/faces/SV_U/.*", re.IGNORECASE)
        uFaceYPat = re.compile(".*/faces/SV_V/.*", re.IGNORECASE)

        uFaceZPat = re.compile(".*/faces/SV_W/.*", re.IGNORECASE)

        uOneFilePath = os.path.join(self.cfd_path, "1", "U")

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

    def defineWalls_multi_h5(self):
        """
        this function just check the phi values to distinguish between
        wall, inlet and outlet
        """

        filename = os.path.join(self.cfd_path, self.datH5file)

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

    def writePhi_multi_h5(self):
        """
        this function works for a multi-block fluent simulation and does
        three things:
        - write the phi values at the end of the simulation
        - scan the phi values to understand which are the wall and which the
          inlet/outlet. It updates the fluentFacesVector vector with the
          results
        """

        print("writing multiblock openFOAM phi files ... ")

        filename = os.path.join(self.cfd_path, self.datH5file)

        phiPat = re.compile(".*/faces/SV_FLUX/.*", re.IGNORECASE)

        phiOneFilePath = os.path.join(self.cfd_path, "1", "phi")

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

    def getFluidFaces(self):
        """
        this function reads the zone topology information relative to
        the faces of the multi block model. Using this info it redefine the
        ranges of the faces to have the internal mesh at the beginning and
        the remaining faces later. Ideally these ranges will help the rest
        of the face definition.
        """

        # filename = os.path.join(self.cfd_path, self.casfile)

        regionPat1 = re.compile(r"\([^()]*?\)")

        nameMatches = regionPat1.findall(self.cfdPostMeshInfoString)

        faceNames = [name.strip("()") for name in nameMatches[1:]]

        faceList = []

        for i, val in enumerate(faceNames):
            newDict = {}
            name = val.split()[0]
            zoneTypeString = val.split()[1]
            ownerName = val.split()[2]
            if len(val.split()) > 3:
                neighbourName = val.split()[3]
            else:
                neighbourName = ""

            newDict["faceName"] = name
            newDict["origOrderID"] = i

            # look for the face ID
            filtVect = [v for v in self.headerDictCas if v["index"] == 39]
            for entry in filtVect:
                if entry["header"].split()[2] == name:
                    newDict["regionID"] = int(entry["header"].split()[0])

            # look for the min and max ID
            filtVect = [v for v in self.headerDictCas if v["index"] == 3013]
            for entry in filtVect:
                headerVec = []
                for val in entry["header"].split():
                    intVal = int(val, 16)
                    headerVec.append(intVal)

                if headerVec[0] == newDict["regionID"]:
                    newDict["minID"] = int(headerVec[1])
                    newDict["maxID"] = int(headerVec[2])
                    newDict["nFaces"] = newDict["maxID"] - newDict["minID"] + 1
                    newDict["zoneType"] = int(headerVec[3])
                    newDict["zoneTypeString"] = zoneTypeString
                    for region in self.regionList:
                        if region["name"] == ownerName:
                            newDict["faceOwner"] = int(region["regionID"])
                    if neighbourName == "":
                        newDict["faceNeighbour"] = 0
                    else:
                        for region in self.regionList:
                            if region["name"] == neighbourName:
                                newDict["faceNeighbour"] = int(region["regionID"])

                    if newDict["faceOwner"] == self.fluid_region_id:
                        newDict["fluid"] = True
                    else:
                        newDict["fluid"] = False

                    # up to this point i read only boundary walls
                    newDict["faceType"] = "wall"

            faceList.append(newDict)

        storedFaces = [v["faceName"] for v in faceList]
        storedRegions = [v["name"] for v in self.regionList]
        stored = storedFaces + storedRegions
        # maxI = max([v["origOrderID"] for v in faceList])
        i += 1

        print(storedFaces)

        filtVect = [v for v in self.headerDictCas if v["index"] == 39]

        for entry in filtVect:
            name = entry["header"].split()[2]

            if name not in stored:
                ownerName = name.split("-")[-1]
                print(ownerName)

                newDict = {}
                newDict["faceName"] = name
                newDict["origOrderID"] = i
                i += 1
                newDict["regionID"] = int(entry["header"].split()[0])

                # look for the min and max ID
                filtVct1 = [v for v in self.headerDictCas if v["index"] == 3013]
                for entry2 in filtVct1:
                    headerVec = []
                    for val in entry2["header"].split():
                        intVal = int(val, 16)
                        headerVec.append(intVal)
                    if headerVec[0] == newDict["regionID"]:
                        newDict["faceType"] = "internal"
                        newDict["minID"] = int(headerVec[1])
                        newDict["maxID"] = int(headerVec[2])
                        newDict["nFaces"] = newDict["maxID"] - newDict["minID"] + 1
                        newDict["zoneType"] = int(headerVec[3])
                        newDict["zoneTypeString"] = "internal"

                        for region in self.regionList:
                            if region["name"] == ownerName:
                                newDict["faceOwner"] = int(region["regionID"])
                                newDict["faceNeighbour"] = newDict["faceOwner"]

                        if newDict["faceOwner"] == self.fluid_region_id:
                            newDict["fluid"] = True
                        else:
                            newDict["fluid"] = False

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

        faceList = sorted(faceList, reverse=False, key=lambda x: x["newOrderID"])

        nWaterOwner = 0

        for face in faceList:
            if face["fluid"]:
                nWaterOwner += face["nFaces"]

        self.faceList = faceList
        self.nWaterFaces = nWaterOwner

        for dic in faceList:
            print(dic)

        return

    def getFluidFaces_h5(self):
        """
        this function reads the zone topology information relative to
        the faces of the multi block model. Using this info it redefine the
        ranges of the faces to have the internal mesh at the beginning and
        the remaining faces later. Ideally these ranges will help the rest
        of the face definition.
        """

        filename = os.path.join(self.cfd_path, self.casH5file)

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

    def createHeaderDictionary(self):
        """
        this function parse all the headers in the cas/dat file and generate
        a dictionary that can be used in the other functions to get the data
        """

        filename = os.path.join(self.cfd_path, self.casfile)

        self.headerDictCas = get_fluent_parse_headers(filename)

        return

    def getFluidCells(self):
        """
        this function reads all the cells in the multi region file and
        stores those that are of the fluid type.
        """

        # get region names

        # regionPat = re.compile(r"\(.*?\)")
        regionPat1 = re.compile(r"\([^()]*?\)")
        filename = os.path.join(self.cfd_path, self.casfile)
        regionNames = get_fluent_parse_regions(filename)
        self.cfdPostMeshInfoString = regionNames

        regionNames = regionPat1.findall(regionNames)

        regionNames = regionNames[0].strip("()").split()

        ## define a list of dictionary with the region info
        regionList = []
        for name in regionNames:
            newDict = {}
            newDict["name"] = name

            # look for the region ID
            filtVect = [v for v in self.headerDictCas if v["index"] == 39]
            for entry in filtVect:
                if entry["header"].split()[2] == name:
                    newDict["regionID"] = int(entry["header"].split()[0])

            # look for the min and max ID
            filtVect = [v for v in self.headerDictCas if v["index"] == 3012]
            for entry in filtVect:
                headerVec = []
                for val in entry["header"].split():
                    intVal = int(val, 16)
                    headerVec.append(intVal)

                if headerVec[0] == newDict["regionID"]:
                    newDict["minID"] = int(headerVec[1])
                    newDict["maxID"] = int(headerVec[2])

            if newDict["name"] == self.fluent_fluid_region_name:
                newDict["fluid"] = True
                self.fluid_region_id = int(newDict["regionID"])
                self.fluid_cellID_min = int(newDict["minID"])
                self.fluid_cellID_max = int(newDict["maxID"])
            else:
                newDict["fluid"] = False

            regionList.append(newDict)

        self.regionList = regionList

        # for v in self.regionList:
        #    print(v)

        return

    def getFluidCells_h5(self):
        """
        this function reads all the cells in the multi region file and
        stores those that are of the fluid type.
        """

        filename = os.path.join(self.cfd_path, self.casH5file)

        regionNamePat = re.compile(r".*/cells/zoneTopology/name\Z", re.IGNORECASE)

        regionIDPat = re.compile(r".*/cells/zoneTopology/id\Z", re.IGNORECASE)

        # partitionIdPat = re.compile(
        #     r".*/cells/partition/.*/partition-ids\Z", re.IGNORECASE
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
            if val == self.fluent_fluid_region_name:
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

    def writeNut_multi_h5(self):
        """
        this file write the turbulent viscosity in the 1 folder, it is
        assumed this values is always calculated: meaning that we always
        convert from a turbulent simulation.
        """
        print("writing multiblock openFOAM nut files ... ")

        filename = os.path.join(self.cfd_path, self.datH5file)

        mutCellPat = re.compile(".*/cells/SV_MU_T/.*", re.IGNORECASE)

        # mutCellValues = get_h5_dataset(filename,mutCellPat)
        mutCellValues = get_h5_dataset_multi(filename, mutCellPat)

        nutOneFilePath = os.path.join(self.cfd_path, "1", "nut")

        filenameCAS = os.path.join(self.cfd_path, self.casH5file)

        ownerPat = re.compile(r".*/faces/c0/\d+\Z", re.IGNORECASE)

        ownerList = get_h5_dataset(filenameCAS, ownerPat)

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

    def getNumCells(self):
        """this function create an attribute to the class that specifies the
        number of internal cells. It does so by reading the U file"""

        folderItems = os.listdir(self.cfd_path)

        folderTimes = [int(itm) for itm in folderItems if check_int(itm)]

        lastTime = max(folderTimes)

        velocityFile = os.path.join(self.cfd_path, str(lastTime), "U")

        internalBlockPat = re.compile(
            "internalField.*?(\d+).*?\(", re.MULTILINE | re.DOTALL
        )

        try:
            inpFile = open(velocityFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open velocity file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            cellNumber = internalBlockPat.findall(text)[0]

        self.num_internal_cells = int(cellNumber)

        return

    def read_density_h5(self):
        """
        this function reads the density stored in the h5 fluent file
        """

        filename = os.path.join(self.cfd_path, self.datH5file)

        cell_density_pat = re.compile(".*/cells/SV_DENSITY/.*", re.IGNORECASE)

        cell_density_path = []

        with h5py.File(filename, "r") as fi:
            paths = get_dataset_keys(fi)

            for path in paths:
                denMatches = cell_density_pat.findall(path)
                if len(denMatches) == 1:
                    cell_density_path.append(path)

            if len(cell_density_path) == 1:
                self.cell_density_path = cell_density_path[0]
                if hasattr(self, "minIDregion"):
                    self.fluent_density = fi[self.cell_density_path][
                        self.minIDregion - 1
                    ]
                else:
                    self.fluent_density = fi[self.cell_density_path][0]
            else:
                raise ValueError("ERROR zero or more than one density datasets found")

        return

    def parse_constant(self):
        """this function parses the constant properties to get the decay
        variable and the others"""

        # common patterns
        dtPat = re.compile("DT\s*DT.*")
        isotope_pat = re.compile("isotope.*")
        lambdaPat = re.compile("lambda\s*lambda.*")
        schPat = re.compile(r"Sct\s*Sct.*")

        cFile = os.path.join(self.fluned_path, "constant", "transportProperties")
        try:
            inpFile = open(cFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open transportProperties file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            dtLines = dtPat.findall(text)
            lambdaLines = lambdaPat.findall(text)
            schLines = schPat.findall(text)
            isotopeLines = isotope_pat.findall(text)
            if len(isotopeLines) != 0:
                vals = isotopeLines[0].strip(" ;").split()
                val = vals[-1]
                self.isotope = val
            else:
                self.isotope = "custom"

            if len(dtLines) != 0:
                vals = dtLines[0].strip(" ;").split()
                val = vals[-1]
                self.molecular_dif = float(val)
            else:
                self.molecular_dif = 0

            if len(lambdaLines) != 0:
                vals = lambdaLines[0].strip(" ;").split()
                val = vals[-1]
                self.decay_constant = float(val)
            else:
                self.decay_constant = 0

            if len(schLines) != 0:
                vals = schLines[0].strip(" ;").split()
                val = vals[-1]
                self.schmidt_number = float(val)
            else:
                self.schmidt_number = 0

        return

    def get_isotope(self):
        """
        using the decay constant this function understand if we are
        considering N-16, N-17 or O-19 and assign the spectrum
        accordingly. If it is not possible a dummy spectrum is
        assigned

        this will be changed by embedding the data in the simulation
        """

        N16_decay_constant = 0.09721559
        N17_decay_constant = 0.1661825
        O19_decay_constant = 0.02578672546
        F20_decay_constant = (0.062093271,)

        dummy_branching_ratio = 1

        if math.isclose(self.decay_constant, N16_decay_constant, rel_tol=1e-3):
            self.isotope = "n16"
        elif math.isclose(self.decay_constant, N17_decay_constant, rel_tol=1e-3):
            self.isotope = "n17"
        elif math.isclose(self.decay_constant, O19_decay_constant, rel_tol=1e-3):
            self.isotope = "o19"
        elif math.isclose(self.decay_constant, F20_decay_constant, rel_tol=1e-3):
            self.isotope = "f20"
        else:
            print("ERROR, isotope not recognized")
            sys.exit()

        if self.isotope != "dummy":
            isotope_database = load_isotopes()
            if self.isotope.lower() not in isotope_database:
                print("ERROR isotope not found in the database")
                sys.exit()
            isotope_data = isotope_database[self.isotope.lower()]
            self.e_bins = isotope_data.e_bins
            self.p_bins = isotope_data.p_bins
            self.e_lines = isotope_data.e_lines
            self.p_lines = isotope_data.p_lines
            self.branching_ratio = isotope_data.branching_ratio
            self.particle_type = isotope_data.emitting_particle
        else:
            self.branching_ratio = dummy_branching_ratio

        return

    def get_time_treatment(self):
        """
        this function parse the fvSchemes file to check if the simulation
        is steady state or transient
        """

        # common patterns
        ddtPat = re.compile("ddtSchemes.*?\{.{1,}?\}", re.MULTILINE | re.DOTALL)

        cFile = os.path.join(self.fluned_path, "system", "fvSchemes")
        try:
            inpFile = open(cFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open fvSchemes file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            ddtText = ddtPat.findall(text)[0]
            if "euler" in ddtText.lower():
                self.time_treatment = "transient"
            elif "steadystate" in ddtText.lower():
                self.time_treatment = "steadystate"
            else:
                print("couldn't open recognise the time scheme")
                sys.exit()

        return

    def get_last_time_step(self):
        """
        this function save as a class attribute the path of the last time step
        """

        folder_itms = os.listdir(self.fluned_path)
        folder_times = [[float(itm), itm] for itm in folder_itms if check_float(itm)]

        folder_times.sort(key=lambda x: x[0])

        self.last_time_step = os.path.join(self.fluned_path, folder_times[-1][1])

        return

    def launch_grad_func_object(self):
        """
        launches the utility to calculate T gradient
        """
        print("launching gradient function object")
        orig_folder = os.getcwd()
        os.chdir(self.fluned_path)
        launch_grad = "postProcess -func 'grad(T)'".split()
        with open("log", "a", encoding="utf-8") as out_file:
            subprocess.Popen(launch_grad, stdout=out_file).wait()
        os.chdir(orig_folder)

        return

    def generate_vtk(self):
        """
        launches the utility to create a vtk file
        """

        print("launching FoamToVTK utility")
        orig_folder = os.getcwd()
        os.chdir(self.fluned_path)
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

    def read_velocities(self):
        """this function reads the velocity from the U file located in the
        zero folder"""

        print("reading velocity values...")

        # common patterns
        internalBlockPat = re.compile(
            r"internalField.*?\((.{1,}?)\n\\s*\)", re.MULTILINE | re.DOTALL
        )

        uFile = os.path.join(self.fluned_path, "0", "U")
        try:
            inpFile = open(uFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open  U file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            numInternalBlocks = internalBlockPat.findall(text)
            internalVels = numInternalBlocks[0].split("\n")[1:]
            internalVels = [val.strip("()") for val in internalVels]

            self.Velocities = np.array(
                [np.array([float(val) for val in v.split()]) for v in internalVels]
            )

        return

    def read_t(self):
        """
        this function reads the T scalar
        """

        print("reading scalar values...")

        # find the last folder
        # folderItms = os.listdir(self.fluned_path)
        # folderTimes=[int(itm) for itm in folderItms if checkInt(itm) == True]
        # lastTime = max(folderTimes)

        # common patterns
        internalBlockPat = re.compile(
            r"internalField.*?\((.{1,}?)\)", re.MULTILINE | re.DOTALL
        )

        tFile = os.path.join(self.last_time_step, "T")

        try:
            inpFile = open(tFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open Volume T file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            numInternalBlocks = internalBlockPat.findall(text)
            internalScalar = numInternalBlocks[0].split("\n")[1:-1]

            print("test", self.num_internal_cells)
            print("test2", len(internalScalar))

            self.TScalar = np.zeros(self.num_internal_cells)
            for i in range(self.num_internal_cells):
                self.TScalar[i] = float(internalScalar[i])

        return

    def create_results_folder(self):
        """
        this function creates the results folder
        """

        resFolder = os.path.join(self.fluned_path, "RESULTS")

        dir_check = os.path.isdir(resFolder)
        if not dir_check:
            os.mkdir(resFolder)

        self.results_folder = resFolder

        return

    def get_vtk_path(self, custom_vtk=False):
        """
        this function looks for the vtk file in the VTK folder generated by the
        vtk openfoam utility and store it for generation of results
        custom values can be provided
        """

        vtkFiles = []

        folder = os.path.join(self.fluned_path, "VTK")

        self.vtk_path = ""

        if custom_vtk:
            filename = os.path.splitext(custom_vtk)[0]
        else:
            filename = ""

        pat_string = r"{}\.vtk\Z".format(filename)

        vtkFilePat = re.compile(pat_string, re.IGNORECASE)

        for filename in os.listdir(folder):
            matchVTKfiles = vtkFilePat.findall(filename)
            if len(matchVTKfiles) == 1:
                vtkFiles.append(filename)

        if len(vtkFiles) == 1:
            self.vtk_path = os.path.join(self.fluned_path, "VTK", vtkFiles[0])
        else:
            print("ERROR zero or more than one vtk files")
            print("found in the VTK folder")
            print("found vtk files: ", vtkFiles)
            sys.exit()

        return

    def write_external_stl(self):
        """
        this function extracts the external surface of the final vtk file and
        write it in a stl file
        """

        stl_path = os.path.join(self.results_folder, "external_surface.stl")

        # read the VTK mesh (handles .vtk, .vtu, .vti, etc.)
        mesh = pv.read(self.vtk_path)

        # extract the exterior faces (drops internal cells)
        surface = mesh.extract_surface().triangulate()

        # write STL (binary by default)
        surface.save(stl_path)

        return

    def read_grad_t(self):
        """this function reads the T gradient"""

        print("reading gradient values...")

        # common patterns
        internalBlockPat = re.compile(
            r"internalField.*?\((.{1,}?)\n\\s*\)", re.MULTILINE | re.DOTALL
        )

        gradFile = os.path.join(self.last_time_step, "grad(T)")
        try:
            inpFile = open(gradFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open grad file")
            sys.exit()
        with inpFile:
            text = inpFile.read()
            numInternalBlocks = internalBlockPat.findall(text)
            internalVels = numInternalBlocks[0].split("\n")[1:]
            internalVels = [val.strip("()") for val in internalVels]
            self.Gradients = np.array(
                [np.array([float(val) for val in v.split()]) for v in internalVels]
            )

        return

    def readPostProcess_flows(self):
        """
        this function reads the post process folders and extract from all
        the folders the flows, area and type of the faces
        """

        print("reconstructing interface properties...")

        postFolder = os.path.join(self.fluned_path, "postProcessing")

        # read the files in the postprocess folder
        postFlowVec = []

        folderItms = os.listdir(postFolder)
        flowFolders = [itm for itm in folderItms if itm[0:8] == "volFlow-"]

        for itm in flowFolders:
            postDic = {}
            postDic["faceID"] = itm[8:]
            postDic["folder"] = itm
            postDic["flowFiles"] = self.getPostFiles(os.path.join(postFolder, itm))
            postDic["areaFile"] = self.postFileArea(postDic["flowFiles"][0])

            timeLists = self.postFileArray(postDic["flowFiles"], "Time")

            flowLists = self.postFileArray(postDic["flowFiles"], "sum(phi)")

            timeListSorted, flowListSorted = merge_continue_runs(
                timeLists,
                flowLists,
            )
            postDic["timeFile"] = timeListSorted
            postDic["phiFlowFile"] = flowListSorted
            postDic["phiFlowFileLast"] = postDic["phiFlowFile"][-1]

            if postDic["phiFlowFileLast"] > 0:
                postDic["typeFile"] = "outlet"
            elif postDic["phiFlowFileLast"] == 0:
                postDic["typeFile"] = "wall"
            elif postDic["phiFlowFileLast"] < 0:
                postDic["typeFile"] = "inlet"
            # print (postDic)
            postFlowVec.append(postDic)

        self.facesPost = postFlowVec

        return

    def getPostFiles(self, fPath):
        """
        this function crawl the folder to reach the data in the post process file
        """

        filePaths = []
        folderItms = os.listdir(fPath)
        folderItms = sorted(folderItms, reverse=False, key=lambda x: float(x))

        for fld in folderItms:
            fPath1 = os.path.join(fPath, fld)
            if os.path.isdir(fPath1):
                fileItms = os.listdir(fPath1)
                completePath = os.path.join(fPath1, fileItms[0])
                filePaths.append(completePath)

        return filePaths

    def postFileArea(self, fPath):
        """
        this function extracts the face area (in m2) contained in the post
        processing file
        """

        try:
            postFile = open(fPath, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open postprocess file")
            sys.exit()
        with postFile:
            lines = postFile.readlines()
            for line in lines:
                line = line.replace("#", "")
                wrds = line.split()
                if "area" in wrds[0].lower():
                    area = float(wrds[-1])
                    break

        return area

    def postFileArray(self, fPathList, name):
        """
        this function extracts a generic array contained in the post
        processing file
        """

        arrayList = []

        for fPath in fPathList:
            array = []
            index = -1

            try:
                postFile = open(fPath, "r", encoding="utf8", errors="ignore")
            except IOError:
                print("couldn't open postprocess file")
                sys.exit()
            with postFile:
                lines = postFile.readlines()
                for line in lines:
                    line = line.replace("#", "")

                    wrds = line.split()

                    if len(wrds) == 0:
                        continue

                    if "time" in wrds[0].lower():
                        index = wrds.index(name)
                        continue

                    if index != -1:
                        array.append(float(wrds[index]))

                arrayList.append(array)

        return arrayList

    def readPostProcess_Tflows(self):
        """
        this function reads the post process folders and extract from all
        the folders the flow of the T scalar
        """

        postFolder = os.path.join(self.fluned_path, "postProcessing")

        # read the files in the postprocess folder
        postTFlowVec = []

        folderItms = os.listdir(postFolder)
        flowFolders = [itm for itm in folderItms if itm[0:9] == "volTFlow-"]

        for itm in flowFolders:
            postDic = {}
            postDic["faceID"] = itm[9:]
            postDic["folder"] = itm
            postDic["flowFiles"] = self.getPostFiles(os.path.join(postFolder, itm))
            timeLists = self.postFileArray(postDic["flowFiles"], "Time")

            flowLists = self.postFileArray(postDic["flowFiles"], "sum(T)")

            timeListSorted, flowListSorted = merge_continue_runs(
                timeLists,
                flowLists,
            )
            postDic["timeFile"] = timeListSorted
            postDic["TFlowFile"] = flowListSorted
            # postDic['TFlowFile']=postFileArray(postDic['flowFiles'],'sum(T)')
            postDic["TFlowFileLast"] = postDic["TFlowFile"][-1]
            # print (postDic)
            postTFlowVec.append(postDic)

        for face in self.facesPost:
            face["TFlowFileLast"] = [
                v["TFlowFileLast"]
                for v in postTFlowVec
                if v["faceID"] == face["faceID"]
            ][0]
            face["TFlowFile"] = [
                v["TFlowFile"] for v in postTFlowVec if v["faceID"] == face["faceID"]
            ][0]
            face["avTFileLast"] = face["TFlowFileLast"] / face["phiFlowFileLast"]
            face["avTFile"] = [
                t / f for t, f in zip(face["TFlowFile"], face["phiFlowFile"])
            ]

        return

    def readPostProcess_Trflows(self):
        """
        this function reads the post process folders and extract from all
        the folders the flow of the Tr scalar
        """

        postFolder = os.path.join(self.fluned_path, "postProcessing")

        # read the files in the postprocess folder
        postTrFlowVec = []

        folderItms = os.listdir(postFolder)
        flowFolders = [itm for itm in folderItms if itm[0:10] == "volTrFlow-"]

        if len(flowFolders) == 0:
            # if the output relative to time is are not present
            # initialize to zero the following variables
            for face in self.facesPost:
                face["avTrFile"] = [0]
                face["avTrGrad"] = [0]
                face["avTrFrac"] = 0
                face["rtdResTime"] = 0
                face["rtdDecRate"] = 0
            return

        for itm in flowFolders:
            postDic = {}
            postDic["faceID"] = itm[10:]
            postDic["folder"] = itm
            postDic["flowFiles"] = self.getPostFiles(os.path.join(postFolder, itm))

            timeLists = self.postFileArray(postDic["flowFiles"], "Time")

            flowLists = self.postFileArray(postDic["flowFiles"], "sum(Tr)")

            timeListSorted, flowListSorted = merge_continue_runs(
                timeLists,
                flowLists,
            )

            postDic["TrFlowFile"] = flowListSorted
            # postDic['TrFlowFile']=postFileArray(
            #        postDic['flowFiles'],'sum(Tr)')
            postDic["TrFlowFileLast"] = postDic["TrFlowFile"][-1]
            # print (postDic)
            postTrFlowVec.append(postDic)

        for face in self.facesPost:
            face["TrFlowFileLast"] = [
                v["TrFlowFileLast"]
                for v in postTrFlowVec
                if v["faceID"] == face["faceID"]
            ][0]
            face["TrFlowFile"] = [
                v["TrFlowFile"] for v in postTrFlowVec if v["faceID"] == face["faceID"]
            ][0]
            face["avTrFileLast"] = face["TrFlowFileLast"] / face["phiFlowFileLast"]
            face["avTrFile"] = [
                t / f for t, f in zip(face["TrFlowFile"], face["phiFlowFile"])
            ]
            face["avTrGrad"] = list(
                (np.gradient(face["avTrFile"], face["timeFile"], edge_order=1))
            )

            dt = face["timeFile"][1] - face["timeFile"][0]

            face["avTrFrac"] = [dt * g for g in face["avTrGrad"]]

            face["rtdResTime"] = sum(
                [t * g for t, g in zip(face["timeFile"], face["avTrFrac"])]
            )

            face["rtdDecRate"] = sum(
                [
                    g * math.exp(-self.decay_constant * t)
                    for t, g in zip(face["timeFile"], face["avTrFrac"])
                ]
            )

        return

    def write_summary(self, arguments):
        """
        write a bunch of results of the simulation
        """

        print("printing final Summary ... ")

        if self.time_treatment == "steadystate":
            self.write_summary_steady(arguments)
        elif self.time_treatment == "transient":
            self.write_summary_transient(arguments)

        return

    def write_summary_transient(self, arguments):
        sumFile = os.path.join(self.results_folder, "SUMMARY.csv")

        inletAtoms = 0
        inletFlow = 0
        outletAtoms = 0
        outletFlow = 0

        self.facesPost = sorted(self.facesPost, key=lambda x: x["faceID"])

        for face in self.facesPost:
            if face["typeFile"] == "inlet":
                inletAtoms += abs(face["TFlowFileLast"])
                inletFlow += abs(face["phiFlowFileLast"])

            if face["typeFile"] == "outlet":
                outletAtoms += abs(face["TFlowFileLast"])
                outletFlow += abs(face["phiFlowFileLast"])

        totalInletConc = inletAtoms / inletFlow
        totalInActivity = totalInletConc * self.decay_constant

        totalOutletConc = outletAtoms / outletFlow
        totalOutActivity = totalOutletConc * self.decay_constant

        cellActivity = np.zeros(len(self.TScalar))

        for i in range(len(self.TScalar)):
            if self.TScalar[i] > 0:
                cellActivity[i] = (
                    self.TScalar[i] * self.Volumes[i] * self.decay_constant
                )
            else:
                cellActivity[i] = 0

        totActivity = sum(cellActivity)
        avgActivity = totActivity / sum(self.Volumes)

        averageVolume = sum(self.Volumes) / len(self.Volumes)

        negCheck = any(n < 0 for n in self.TScalar)

        with open(sumFile, "w") as fw:
            fw.write("FLUNED SIMULATION SUMMARY\n")
            fw.write("CASE,{},\n".format(self.case))
            fw.write("TRANSIENT SIMULATION\n")
            fw.write("N ELEMENTS,{},\n".format(len(self.TScalar)))
            fw.write("ISOTOPE,{},\n".format(self.isotope.upper()))
            fw.write("DECAY CONSTANT,{:e},\n".format(self.decay_constant))
            fw.write("MOL DIFFUSION,{:e},\n".format(self.molecular_dif))
            fw.write("TURB SCHMIDT N,{:f},\n".format(self.schmidt_number))
            fw.write("\n")
            fw.write("\n")
            fw.write("QUALITY\n")
            fw.write("AVG VOL [m3],{:e},\n".format(averageVolume))

            if negCheck:
                fw.write("WARNING some elements are negative")

            if arguments.cdgs:
                fw.write("\n")
                fw.write("\n")
                fw.write("SOURCE SAMPLING\n")
                fw.write("SAMPLING RESOLUTION [m],{:f},\n".format(self.precision))
                fw.write(
                    "SAMPLED VOXELS [#],{:d},\n".format(
                        self.xInts * self.yInts * self.zInts
                    )
                )
                fw.write(
                    "VTK EMISSION RATE [#/s],{:e},\n".format(self.originalEmissionRate)
                )
                sampString = "SAMPLED EMISSION RATE (UNSCALED) [#/s],{:e},\n"
                fw.write(sampString.format(self.unscaledEmissionRate))

            fw.write("\n")
            fw.write("\n")
            fw.write("ACTIVATION\n")
            fw.write("INLET ATOMS FINAL [#/s],{:.5e},\n".format(inletAtoms))
            fw.write("OUTLET ATOMS FINAL [#/s],{:.5e},\n".format(outletAtoms))
            fw.write("TOT IN ACTIVITY FINAL [Bq/m3],{:.5e},\n".format(totalInActivity))
            fw.write(
                "TOT OUT ACTIVITY FINAL [Bq/m3],{:.5e},\n".format(totalOutActivity)
            )
            fw.write("TOT CASE ACTIVITY FINAL [Bq],{:.5e},\n".format(totActivity))
            fw.write("TOT AVG ACTIVITY FINAL [Bq/m3],{:.5e},\n".format(avgActivity))
            if totalInActivity != 0:
                normAvgActivity = avgActivity / totalInActivity
                stri = "INLET-NORMALIZED AVG ACTIVITY FINAL [Bq/m3],{:.5e},\n"
                fw.write(stri.format(normAvgActivity))

            if inletAtoms != 0:
                reductionRate = outletAtoms / inletAtoms
                fw.write("OUT/IN RATIO FINAL,{:.5e},\n".format(reductionRate))

            fw.write("\n")
            fw.write("\n")
            fw.write("FACES\n")

            for face in self.facesPost:
                if face["typeFile"] != "wall":
                    fw.write(face["faceID"])
                    fw.write("\n")
                    fw.write("TYPE,{},\n".format(face["typeFile"]))
                    fw.write("AREA [m2],{:.5e},\n".format(face["areaFile"]))
                    fw.write(
                        "FLUID FLOW FINAL [m3/s],{:.5e},\n".format(
                            face["phiFlowFileLast"]
                        )
                    )
                    fw.write(
                        "ATOM FLOW FINAL [#/s],{:.5e},\n".format(face["TFlowFileLast"])
                    )
                    fw.write(
                        "ATOM CONC FINAL [#/m3],{:.5e},\n".format(face["avTFileLast"])
                    )
                    fw.write(
                        "SPECIFIC ACTIVITY FINAL [Bq/m3],{:.5e},\n".format(
                            abs(face["avTFileLast"]) * self.decay_constant
                        )
                    )
                    if face["typeFile"] == "outlet":
                        fluxFrac = abs(face["phiFlowFileLast"] / outletFlow)
                        fw.write(
                            "AVG RTD RES T [s],{:.5f},\n".format(
                                abs(face["rtdResTime"])
                            )
                        )
                        fw.write(
                            "RTD FACE RED RATE,{:.5f},\n".format(
                                fluxFrac * (face["rtdDecRate"])
                            )
                        )
                    fw.write("\n")

        for face in self.facesPost:
            if face["typeFile"] == "wall":
                continue

            # write the flowing atoms at inlet-outlet
            faceSummary = (
                "face-atom-flow-" + face["faceID"] + "-" + face["typeFile"] + ".csv"
            )
            sumFile1 = os.path.join(self.results_folder, faceSummary)
            with open(sumFile1, "w") as fw:
                for t, c in zip(face["timeFile"], face["TFlowFile"]):
                    fw.write("{:.3f},{:.5e},\n".format(t, c))

            # write the concentration at inlet-outlet
            faceSummary = (
                "face-conc-" + face["faceID"] + "-" + face["typeFile"] + ".csv"
            )
            sumFile1 = os.path.join(self.results_folder, faceSummary)
            with open(sumFile1, "w") as fw:
                for t, c in zip(face["timeFile"], face["avTFile"]):
                    fw.write("{:.3f},{:.5e},\n".format(t, c))

            # write the specific activity at inlet-outlet
            faceSummary = (
                "face-specific-activity-"
                + face["faceID"]
                + "-"
                + face["typeFile"]
                + ".csv"
            )
            sumFile1 = os.path.join(self.results_folder, faceSummary)
            with open(sumFile1, "w") as fw:
                for t, c in zip(face["timeFile"], face["avTFile"]):
                    fw.write("{:.3f},{:.5e},\n".format(t, c * self.decay_constant))

            # write the specific time-conc at inlet-outlet
            faceSummary = (
                "face-fictitious-time-"
                + face["faceID"]
                + "-"
                + face["typeFile"]
                + ".csv"
            )
            sumFile1 = os.path.join(self.results_folder, faceSummary)
            with open(sumFile1, "w") as fw:
                for t, c in zip(face["timeFile"], face["avTrFile"]):
                    fw.write("{:.3f},{:.5e},\n".format(t, c))

            # write the RTD at inlet-outlet
            faceSummary = (
                "face-RTD-raw-" + face["faceID"] + "-" + face["typeFile"] + ".csv"
            )
            sumFile1 = os.path.join(self.results_folder, faceSummary)
            with open(sumFile1, "w") as fw:
                for t, g in zip(face["timeFile"], face["avTrGrad"]):
                    fw.write("{:.3f},{:.5e},\n".format(t, g))

        return

    def write_summary_steady(self, arguments):
        """
        This function writes the summary file in the RESULTS/ folder
        """

        summary_file = os.path.join(self.results_folder, "SUMMARY.csv")

        inletAtoms = 0
        inletFlow = 0
        outletAtoms = 0
        outletFlow = 0

        self.facesPost = sorted(self.facesPost, key=lambda x: x["faceID"])

        for face in self.facesPost:
            if face["typeFile"] == "inlet":
                inletAtoms += abs(face["TFlowFileLast"])
                inletFlow += abs(face["phiFlowFileLast"])

            if face["typeFile"] == "outlet":
                outletAtoms += abs(face["TFlowFileLast"])
                outletFlow += abs(face["phiFlowFileLast"])

        totalInletConc = inletAtoms / inletFlow
        totalInActivity = totalInletConc * self.decay_constant

        totalOutletConc = outletAtoms / outletFlow
        totalOutActivity = totalOutletConc * self.decay_constant

        cellActivity = np.zeros(len(self.TScalar))

        for i, _ in enumerate(self.TScalar):
            if self.TScalar[i] > 0:
                cellActivity[i] = (
                    self.TScalar[i] * self.Volumes[i] * self.decay_constant
                )
            else:
                cellActivity[i] = 0

        totActivity = sum(cellActivity)
        avgActivity = totActivity / sum(self.Volumes)

        averageVolume = sum(self.Volumes) / len(self.Volumes)

        # parameter 1 based on residence times
        # averageCellReTime = sum(self.resTimes)/len(self.resTimes)

        # qualityPar = [t*c*lbda/(ln2*totS) for t,c
        #              in zip(self.resTimes,self.TScalar)]

        # avgMeshQualParameter = sum(qualityPar)/len(qualityPar)

        if arguments.check:
            lbda = self.decay_constant
            ln2 = np.log(2)

            totS = sum(self.TScalar)

            # parameter 2 based on transit times

            transitTimes = [
                (vol ** (1 / 3)) / mod(v)
                for vol, v in zip(self.Volumes, self.Velocities)
            ]

            averageTransitTime = sum(transitTimes) / len(transitTimes)

            qualityPar2 = [
                t * c * lbda / (ln2 * totS) for t, c in zip(transitTimes, self.TScalar)
            ]

            avgMeshQualParameter2 = sum(qualityPar2) / len(qualityPar2)

            # parameter 3 based on scalar gradient

            gradients = [mod(g) for g in self.Gradients]

            # avgGrad = sum(gradients) / len(gradients)

            gradientDist = [
                g * (vol ** (1 / 3)) for vol, g in zip(self.Volumes, gradients)
            ]

            avgGradDist = sum(gradientDist) / len(gradientDist)

            qualityPar3 = [gd * c / totS for gd, c in zip(gradientDist, self.TScalar)]

            avgMeshQualParameter3 = sum(qualityPar3) / len(qualityPar3)

        negCheck = any(n < 0 for n in self.TScalar)

        with open(summary_file, "w") as fw:
            fw.write("FLUNED SIMULATION SUMMARY\n")
            fw.write("CASE,{},\n".format(self.case))
            fw.write("STEADY STATE SIMULATION\n")
            fw.write("N ELEMENTS,{},\n".format(len(self.TScalar)))
            fw.write("ISOTOPE,{},\n".format(self.isotope.upper()))
            fw.write("DECAY CONSTANT,{:e},\n".format(self.decay_constant))
            fw.write("MOL DIFFUSION,{:e},\n".format(self.molecular_dif))
            fw.write("TURB SCHMIDT N,{:f},\n".format(self.schmidt_number))
            fw.write("\n")
            fw.write("\n")
            fw.write("QUALITY\n")
            fw.write("AVG VOL [m3],{:e},\n".format(averageVolume))
            # fw.write("AVG CELL RES T [s],{:f},\n".format(averageCellReTime))
            # fw.write("AVG CELL Q1,{:e},\n".format(avgMeshQualParameter))
            if arguments.check:
                fw.write("AVG CELL TRNST T [s],{:f},\n".format(averageTransitTime))
                fw.write("AVG CELL Q2,{:e},\n".format(avgMeshQualParameter2))
                fw.write("AVG CELL GRADxLEN,{:e},\n".format(avgGradDist))
                fw.write("AVG CELL Q3,{:e},\n".format(avgMeshQualParameter3))
            if negCheck:
                fw.write("WARNING some elements are negative")
            if arguments.cdgs:
                fw.write("\n")
                fw.write("\n")
                fw.write("SOURCE SAMPLING\n")
                fw.write("SAMPLING RESOLUTION [m],{:f},\n".format(self.precision))
                fw.write(
                    "SAMPLED VOXELS [#],{:d},\n".format(
                        self.xInts * self.yInts * self.zInts
                    )
                )
                fw.write(
                    "VTK EMISSION RATE [#/s],{:e},\n".format(self.originalEmissionRate)
                )
                sampString = "SAMPLED EMISSION RATE (UNSCALED) [#/s],{:e},\n"
                fw.write(sampString.format(self.unscaledEmissionRate))
            fw.write("\n")
            fw.write("\n")
            fw.write("ACTIVATION\n")
            fw.write("INLET ATOMS [#/s],{:.5e},\n".format(inletAtoms))
            fw.write("OUTLET ATOMS [#/s],{:.5e},\n".format(outletAtoms))
            fw.write("TOT IN ACTIVITY [Bq/m3],{:.5e},\n".format(totalInActivity))
            fw.write("TOT OUT ACTIVITY[Bq/m3],{:.5e},\n".format(totalOutActivity))
            fw.write("TOT CASE ACTIVITY[Bq],{:.5e},\n".format(totActivity))
            fw.write("TOT AVG ACTIVITY[Bq/m3],{:.5e},\n".format(avgActivity))
            if totalInActivity != 0:
                normAvgActivity = avgActivity / totalInActivity
                string = "INLET-NORMALIZED AVG ACTIVITY[Bq/m3],{:.5e},\n"
                fw.write(string.format(normAvgActivity))

            if inletAtoms != 0:
                reductionRate = outletAtoms / inletAtoms
                fw.write("OUT/IN RATIO,{:.5e},\n".format(reductionRate))

            fw.write("\n")
            fw.write("\n")
            fw.write("FACES\n")

            for face in self.facesPost:
                if face["typeFile"] != "wall":
                    fw.write(face["faceID"])
                    fw.write("\n")
                    fw.write("TYPE,{},\n".format(face["typeFile"]))
                    fw.write("AREA [m2],{:.5e},\n".format(face["areaFile"]))
                    fw.write(
                        "FLUID FLOW [m3/s],{:.5e},\n".format(face["phiFlowFileLast"])
                    )
                    fw.write("ATOM FLOW [#/s],{:.5e},\n".format(face["TFlowFileLast"]))
                    fw.write("ATOM CONC [#/m3],{:.5e},\n".format(face["avTFileLast"]))
                    fw.write(
                        "SPECIFIC ACTIVITY [Bq/m3],{:.5e},\n".format(
                            abs(face["avTFileLast"]) * self.decay_constant
                        )
                    )
                    fw.write(
                        "AVG RES T [s],{:.5f},\n".format(abs(face["avTrFileLast"]))
                    )
                    fw.write("\n")

        return

    def calculateSamplingCoordinates(self):
        """
        self explainatory
        """
        mesh = pv.read(self.vtk_path)

        vtkBoundaries = mesh.bounds

        xBounds = [
            math.floor(vtkBoundaries[0] * 100),
            math.ceil(vtkBoundaries[1] * 100),
        ]
        yBounds = [
            math.floor(vtkBoundaries[2] * 100),
            math.ceil(vtkBoundaries[3] * 100),
        ]
        zBounds = [
            math.floor(vtkBoundaries[4] * 100),
            math.ceil(vtkBoundaries[5] * 100),
        ]

        samplingResolution = self.precision

        xInts = math.ceil((xBounds[1] - xBounds[0]) / (samplingResolution * 100))
        yInts = math.ceil((yBounds[1] - yBounds[0]) / (samplingResolution * 100))
        zInts = math.ceil((zBounds[1] - zBounds[0]) / (samplingResolution * 100))

        self.xInts = xInts
        self.yInts = yInts
        self.zInts = zInts

        xNodes = [(xBounds[0] + samplingResolution * 100 * i) for i in range(xInts + 1)]
        yNodes = [(yBounds[0] + samplingResolution * 100 * i) for i in range(yInts + 1)]
        zNodes = [(zBounds[0] + samplingResolution * 100 * i) for i in range(zInts + 1)]

        self.xNodes = xNodes
        self.yNodes = yNodes
        self.zNodes = zNodes

        voxelVector = []
        sampleCoordinates = []
        id = 1

        for i in range(xInts):
            # xVoxelNodes = [xNodes[i], xNodes[i + 1]]
            xVoxelCenter = (xNodes[i + 1] + xNodes[i]) / 2
            for j in range(yInts):
                # yVoxelNodes = [yNodes[j], yNodes[j + 1]]
                yVoxelCenter = (yNodes[j + 1] + yNodes[j]) / 2
                for k in range(zInts):
                    # zVoxelNodes = [zNodes[k], zNodes[k + 1]]
                    zVoxelCenter = (zNodes[k + 1] + zNodes[k]) / 2
                    newDic = {}
                    newDic["ID"] = id
                    newDic["index"] = [i, j, k]
                    newDic["centCoords"] = [xVoxelCenter, yVoxelCenter, zVoxelCenter]
                    sampleCoordinates.append(
                        [cord / 100 for cord in newDic["centCoords"]]
                    )
                    voxelVector.append(newDic)
                    id += 1

        self.voxelVector = voxelVector
        self.sampleCoordinates = sampleCoordinates
        self.voxelVolume = (self.precision * 100) ** 3  # cubic cc

        return

    def get_total_decays(self):
        """
        this function reads the vtk and calulates the total conc
        """

        print("Calculating total decay rate ... ")

        decay_rate_elements = self.get_decays(
            self.vtk_path,
            self.dataset,
            self.decay_constant,
            self.branching_ratio,
            self.scaling,
        )

        self.originalEmissionRate = sum(decay_rate_elements)

        return

    def get_decays(
        self, vtkFile, datasetName, decay_constant, branching_ratio, scaling
    ):
        """
        this function reads the vtk and calulates the total conc
        """

        checkFile = os.path.isfile(vtkFile)
        if not checkFile:
            print("ERROR vtk file not found")
            sys.exit()

        # read the vtk file with an unstructured grid
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(vtkFile)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()

        elementConc = data.GetCellData().GetArray(datasetName)
        elementConcArray = VN.vtk_to_numpy(elementConc)

        # pyvista part
        mesh = pv.read(vtkFile)
        mesh = mesh.compute_cell_sizes()
        volume_data = mesh.cell_data["Volume"]

        decay_rate_elements = [
            conc * vol * decay_constant * branching_ratio * scaling
            for conc, vol in zip(elementConcArray, volume_data)
        ]

        return decay_rate_elements

    def writeCDGS(self):
        """
        this function write the sample CDGS file
        """

        cdgsFile = self.vtk_path[:-3] + "CDGS"

        with open(cdgsFile, "w") as fw:
            fw.write("num_meshes 1\n")
            fw.write("global_source {:e}\n".format(self.scaledEmissionRate))
            fw.write("mesh_id 1\n")
            fw.write("cdgs from vtk\n")
            fw.write("Cooling_time 0.0\n")
            fw.write("total_source {:e}\n".format(self.scaledEmissionRate))
            fw.write("energy_type {}\n".format("bins"))
            fw.write("energy_boundaries {:d}\n".format(len(self.e_bins)))
            # WRITE SPECTRUM BINS
            specString = formatValues(self.e_bins)
            fw.write(specString)

            fw.write("mesh_type rec\n")
            fw.write(
                "mesh_boundaries {:d} {:d} {:d}\n".format(
                    self.xInts + 1, self.yInts + 1, self.zInts + 1
                )
            )
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
            specErrorString = formatValues([0] * (len(self.p_bins)))

            for vox in self.voxelVector:
                if vox["emission"] > 0:
                    fw.write(
                        voxelString1.format(
                            vox["ID"], vox["emission"], self.voxelVolume
                        )
                    )
                    fw.write(voxelString2.format(vox["emission"]))
                    emittingSpectrum = [val * vox["emission"] for val in self.p_bins]
                    spectrumString = formatValues(emittingSpectrum)
                    fw.write(spectrumString)
                    fw.write(specErrorString)

            fw.write("end_source_data\n")

        return

    def writeOpenMCSource(self):
        """
        this function write the mesh-based radiation source file
        """

        from lxml import etree as et

        openmc_source_file = os.path.join(
            self.results_folder, "structgrid_fluned_source.xml"
        )
        openmc_source_file_name = os.path.basename(openmc_source_file)
        openmc_source_import_commands = os.path.join(
            self.results_folder, "openmc_commands.txt"
        )

        # now it is hardcoded, later I will find a better way to handle this
        mesh_id = 100

        # dimensions = [self.xInts, self.yInts, self.zInts]

        mesh_lower_left = [
            min(self.xNodes),
            min(self.yNodes),
            min(self.zNodes),
        ]
        mesh_upper_right = [max(self.xNodes), max(self.yNodes), max(self.zNodes)]

        # create root element
        root = et.Element("source")

        # create sublement with the mesh source
        source_mesh = et.SubElement(
            root,
            "source",
            type="mesh",
            strength=str(self.scaledEmissionRate),
            mesh=str(mesh_id),
        )

        arr = np.asarray(self.voxelVector, dtype=object)  # keep objects untouched
        reordered = arr.reshape((self.xInts, self.yInts, self.zInts), order="C").ravel(
            order="F"
        )

        ebins_temp = [e * 1e6 for e in self.e_lines]  # convert from MeV to eV
        energy_parameters = [*ebins_temp, *self.p_lines]
        particle_type = self.particle_type

        for i, voxel in enumerate(reordered):
            sub_source = et.SubElement(
                source_mesh,
                "source",
                type="independent",
                strength=str(voxel["emission"]),
                particle=particle_type,
            )

            et.SubElement(sub_source, "angle", type="isotropic")
            energy = et.SubElement(sub_source, "energy", type="discrete")
            params = et.SubElement(energy, "parameters")
            params.text = " ".join(map(str, energy_parameters))

        # create mesh element with id attribute
        mesh = et.SubElement(root, "mesh", id=str(mesh_id))

        # add child elements with text content
        dimension = et.SubElement(mesh, "dimension")
        dimension.text = "{} {} {}".format(self.xInts, self.yInts, self.zInts)

        lower_left = et.SubElement(mesh, "lower_left")
        lower_left.text = "{} {} {}".format(
            mesh_lower_left[0], mesh_lower_left[1], mesh_lower_left[2]
        )

        upper_right = et.SubElement(mesh, "upper_right")
        upper_right.text = "{} {} {}".format(
            mesh_upper_right[0], mesh_upper_right[1], mesh_upper_right[2]
        )

        # write to file with xml declaration
        tree = et.ElementTree(root)
        tree.write(
            openmc_source_file,
            encoding="utf-8",
            pretty_print=True,
            xml_declaration=True,
        )

        with open(openmc_source_import_commands, "w") as f:
            f.write("from lxml import etree\n")
            f.write(
                'source_root = etree.parse("{}").getroot()\n'.format(
                    openmc_source_file_name
                )
            )
            f.write("mesh_element = source_root.find('mesh')\n")
            f.write("source_element = source_root.find('source')\n")
            f.write("mesh_geo = openmc.RegularMesh().from_xml_element(mesh_element)\n")
            f.write(
                "mesh_source = openmc.MeshSource.from_xml_element(source_element, {100:mesh_geo})\n"
            )

        return

    def write_sampled_source_vtk(self):
        """
        this function write the sampled source in a vtk file
        """

        vtkFile = os.path.join(self.results_folder, "sampled_source.vtk")

        dims = (self.xInts + 1, self.yInts + 1, self.zInts + 1)

        spacing = (
            (self.xNodes[-1] - self.xNodes[0]) / self.xInts,
            (self.yNodes[-1] - self.yNodes[0]) / self.yInts,
            (self.zNodes[-1] - self.zNodes[0]) / self.zInts,
        )

        origin = (self.xNodes[0], self.yNodes[0], self.zNodes[0])

        raw = np.asarray([vox["emission"] for vox in self.voxelVector], dtype=float)
        reordered = raw.reshape((self.xInts, self.yInts, self.zInts), order="C").ravel(
            order="F"
        )

        grid = pv.ImageData(dimensions=dims, spacing=spacing, origin=origin)
        grid.cell_data["emission_rate"] = reordered
        grid.save(vtkFile)

        return

    def write_openmc_um_source(self):
        """
        This function write the unstructured mesh-based radiation source file
        based on the vtk file
        """

        print("writing openmc unstructured mesh source file ..")

        from lxml import etree as et
        import meshio

        openmc_source_file = os.path.join(self.results_folder, "um_fluned_source.xml")
        h5m_basename = "um_geometry.h5m"
        openmc_source_mesh_file = os.path.join(self.results_folder, h5m_basename)
        vtk_intermediate_source_file = os.path.join(self.results_folder, "um_temp.vtk")
        vtk_intermediate_source_file_2 = os.path.join(
            self.results_folder, "um_temp2.vtk"
        )
        openmc_source_file_name = os.path.basename(openmc_source_file)
        openmc_source_import_commands = os.path.join(
            self.results_folder, "openmc_um_commands.txt"
        )

        # now it is hardcoded, later I will find a better way to handle this
        mesh_id = 100

        if self.vtk_path.lower().endswith(".vtu"):
            reader = vtk.vtkXMLUnstructuredGridReader()
        else:  # legacy ASCII/Binary *.vtk
            reader = vtk.vtkUnstructuredGridReader()
            reader.ReadAllVectorsOn()
            reader.ReadAllScalarsOn()
        reader.SetFileName(self.vtk_path)
        reader.Update()
        mesh = reader.GetOutput()
        mesh.GetPointData().Initialize()  # remove all point data
        cell_data = mesh.GetCellData()
        for i in reversed(range(cell_data.GetNumberOfArrays())):  # iterate safely
            if cell_data.GetArrayName(i) != "T":
                cell_data.RemoveArray(i)
        sx, sy, sz = 100, 100, 100

        tfm = vtk.vtkTransform()
        tfm.Scale(sx, sy, sz)

        tfilter = vtk.vtkTransformFilter()
        tfilter.SetTransform(tfm)
        tfilter.SetInputData(mesh)
        tfilter.Update()
        mesh_scaled = tfilter.GetOutput()

        tri = vtk.vtkDataSetTriangleFilter()
        tri.SetInputData(mesh_scaled)
        tri.SetTetrahedraOnly(True)
        tri.Update()
        mesh_tet = tri.GetOutput()

        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileTypeToBinary()
        writer.SetFileName(vtk_intermediate_source_file)
        if vtk.VTK_MAJOR_VERSION < 6:
            writer.SetInput(mesh_tet)
        else:
            writer.SetInputData(mesh_tet)
        writer.Write()

        # make a vtk file with no cell data to convert to h5m
        mesh_tet.GetCellData().Initialize()  # remove all cell data
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileTypeToBinary()
        writer.SetFileName(vtk_intermediate_source_file_2)
        if vtk.VTK_MAJOR_VERSION < 6:
            writer.SetInput(mesh_tet)
        else:
            writer.SetInputData(mesh_tet)
        writer.Write()

        meshio_object = meshio.read(vtk_intermediate_source_file_2)
        meshio.write(openmc_source_mesh_file, meshio_object)

        decay_rate_elements = self.get_decays(
            vtk_intermediate_source_file,
            self.dataset,
            self.decay_constant,
            self.branching_ratio,
            self.scaling,
        )

        print("length of the decay rate elements: ", len(decay_rate_elements))

        print(
            "consistency check, total emission rate from tetras in the vtk file: ",
            sum(decay_rate_elements),
        )
        print(
            "consistency check, total emission rate from the vtk file: ",
            self.originalEmissionRate,
        )

        with open(openmc_source_import_commands, "w") as f:
            f.write("from lxml import etree\n")
            f.write("parser = etree.XMLParser(huge_tree=True)\n")
            f.write(
                'source_root = etree.parse("{}", parser=parser).getroot()\n'.format(
                    openmc_source_file_name
                )
            )
            f.write("mesh_element = source_root.find('mesh')\n")
            f.write(
                "mesh_geo = openmc.UnstructuredMesh.from_xml_element(mesh_element)\n"
            )
            f.write("source_element = source_root.find('source')\n")
            f.write(
                "source = openmc.IndependentSource.from_xml_element(source_element, {100:mesh_geo})\n"
            )

        # create root element
        root = et.Element("source")

        # create sublement with the mesh source
        source_mesh = et.SubElement(
            root,
            "source",
            type="independent",
            particle=self.particle_type,
            strength=str(self.originalEmissionRate),
        )

        space = et.SubElement(
            source_mesh,
            "space",
            type="mesh",
            mesh_id=str(mesh_id),
            volume_normalized="False",
        )
        strengths = et.SubElement(space, "strengths")
        strengths.text = " ".join(
            map(str, [decay / 1e6 for decay in decay_rate_elements])
        )  # adjust that the decay rate has been calculated after scaling the vtk file

        # angle = et.SubElement(source_mesh, "angle", type="isotropic")

        ebins_temp = [e * 1e6 for e in self.e_lines]  # convert from MeV to eV
        energy_parameters = [*ebins_temp, *self.p_lines]
        energy = et.SubElement(source_mesh, "energy", type="discrete")
        params = et.SubElement(energy, "parameters")
        params.text = " ".join(map(str, energy_parameters))

        # create mesh element with id attribute
        mesh = et.SubElement(
            root,
            "mesh",
            id=str(mesh_id),
            name="source_mesh",
            type="unstructured",
            library="moab",
        )

        # add child elements with text content
        filename = et.SubElement(mesh, "filename")
        filename.text = h5m_basename

        # write to file with xml declaration
        tree = et.ElementTree(root)
        tree.write(
            openmc_source_file,
            encoding="utf-8",
            pretty_print=True,
            xml_declaration=True,
        )

        return

    def sampleCoordinatesVTK(self, vtkFile, datasetName, coordinates):
        """this function reads the vtk and sample the reaction rates"""

        checkFile = os.path.isfile(vtkFile)
        if not checkFile:
            print("ERROR vtk file not found")
            sys.exit()

        print("Sampling vtk ... ")

        # read the vtk file with an unstructured grid
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(vtkFile)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()

        # define probe

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(len(coordinates))

        for i, val in enumerate(coordinates):
            points.SetPoint(i, val[0], val[1], val[2])

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)

        # Perform the interpolation
        probeFilter = vtk.vtkProbeFilter()
        probeFilter.SetSourceData(data)
        probeFilter.SetInputData(polydata)
        probeFilter.Update()

        vtkArray = probeFilter.GetOutput().GetPointData().GetArray(datasetName)

        concentrations = VN.vtk_to_numpy(vtkArray)

        for i in range(len(concentrations)):  # remove eventual tiny negative conc
            if concentrations[i] < 0:
                concentrations[i] = 0

        return concentrations

    def sampleCoordinatesValues(self):
        """
        add the sampled value to the sample coordinates vector
        """

        vtkFile = self.vtk_path
        vtkDataSet = self.dataset

        sampledRates = self.sampleCoordinatesVTK(
            vtkFile, vtkDataSet, self.sampleCoordinates
        )

        totalEmission = 0

        voxelVolume = self.voxelVolume

        for voxel, concentration in zip(self.voxelVector, sampledRates):
            voxel["emission"] = (
                concentration
                * voxelVolume
                * self.decay_constant
                * self.branching_ratio
                * self.scaling
                * 1e-06
            )  # atoms per m3 to cm3

            totalEmission += voxel["emission"]

        self.unscaledEmissionRate = totalEmission

        ratioVtkSampling = self.originalEmissionRate / totalEmission

        totalEmissionScaled = 0

        for voxel, concentration in zip(self.voxelVector, sampledRates):
            voxel["emission"] = voxel["emission"] * ratioVtkSampling

            totalEmissionScaled += voxel["emission"]

        self.scaledEmissionRate = totalEmissionScaled

        return
