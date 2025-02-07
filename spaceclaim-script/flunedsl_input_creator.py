import os
import re
import sys

# tested with v24

# start by selecting the spaceclaim folder that contains the lines.

defaultTemp = 298  # Kelvin degrees

defaultPressure = 2  # MPa

first_surf_MCNP = 1

first_cell_MCNP = 1

def check_vector_orientation(vec1, vec2):
    """
    this function returns true if the two vectors are oriented in the same
    positive semi-space, false otherwise
    """

    dot_product = vec1.X * vec2.X + vec1.Y * vec2.Y + vec1.Z * vec2.Z

    if dot_product > 0:
        return True
    else:
        return False


def calculate_center_box(selection, padding_factor=0.1):


    options = WorkpieceOptions()
    options.WorkpieceType = WorkpieceType.Box
    options.WorkpieceCushion = BoxWorkpieceCushion(0)

    resultBox = Workpiece.Create(selection, options)
    box = resultBox.CreatedBodies[0]
    boxCompo = resultBox.GetCreated[IComponent]()
    BoxEdges = box.Edges
    BoxFaces = box.Faces

    lowZ = BoxFaces[0].GetFacePoint(0.5, 0.5)[2]
    highZ = BoxFaces[0].GetFacePoint(0.5, 0.5)[2]
    for face in BoxFaces:
        faceZ = face.GetFacePoint(0.5, 0.5)[2]
        if faceZ <= lowZ:
            lowZ = faceZ
            lowFace = face
        if faceZ >= highZ:
            highZ = faceZ
            highFace = face

    Vertices = []
    for edge in BoxEdges:
        Vertices.append(edge.EvalStart().Point)
        Vertices.append(edge.EvalEnd().Point)
    xpoints = [val[0] for val in Vertices]
    ypoints = [val[1] for val in Vertices]
    zpoints = [val[2] for val in Vertices]
    xmin = min(xpoints)
    xmax = max(xpoints)
    ymin = min(ypoints)
    ymax = max(ypoints)
    zmin = min(zpoints)
    zmax = max(zpoints)
    # print xextent
    # print yextent
    # print zextent
    x_extent = xmax - xmin
    y_extent = ymax - ymin
    z_extent = zmax - zmin


    x_extent = x_extent*(1+padding_factor)
    y_extent = y_extent*(1+padding_factor)
    z_extent = z_extent*(1+padding_factor)


    xcenter = (xmin + xmax)/2
    ycenter = (ymin + ymax)/2
    zcenter = (zmin + zmax)/2


    x_high = xcenter + x_extent/2
    x_low = xcenter - x_extent/2

    y_high = ycenter + y_extent/2
    y_low = ycenter - y_extent/2

    z_high = zcenter + z_extent/2
    z_low = zcenter - z_extent/2


    box_vector = [[x_low, x_high], [y_low, y_high], [z_low, z_high]]

    # delete box
    Delete.Execute(Selection.Create(boxCompo))

    return box_vector




def getParentComponents(body):
    level = body
    componentsVector = List[IDocObject]()

    while True:

        if level.Parent.Parent:
            level = level.Parent.Parent
            # print level.GetType()
            # print level.GetName()
            componentsVector.Add(level)
        else:

            level = level.Parent
            # print level.GetType()
            # print level.GetName()
            componentsVector.Add(level)
            break

    return componentsVector


def calculateDistance(point1, point2):
    distance = math.sqrt(
        (point1[0] - point2[0]) ** 2
        + (point1[1] - point2[1]) ** 2
        + (point1[2] - point2[2]) ** 2
    )
    return distance


def get_components(comp):
    """
    this function return the list of the
    first level components. If there are
    loose bodies in the second or first
    level it will produce an error
    """
    level_one_comp = comp.GetComponents()

    level_one_bodies = comp.GetBodies()
    if len(level_one_bodies) != 0:
        print("error, loose bodies at root level")
        sys.exit()

    return level_one_comp


def filter_pipes(bodies_list):
    real_pipes = List[IDesignBody]()
    for body in bodies_list:
        n_faces = 0
        n_cyls = 0
        body_faces = body.Faces
        is_pipe = True
        for face in body_faces:
            if isinstance(face.Shape.Geometry, Plane):
                n_faces += 1
            elif isinstance(face.Shape.Geometry, Cylinder):
                n_cyls += 1
            else:
                is_pipe = False
                break

        if is_pipe and n_faces == 2 and n_cyls == 1:
            real_pipes.Add(body)


    return real_pipes


def get_lines_info(comps_list):
    """
    this function create a list of dictionaries where each of these has the
    properties of one fluid cylinder
    """

    node_list = []

    line_dict = {}

    for comp in comps_list:

        bodies = comp.GetAllBodies()

        # find the LINE/TANK name - fill the LINE/TANK dictionaries
        body_parent_vector = getParentComponents(bodies[0])
        pat_line = re.compile("\A(\d+)_((LINE)|(TANK)|(CFD)|(STL))_", re.IGNORECASE)
        line_parents_pat = re.compile("\((.+?)\)")

        line_temp_pat = re.compile("_T(\d+(?:\.\d+)?)_?", re.IGNORECASE)
        line_pres_pat = re.compile("_P(\d+(?:\.\d+)?)_?", re.IGNORECASE)
        line_ext_pat = re.compile("ext(\d+)", re.IGNORECASE)

        # look for line/tank identifier
        for parent in body_parent_vector:
            parent_name = parent.GetName()
            water_line = pat_line.findall(parent_name)
            if water_line:
                # parentLine = parent
                water_line_name = parent.GetName()
                break


        if water_line:

            line_number = int(water_line[0][0])
            line_type = water_line[0][1].upper()
            line_temps = line_temp_pat.findall(water_line_name)
            line_pressures = line_pres_pat.findall(water_line_name)

            if len(line_temps) == 1:
                line_temp = float(line_temps[0])
            else:
                line_temp = defaultTemp

            if len(line_pressures) == 1:
                line_press = float(line_pressures[0])
            else:
                line_press = defaultPressure

            # add the line in the dictionary
            if line_number not in line_dict:

                line_dict[line_number] = {}
                line_dict[line_number]["fathers"] = []
                line_dict[line_number]["temperature"] = line_temp
                line_dict[line_number]["type"] = line_type
                line_dict[line_number]["pressure"] = line_press
                parent_groups = line_parents_pat.findall(water_line_name)
                if len(parent_groups) != 0:
                    for element in parent_groups:
                        definition = element.split(",")
                        while "" in definition:
                            definition.remove("")

                        if len(definition) == 1:

                            if "ext" in definition[0].lower():
                                exts = line_ext_pat.findall(definition[0])
                                parent_type = "ext"
                                parent_number = int(exts[0])
                                parent_fraction = 1.00
                            # print definition
                            else:
                                parent_type = "int"
                                parent_number = int(definition[0])
                                parent_fraction = 1.00
                        elif len(definition) == 2:
                            if "ext" in definition[0].lower():
                                exts = line_ext_pat.findall(definition[0])
                                parent_type = "ext"
                                parent_number = int(exts[0])
                                parent_fraction = float(definition[1])
                            else:
                                parent_type = "int"
                                parent_number = int(definition[0])
                                parent_fraction = float(definition[1])
                        else:
                            print("error in the parent line definition")
                            print(water_line_name)
                            sys.exit()
                        line_dict[line_number]["fathers"].append(
                            [parent_type, parent_number, parent_fraction]
                        )

                else:
                    print("error NO parent defined in the line")
                    print(water_line_name)
                    sys.exit()
        else:
            print("could not find the _LINE_ or _TANK_ or _CFD_ pattern")
            sys.exit()

        if line_type == "LINE":

            filtered_pipes = filter_pipes(bodies)

            if len(filtered_pipes) != len(bodies):
                print("error: in the line: ", water_line_name )
                print("not all pipes are cylinders")
                sys.exit()

        for b in bodies:

            # get cylinder information
            node_dict = {}

            # assign the component type
            if line_type == "TANK":
                node_dict["Type"] = "tank"
            else:
                node_dict["Type"] = "pipe"

            # assign the body type (tank components are allowed to be
            # bodies which are not cylinders)
            if line_type == "LINE":
                node_dict["body_type"] = "pipe"
            if line_type == "TANK":
                check_tank = filter_pipes([b])
                if len(check_tank) == 1:
                    node_dict["body_type"] = "tank-cyl"
                else:
                    node_dict["body_type"] = "tank"
            if line_type == "CFD":
                node_dict["body_type"] = "cfd"

            if line_type == "STL":
                node_dict["body_type"] = "stl"
                node_dict["body_object"] = b


            node_dict["Name"] = int(b.GetName())
            node_dict["Volume"] = b.MassProperties.Mass * 1000000
            node_dict["GroupName"] = line_type + "_" + str(line_number)
            node_dict["GroupLine"] = line_number
            node_dict["Temperature"] = float(line_temp)
            node_dict["Pressure"] = float(line_press)

            # calculate cyl geometry parameters
            if node_dict["body_type"] in ["pipe","tank-cyl"]:
                centers = []
                plane_normals = []
                d_values = []
                for face in b.Faces:
                    if isinstance(face.Shape.Geometry, Plane):
                        c_point = face.GetFacePoint(0.5, 0.5)
                        centers.append(c_point)
                        n_vec = face.GetFaceNormal(0.5, 0.5).UnitVector
                        plane_normals.append(n_vec)
                        d_val = (
                            n_vec.X * c_point.X * 100
                            + n_vec.Y * c_point.Y * 100
                            + n_vec.Z * c_point.Z * 100
                        )
                        d_values.append(d_val)

                    if isinstance(face.Shape.Geometry, Cylinder):
                        radius = face.Shape.Geometry.Radius
                        axis = face.Shape.Geometry.Axis.Direction.UnitVector

                node_dict["plane_centers"] = centers
                node_dict["plane_normals"] = plane_normals
                node_dict["plane_d_values"] = d_values

                node_dict["Radius"] = radius * 100
                node_dict["Axis"] = axis


                if check_vector_orientation(axis, plane_normals[0]):
                    new_axis = Direction.Create(-axis.X, -axis.Y, -axis.Z)
                    node_dict["Axis"] = new_axis

                node_dict["Length"] = ((calculateDistance(
                                    centers[0],centers[1])) * 100)

                # formula was obtained by expanding:
                #https://math.stackexchange.com/questions/1732385/cartesian-equation-cylinder-along-a-line

                # these coefficients are used to define the cylinder as GQ
                # surface for the MCNP input - However now the MCNP file is
                # going to be generated with GEOUNED so it will be removed

                x1 = centers[0].X * 100
                y1 = centers[0].Y * 100
                z1 = centers[0].Z * 100

                x2 = centers[1].X * 100
                y2 = centers[1].Y * 100
                z2 = centers[1].Z * 100

                r = radius * 100
                r2 = r ** 2

                x2mx1 = x2 - x1
                y2my1 = y2 - y1
                z2mz1 = z2 - z1

                xdiff2 = x2mx1 ** 2
                ydiff2 = y2my1 ** 2
                zdiff2 = z2mz1 ** 2

                a = y1*z2 - y2*z1
                b = x2*z1 - x1*z2
                c = x1*y2 - x2*y1

                a2 = a**2
                b2 = b**2
                c2 = c**2

                cyl_coefficients = {}

                cyl_coefficients["a"] = ydiff2 + zdiff2
                cyl_coefficients["b"] = xdiff2 + zdiff2
                cyl_coefficients["c"] = xdiff2 + ydiff2

                cyl_coefficients["d"] = -2 * x2mx1 * y2my1
                cyl_coefficients["e"] = -2 * y2my1 * z2mz1
                cyl_coefficients["f"] = -2 * x2mx1 * z2mz1

                cyl_coefficients["g"] = 2*(z2mz1*b - y2my1*c)
                cyl_coefficients["h"] = 2*(x2mx1*c - z2mz1*a)
                cyl_coefficients["j"] = 2*(y2my1*a - x2mx1*b)

                cyl_coefficients["k"] = (-r2*(xdiff2 + ydiff2 + zdiff2) +
                                         a2 +
                                         b2 +
                                         c2
                                            )

                node_dict["cyl_coefficients"] = cyl_coefficients
                # print (cyl_coefficients)



            node_list.append(node_dict)

    # find the last cell of the line, for each line
    for line in sorted(line_dict.keys()):
        cellVector = ([ element["Name"]
                        for element in node_list
                        if element["GroupLine"] == line ])

        if len(cellVector) == 0:
            print("error could not find the line:")
            print(line)
            sys.exit()

        line_dict[line]["lastCell"] = max(cellVector)

    # print CylList
    # print lineDictionary
    sort_list = sorted(node_list, key=lambda k: int(k["Name"]))

    # assign MCNP cell and surface numbers
    cell_count = first_cell_MCNP
    surf_count = first_surf_MCNP

    # a cell place-holder is set
    for element in sort_list:
        if element["body_type"] in ["pipe","tank-cyl"]:
            element["cell_number"] = cell_count
            cell_count += 1
            element["surf_number"] = surf_count
            surf_count += 3
        else:
            element["cell_number"] = cell_count
            element["surf_number"] = 1
            cell_count += 1

    return sort_list, line_dict



def write_circuit_mcnp(body_list, compList, name):
    """
    write the MCNP of the source circuit - only the cylinders


    NOT USED ANYMORE - MCNP input is generated with GEOUNED
    """
    box_list = []

    for comp in compList:

        bodies = comp.GetAllBodies()

        box_list.extend(bodies)

    selBodies = BodySelection.Create(box_list)

    box_vector = calculate_center_box(selBodies)

    # print (box_vector)

    try:
        path = GetRootPart().Document.Path[0:-6] + "-" + name + "-mcnp.dat"
    except:
        path = "D:\\pipes-error.dat"



    last_cell = max([element["cell_number"] for element in body_list])
    last_surf = max([element["surf_number"] for element in body_list])

    container_cell = last_cell + 1
    row_cell = last_cell + 2

    px_surf_low  = last_surf + 3
    px_surf_high = last_surf + 4
    py_surf_low  = last_surf + 5
    py_surf_high = last_surf + 6
    pz_surf_low  = last_surf + 7
    pz_surf_high = last_surf + 8

    with open(path, "w") as infile:

        infile.write("C MCNP model of the source circuit {:}\n".format(name))
        infile.write("C CELL CARDS\n")

        cyl_string = "{:<9d} {:d} -1.0 ( -{:d} -{:d} -{:d} ) \n"
        tag_string = "            IMP:N=1\n"
        tag_string_1 = "            IMP:N=0\n"
        comment_string = "C {:} --- {:} --- {:}\n"

        container_string = "{:<9d} 2 -0.0012 ({:d} -{:d} {:d} -{:d} {:d} -{:d})\n"
        neg_string = "         "
        row_string = "{:<9d} 0 (-{:d}:{:d}:-{:d}:{:d}:-{:d}:{:d})\n"
        for element in body_list:

            if element["body_type"] not in  ["pipe","tank-cyl"]:
                continue

            infile.write(
                cyl_string.format(
                    element["cell_number"],
                    1,
                    element["surf_number"],
                    element["surf_number"] + 1,
                    element["surf_number"] + 2,
                )
            )
            infile.write(tag_string)
            infile.write(comment_string.format(name,
                         element["GroupName"],
                         element["Name"])
                        )

        # write container cell
        infile.write(
            container_string.format(
                container_cell,
                px_surf_low,
                px_surf_high,
                py_surf_low,
                py_surf_high,
                pz_surf_low,
                pz_surf_high,
            )
        )
        # rough negation
        for cell in body_list:
            if cell["body_type"] not in  ["pipe","tank-cyl"]:
                continue

            if len(neg_string + " #" + str(cell["cell_number"])) > 79:
                infile.write(neg_string + "\n")
                neg_string = "         " + "#" + str(cell["cell_number"])
            else:
                neg_string += " #" + str(cell["cell_number"])

        infile.write(neg_string + "\n")

        infile.write(tag_string)
        infile.write("C container cell\n")

        # rest of world cell
        infile.write(
            row_string.format(
                row_cell,
                px_surf_low,
                px_surf_high,
                py_surf_low,
                py_surf_high,
                pz_surf_low,
                pz_surf_high,
            )
        )
        infile.write(tag_string_1)
        infile.write("C rest of world cell\n")

        infile.write("\n")
        infile.write("C SURFACE CARDS\n")
        cyl_surf_string_1 = "{:<9d} GQ {:<+20.15e} {:<+20.15e} {:<20.15e}\n"
        cyl_surf_string_2 = "             {:<+20.15e} {:<+20.15e} {:<20.15e}\n"
        cyl_surf_string_3 = "             {:<+20.15e}\n"
        plane_surf_string_1 = "{:<9d} P  {:<+20.15e} {:<+20.15e}\n"
        plane_surf_string_2 = "             {:<+20.15e} {:<+20.15e}\n"
        container_surf_string = "{:<9d} {} {:<+20.15e}\n"
        for element in body_list:
            if element["body_type"] not in  ["pipe","tank-cyl"]:
                continue

            infile.write(
                cyl_surf_string_1.format(
                    element["surf_number"],
                    element["cyl_coefficients"]["a"],
                    element["cyl_coefficients"]["b"],
                    element["cyl_coefficients"]["c"],
                )
            )
            infile.write(
                cyl_surf_string_2.format(
                    element["cyl_coefficients"]["d"],
                    element["cyl_coefficients"]["e"],
                    element["cyl_coefficients"]["f"],
                )
            )
            infile.write(
                cyl_surf_string_2.format(
                    element["cyl_coefficients"]["g"],
                    element["cyl_coefficients"]["h"],
                    element["cyl_coefficients"]["j"],
                )
            )
            infile.write(
                cyl_surf_string_3.format(
                    element["cyl_coefficients"]["k"],
                )
            )
            infile.write(
                plane_surf_string_1.format(
                    element["surf_number"] + 1,
                    element["plane_normals"][0].X,
                    element["plane_normals"][0].Y,
                )
            )
            infile.write(
                plane_surf_string_2.format(
                    element["plane_normals"][0].Z,
                    element["plane_d_values"][0],
                )
            )

            infile.write(
                plane_surf_string_1.format(
                    element["surf_number"] + 2,
                    element["plane_normals"][1].X,
                    element["plane_normals"][1].Y,
                )
            )
            infile.write(
                plane_surf_string_2.format(
                    element["plane_normals"][1].Z,
                    element["plane_d_values"][1],
                )
            )
        infile.write(
            container_surf_string.format(
                px_surf_low,
                "PX",
                box_vector[0][0] * 100,
            )
        )
        infile.write(
            container_surf_string.format(
                px_surf_high,
                "PX",
                box_vector[0][1] * 100,
            )
        )
        infile.write(
            container_surf_string.format(
                py_surf_low,
                "PY",
                box_vector[1][0] * 100,
            )
        )
        infile.write(
            container_surf_string.format(
                py_surf_high,
                "PY",
                box_vector[1][1] * 100,
            )
        )
        infile.write(
            container_surf_string.format(
                pz_surf_low,
                "PZ",
                box_vector[2][0] * 100,
            )
        )
        infile.write(
            container_surf_string.format(
                pz_surf_high,
                "PZ",
                box_vector[2][1] * 100,
            )
        )
        water_mat_definition = """C
c ***********************************
c *  WATER
c *  mass density [g/cc] - 1.0
c *  volume fraction [-] - 100
c *
c ***********************************
c *
M1       1001.31c 6.68466E-002 $ H  1 amount(-)  2.0000 ab(-) 99.99
         1002.31c 7.68824E-006 $ H  2 amount(-)  2.0000 ab(-)  0.02
         8016.31c 3.33459E-002 $ O 16 amount(-)  1.0000 ab(-) 99.757
         8017.31c 1.27023E-005 $ O 17 amount(-)  1.0000 ab(-)  0.0368
         8018.31c 6.85256E-005 $ O 18 amount(-)  1.0000 ab(-)  0.205
c *
c *  t.a.d. = 1.00281E-001
c *  eff.density = 1.00000E+000
C
c ***********************************
"""
        air_mat_definition = """C
c ***********************************
c *  DRY AIR
c *  mass density [g/cc] - 1.2E-003
c *  volume fraction [-] - 100
c *
c ***********************************
c
M2   6012.31c  7.38078E-09  $  C  12  WEIGHT(%)  0.0124  AB.(%)  98.93
      6013.31c  7.98285E-11  $  C  13  WEIGHT(%)  0.0124  AB.(%)  1.07
      7014.31c  3.88226E-05  $  N  14  WEIGHT(%)  75.5268  AB.(%)  99.632
      7015.31c  1.43395E-07  $  N  15  WEIGHT(%)  75.5268  AB.(%)  0.368
      8016.31c  1.04433E-05  $  O  16  WEIGHT(%)  23.1781  AB.(%)  99.757
      8017.31c  3.97814E-09  $  O  17  WEIGHT(%)  23.1781  AB.(%)  0.038
      8018.31c  2.14610E-08  $  O  18  WEIGHT(%)  23.1781  AB.(%)  0.205
      18036.31c  7.80801E-10  $  Ar  36  WEIGHT(%)  1.2827  AB.(%)  0.3365
      18038.31c  1.46647E-10  $  Ar  38  WEIGHT(%)  1.2827  AB.(%)  0.0632
      18040.31c  2.31109E-07  $  Ar  40  WEIGHT(%)  1.2827  AB.(%)  99.6003
c
c
c ***********************************
"""


        infile.write("\n")
        infile.write("C DATA CARDS\n")
        infile.write(water_mat_definition)
        infile.write(air_mat_definition)
        infile.write("MODE N\n")
        infile.write("SDEF\n")
        infile.write("NPS 10\n")

    return 0

def write_stl_nodes(node_list):
    """
    this function write the stl of the nodes represented with a stl
    """

    doc_path = GetRootPart().Document.Path
    path_base = os.path.dirname(doc_path) + '\\'

    for node in node_list:
        if node["body_type"] == 'stl':
            file_name = path_base + node["GroupName"] + "_" + str(node["Name"]) + ".stl"
            comp_bodies = BodySelection.Create([node["body_object"]])


            # copy the bodies to a new document
            DocumentHelper.CreateNewDocument()
            paste_result = Copy.Execute(comp_bodies)

            # Save File
            options = ExportOptions.Create()
            DocumentSave.Execute(file_name, options)
            options = STLExportOptions()

            options.ExportStlUnits = LengthUnits.M

            STLFile.Export(file_name, options)

            # close new design
            DocumentHelper.CloseDocument()




def write_pipes_dat(pipe_list, name):
    """
    this function write the pipes.dat file
    """

    path = GetRootPart().Document.Path[0:-6] + "-" + name + "-nodes.dat"

    with open(path, "w") as infile:


        header_data = [
            "Cell_Number",
            "Group_Name",
            "Node_Type",
            "Node_Type_Method",
            "Activity_Scaling",
            "Temperature_K",
            "Pressure_MPa",
            "MCNP_Cell",
            "Volume_cm3",
            "Reaction_Rate_cm-3",
            "Axis.X",
            "Axis.Y",
            "Axis.Z",
            "Origin.X_cm",
            "Origin.Y_cm",
            "Origin.Z_cm",
            "Length_cm",
            "Radius_cm",
        ]

        header_str = "{: >20}," * len(header_data) + "\n"
        data_str = (
            "{:>20d}," +
            "{:>20s}," +
            "{:>20s}," +
            "{:>+20d}," +
            "{:>+20d}," +
            "{:>20.5f}," +
            "{:>20.5f}," +
            "{:>20d}," +
            "{:>+20.5e}," +
            "{:>20d}," +
            "{:>+20.5e}," +
            "{:>+20.5e}," +
            "{:>+20.5e}," +
            "{:>+20.5e}," +
            "{:>+20.5e}," +
            "{:>+20.5e}," +
            "{:>+20.5e}," +
            "{:>+20.5e}," +
            "\n"
           )


        infile.write(header_str.format(*header_data))

        for element in pipe_list:

            if element['body_type'] == 'tank':
                write_str = data_str.format(
                    int(element["Name"]),
                    element["GroupName"],
                    element["body_type"],
                    int(-1),
                    int(1),
                    element["Temperature"],
                    element["Pressure"],
                    int(element["cell_number"]),
                    float(element["Volume"]),
                    int(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                )

            elif element['body_type'] == 'cfd':
                write_str = data_str.format(
                    int(element["Name"]),
                    element["GroupName"],
                    element["body_type"],
                    int(0),
                    int(1),
                    element["Temperature"],
                    element["Pressure"],
                    int(element["cell_number"]),
                    float(element["Volume"]),
                    int(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                )

            elif element['body_type'] == 'stl':
                write_str = data_str.format(
                    int(element["Name"]),
                    element["GroupName"],
                    element["body_type"],
                    int(-1),
                    int(1),
                    element["Temperature"],
                    element["Pressure"],
                    int(element["cell_number"]),
                    float(element["Volume"]),
                    int(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                    float(0),
                )

            elif element['body_type'] in ['tank-cyl','pipe']:
                write_str = data_str.format(
                    int(element["Name"]),
                    element["GroupName"],
                    element["body_type"],
                    int(-1),
                    int(1),
                    element["Temperature"],
                    element["Pressure"],
                    int(element["cell_number"]),
                    float(element["Volume"]),
                    int(0),
                    float(element["Axis"].X),
                    float(element["Axis"].Y),
                    float(element["Axis"].Z),
                    float(element["plane_centers"][0].X*100),
                    float(element["plane_centers"][0].Y*100),
                    float(element["plane_centers"][0].Z*100),
                    float(element["Length"]),
                    float(element["Radius"]),
                )

            infile.write(write_str)


def write_links_dat(nodes_list, lines_dict, name):
    """
    this function write the links.dat file that describes the connection
    between the circuit nodes
    """

    path = GetRootPart().Document.Path[0:-6] + "-" + name + "-links.dat"

    with open(path, "w") as infile:

        # write header
        header_data = ["Cell_Num", "n_Parents", "Parents"]
        header_str = "{: >10}," * len(header_data) + "\n"
        infile.write(header_str.format(*header_data))

        for line in sorted(lines_dict.keys()):

            fathers = lines_dict[line]["fathers"]

            # define and use a temporary list
            temp_line = [el for el in nodes_list if el["GroupLine"] == int(line)]
            sorted_temp_list = sorted(temp_line, key=lambda k: k["Name"])

            for i, element in enumerate(sorted_temp_list):
                if i == 0:
                    cell_str = "{: >10d},{: >10d},"
                    cell_vec = [element["Name"], len(fathers)]
                    for father in fathers:
                        cell_str += "{: >10},{: >10d},{: >10.5f},"
                        cell_vec.append(father[0])
                        if father[0] == "ext":
                            cell_vec.append(father[1])
                        if father[0] == "int":
                            fatherID = lines_dict[father[1]]["lastCell"]
                            cell_vec.append(fatherID)
                        cell_vec.append(father[2])

                    cell_str += "\n"
                    infile.write(cell_str.format(*cell_vec))

                else:
                    write_str = "{:>10d},{:>10d},\n".format(element["Name"], 1)
                    infile.write(write_str)
    return


def write_circuit_stl(CompList, name):
    """
    this function generate the stl of the geometry - SCALE IN METERS
    """

    path = GetRootPart().Document.Path[0:-6] + "-" + name + "-geometry_meters.stl"


    # select all bodies
    bodies_list = []
    for comp in CompList:
        bodies = comp.GetAllBodies()

        bodies_list.extend(bodies)
    comp_bodies = BodySelection.Create(bodies_list)


    # copy the bodies to a new document
    DocumentHelper.CreateNewDocument()
    paste_result = Copy.Execute(comp_bodies)


    # Scale
    bodies_scale = GetRootPart().GetAllBodies()
    new_sel = Selection.Create(bodies_scale)
    result = Scale.Execute(new_sel, Point.Create(0, 0, 0), 0.001, True)

    # Save File
    options = ExportOptions.Create()
    DocumentSave.Execute(path, options)

    # close new design
    DocumentHelper.CloseDocument()

    return 0


def write_inlets_dat(linesDict, name):
    """
    this function writes the file for the definition of the external inlets
    """

    path = GetRootPart().Document.Path[0:-6] + "-" + name + "-inlets.dat"

    with open(path, "w") as infile:

        # write header
        header_data = [
            "Ext_Source_ID",
            "Type",
            "Type-Value",
            "Isotope",
            "Isotope-Value",
        ]
        header_string = "{: >15}," * len(header_data) + "\n"
        infile.write(header_string.format(*header_data))

        for line in linesDict:
            for father in linesDict[line]["fathers"]:
                if father[0] == "ext":
                    cell_str = "{: >15}," * len(header_data)
                    inlet_id = str(father[1])
                    inlet_type = "source-node"

                    cell_vec = [
                        inlet_id,
                        inlet_type,
                        10,
                        "N16",
                        100,
                    ]
                    cell_str += "\n"
                    infile.write(cell_str.format(*cell_vec))

    return

def write_circuit_mcnp_test_file(body_list, compList, name):
    """
    write the MCNP for the enclosure of the cdgs
    """

    box_list = []

    voxel_side = 10 # cm

    for comp in compList:

        bodies = comp.GetAllBodies()

        box_list.extend(bodies)

    selBodies = BodySelection.Create(box_list)

    box_vector = calculate_center_box(selBodies, padding_factor=0.1)
    box_vector_fmesh = calculate_center_box(selBodies, padding_factor=0.05)


    # print (box_vector)

    path = GetRootPart().Document.Path[0:-6] + "-" + name + "-mcnp_source_test.i"



    last_cell = max([element["cell_number"] for element in body_list])
    last_surf = max([element["surf_number"] for element in body_list])

    n_cylinders = len(body_list)

    container_cell = last_cell + 1
    row_cell = last_cell + 2

    px_surf_low  = last_surf + 3
    px_surf_high = last_surf + 4
    py_surf_low  = last_surf + 5
    py_surf_high = last_surf + 6
    pz_surf_low  = last_surf + 7
    pz_surf_high = last_surf + 8

    with open(path, "w") as infile:

        infile.write("C MCNP test model for the source circuit {:}\n".format(name))
        infile.write("C CELL CARDS\n")

        #cyl_string = "{:<9d} {:d} -1.0 ( -{:d} -{:d} -{:d} ) \n"
        tag_string = "            IMP:N=1\n"
        tag_string_1 = "            IMP:N=0\n"
        #comment_string = "C {:} --- {:} --- {:}\n"

        container_string = "{:<9d} 2 -0.0012 ({:d} -{:d} {:d} -{:d} {:d} -{:d})\n"
        row_string = "{:<9d} 0 (-{:d}:{:d}:-{:d}:{:d}:-{:d}:{:d})\n"

        # write container cell
        infile.write(
            container_string.format(
                container_cell,
                px_surf_low,
                px_surf_high,
                py_surf_low,
                py_surf_high,
                pz_surf_low,
                pz_surf_high,
            )
        )


        infile.write(tag_string)
        infile.write("C container cell\n")

        # rest of world cell
        infile.write(
            row_string.format(
                row_cell,
                px_surf_low,
                px_surf_high,
                py_surf_low,
                py_surf_high,
                pz_surf_low,
                pz_surf_high,
            )
        )
        infile.write(tag_string_1)
        infile.write("C rest of world cell\n")

        infile.write("\n")
        infile.write("C SURFACE CARDS\n")
        #cyl_surf_string_1 = "{:<9d} GQ {:<+20.15e} {:<+20.15e} {:<20.15e}\n"
        #cyl_surf_string_2 = "             {:<+20.15e} {:<+20.15e} {:<20.15e}\n"
        #cyl_surf_string_3 = "             {:<+20.15e}\n"
        #plane_surf_string_1 = "{:<9d} P  {:<+20.15e} {:<+20.15e}\n"
        #plane_surf_string_2 = "             {:<+20.15e} {:<+20.15e}\n"
        container_surf_string = "{:<9d} {} {:<+20.15e}\n"

        infile.write(
            container_surf_string.format(
                px_surf_low,
                "PX",
                box_vector[0][0] * 100,
            )
        )
        infile.write(
            container_surf_string.format(
                px_surf_high,
                "PX",
                box_vector[0][1] * 100,
            )
        )
        infile.write(
            container_surf_string.format(
                py_surf_low,
                "PY",
                box_vector[1][0] * 100,
            )
        )
        infile.write(
            container_surf_string.format(
                py_surf_high,
                "PY",
                box_vector[1][1] * 100,
            )
        )
        infile.write(
            container_surf_string.format(
                pz_surf_low,
                "PZ",
                box_vector[2][0] * 100,
            )
        )
        infile.write(
            container_surf_string.format(
                pz_surf_high,
                "PZ",
                box_vector[2][1] * 100,
            )
        )
        water_mat_definition = """C
c ***********************************
c *  WATER
c *  mass density [g/cc] - 1.0
c *  volume fraction [-] - 100
c *
c ***********************************
c *
M1       1001.31c 6.68466E-002 $ H  1 amount(-)  2.0000 ab(-) 99.99
         1002.31c 7.68824E-006 $ H  2 amount(-)  2.0000 ab(-)  0.02
         8016.31c 3.33459E-002 $ O 16 amount(-)  1.0000 ab(-) 99.757
         8017.31c 1.27023E-005 $ O 17 amount(-)  1.0000 ab(-)  0.0368
         8018.31c 6.85256E-005 $ O 18 amount(-)  1.0000 ab(-)  0.205
c *
c *  t.a.d. = 1.00281E-001
c *  eff.density = 1.00000E+000
C
c ***********************************
"""
        air_mat_definition = """C
c ***********************************
c *  DRY AIR
c *  mass density [g/cc] - 1.2E-003
c *  volume fraction [-] - 100
c *
c ***********************************
c
M2   6012.31c  7.38078E-09  $  C  12  WEIGHT(%)  0.0124  AB.(%)  98.93
      6013.31c  7.98285E-11  $  C  13  WEIGHT(%)  0.0124  AB.(%)  1.07
      7014.31c  3.88226E-05  $  N  14  WEIGHT(%)  75.5268  AB.(%)  99.632
      7015.31c  1.43395E-07  $  N  15  WEIGHT(%)  75.5268  AB.(%)  0.368
      8016.31c  1.04433E-05  $  O  16  WEIGHT(%)  23.1781  AB.(%)  99.757
      8017.31c  3.97814E-09  $  O  17  WEIGHT(%)  23.1781  AB.(%)  0.038
      8018.31c  2.14610E-08  $  O  18  WEIGHT(%)  23.1781  AB.(%)  0.205
      18036.31c  7.80801E-10  $  Ar  36  WEIGHT(%)  1.2827  AB.(%)  0.3365
      18038.31c  1.46647E-10  $  Ar  38  WEIGHT(%)  1.2827  AB.(%)  0.0632
      18040.31c  2.31109E-07  $  Ar  40  WEIGHT(%)  1.2827  AB.(%)  99.6003
c
c
c ***********************************
"""


        infile.write("\n")
        infile.write("C DATA CARDS\n")
        infile.write(water_mat_definition)
        infile.write(air_mat_definition)
        infile.write("MODE P\n")
        infile.write("NPS 10\n")
        infile.write("IDUM {:d}  0 0  1 0\n".format(n_cylinders))
        infile.write("FILES 22  cdgs\n")
        infile.write("C CDGS test meshtal - FLUX \n")
        xmin = box_vector_fmesh[0][0] * 100
        xmax = box_vector_fmesh[0][1] *100
        xints = int(math.ceil((xmax - xmin) / voxel_side))
        ymin = box_vector_fmesh[1][0] * 100
        ymax = box_vector_fmesh[1][1] *100
        yints = int(math.ceil((ymax - ymin) / voxel_side))
        zmin = box_vector_fmesh[2][0] * 100
        zmax = box_vector_fmesh[2][1] *100
        zints = int(math.ceil((zmax - zmin) / voxel_side))
        infile.write("FMESH4:P GEOM=REC ORIGIN {:.1f} {:.1f} {:.1f}\n".format(xmin, ymin, zmin))
        infile.write("        IMESH   {:.1f}\n".format(xmax))
        infile.write("        IINTS   {:d}\n".format(xints))
        infile.write("        JMESH   {:.1f}\n".format(ymax))
        infile.write("        JINTS   {:d}\n".format(yints))
        infile.write("        KMESH   {:.1f}\n".format(zmax))
        infile.write("        KINTS   {:d}\n".format(zints))
        infile.write("C CDGS test meshtal - DOSE \n")
        infile.write("FMESH14:P GEOM=REC ORIGIN {:.1f} {:.1f} {:.1f}\n".format(xmin, ymin, zmin))
        infile.write("        IMESH   {:.1f}\n".format(xmax))
        infile.write("        IINTS   {:d}\n".format(xints))
        infile.write("        JMESH   {:.1f}\n".format(ymax))
        infile.write("        JINTS   {:d}\n".format(yints))
        infile.write("        KMESH   {:.1f}\n".format(zmax))
        infile.write("        KINTS   {:d}\n".format(zints))
        energy_function_string = """C
DE14   0.01 0.015 0.02 0.03  0.04 0.05
       0.06 0.07  0.08 0.10  0.15 0.20
       0.3  0.4   0.5  0.6   0.8  1.0
       2.0  4.0   6.0  8.0  10.0
DF14   0.0485  0.1254  0.2050  0.2999  0.3381 0.3572
       0.3780  0.4066  0.4399  0.5172  0.7523 1.0041
       1.5083  1.9958  2.4657  2.9082  3.7269 4.4834
       7.4896 12.0153 15.9873 19.9191 23.7600
"""
        infile.write("FM14        1\n")
        infile.write(energy_function_string)


    return 0

selectedComp = Selection.GetActive().GetItems[IComponent]()[0]
comp_name = selectedComp.GetName()

# get the component list below the selected components
print("get components ..")
comp_list = get_components(selectedComp)


# get data on the water cylinders
print("get cylinders data ..")
body_list, lines_dict = get_lines_info(comp_list)

# print linesDictionary

# order the list according to the body name
sorted_body_list = sorted(body_list, key=lambda k: int(k["Name"]))


# write pipes.dat
print("write pipes.DAT")
write_pipes_dat(sorted_body_list, comp_name)

# write links.dat
print("write links.DAT")
write_links_dat(body_list, lines_dict, comp_name)

# write  Inlets.dat
print("write inlets.DAT")
write_inlets_dat(lines_dict, comp_name)

#write stl nodes
write_stl_nodes(sorted_body_list)

# write geometry file in STL format - scale in meters
print("write STL file")
write_circuit_stl(comp_list, comp_name)

## write MCNP file of the source cylinders
## at the moment we write the MCNP file with GEOUNED
##print("write MCNP file")
##write_circuit_mcnp(sorted_body_list, comp_list, compName)

# write MCNP file to test the cdgs source
print("write MCNP enclosure file")
write_circuit_mcnp_test_file(sorted_body_list, comp_list, comp_name)


print("finished!")
