import os
import copy


def get_post_process_list(path, dir_prefix, face_name, col_name):
    """
    this function gets the post process flow from the simulation folder
    """

    dir_name = dir_prefix + face_name

    flow_folders = [itm for itm in os.listdir(path) if itm == dir_name]

    if len(flow_folders) != 1:
        raise ValueError("Error with the number of flow folders")

    flow_files = get_post_files(os.path.join(path, flow_folders[0]))

    area_m2 = extract_area_field(flow_files[0])

    time_lists = post_file_array(flow_files, "Time")
    flow_lists = post_file_array(flow_files, col_name)

    time_list_sorted, flow_list_sorted = merge_continue_runs(
        time_lists,
        flow_lists,
    )

    return flow_list_sorted, time_list_sorted, area_m2


def merge_continue_runs(time_lists, data_lists):
    """
    this function takes time series segments and join them in a single one -
    the overlapping sections are removed
    """

    if len(time_lists) == 0 or len(data_lists) == 0:
        raise ValueError("Error with the number of post process files")

    if len(time_lists) != len(data_lists):
        print("ERROR mismatch in length of data series")
        print("number of time lists: ", len(time_lists))
        print("number of data lists: ", len(data_lists))
        raise ValueError("Error with the number of post process files")

    tot_len_time_lists = sum([len(x) for x in time_lists])
    tot_len_data_lists = sum([len(x) for x in data_lists])

    if tot_len_time_lists != tot_len_data_lists:
        print("ERROR mismatch in length of data series")
        print("number of time points: ", tot_len_time_lists)
        print("number of data points: ", tot_len_data_lists)
        raise ValueError("Error with the number of post process data points")

    temp_list = list(zip(time_lists, data_lists))

    # time_lists_sorted = sorted(time_lists,
    #                           key=lambda x:x[0])

    sorted_temp_list = sorted(temp_list, key=lambda x: x[0][0])

    time_lists_sorted, data_lists_sorted = zip(*sorted_temp_list)
    # data_lists_sorted = sorted(data_lists,
    #                           key=lambda x:x[0])

    time_series = time_lists_sorted[0]
    data_series = data_lists_sorted[0]

    if len(time_lists) > 1:
        for i, data_list in enumerate(time_lists_sorted):
            time_series_temp = copy.deepcopy(time_series)
            if i == 0:
                continue
            appended = False
            for j, val in enumerate(time_series_temp):
                if val >= time_lists_sorted[i][0]:
                    time_series = time_series[0:j]
                    time_series.extend(data_list)
                    data_series = data_series[0:j]
                    data_series.extend(data_lists_sorted[i])
                    appended = True
                    break

            if not appended:
                time_series.extend(data_list)
                data_series.extend(data_lists_sorted[i])

    return time_series, data_series


def get_post_files(path):
    """
    this function crawl the folder to reach the data in the post process file
    """

    file_paths = []
    folder_itms = os.listdir(path)
    folder_itms = sorted(folder_itms, reverse=False, key=lambda x: float(x))

    for fld in folder_itms:
        path_1 = os.path.join(path, fld)
        if os.path.isdir(path_1):
            file_itms = os.listdir(path_1)
            complete_path = os.path.join(path_1, file_itms[0])
            file_paths.append(complete_path)

    return file_paths


def post_file_array(path_list, name):
    """
    Extract a generic array contained in the post-processing files.
    If any file cannot be opened, immediately raise an exception
    instead of silently skipping it.
    """
    array_list = []

    for f_path in path_list:
        array = []
        index = -1

        try:
            with open(f_path, "r", encoding="utf8", errors="ignore") as post_file:
                for line in post_file:  # stream lines, don't load whole file
                    line = line.replace("#", "")
                    wrds = line.split()

                    if not wrds:
                        continue
                    if "time" in wrds[0].lower():
                        index = wrds.index(name)
                        continue
                    if index != -1:
                        array.append(float(wrds[index]))
        except IOError as e:
            raise IOError(f"couldn't open post-process file: {f_path}") from e

        array_list.append(array)

    return array_list


def extract_area_field(file_path):
    """
    This function reads a postprocess file and
    if the area value is present it returns it
    """

    area = 0

    try:
        postFile = open(file_path, "r", encoding="utf8", errors="ignore")
    except IOError:
        raise FileNotFoundError("couldn't open postprocess file")
    with postFile:
        lines = postFile.readlines()
        for line in lines:
            line = line.replace("#", "")
            wrds = line.split()
            if "area" in wrds[0].lower():
                area = float(wrds[-1])
                break

    return area


class PatchClass:
    """
    this class represents a patch of the mesh of an openfoam simulation
    """

    def __init__(self, patch_dict, path):
        """
        Constructs all the necessary attributes for the PatchClass object.

        Parameters:
        -----------
        patch_dict : dict
            A dictionary with the patch information.
        path : str
            The complete path to the simulation folder.
        """

        self.simulation_path = path
        self.post_process_path = os.path.join(self.simulation_path, "postProcessing")
        self.face_id = patch_dict["face_id"]
        self.face_elements_n = patch_dict["face_elements_n"]
        self.face_first_element = patch_dict["face_first_element"]
        self.face_type = patch_dict["type"]
        self.face_sum_phis = patch_dict["sum_phis"]
        self.post_process_flow = []
        self.post_process_time = []
        self.post_process_t_flow = []
        self.post_process_ta_flow = []
        self.post_process_td_flow = []
        self.post_process_tr_flow = []
        self.t_conc_atoms_m = []
        self.ta_conc_atoms_m3 = []
        self.td_conc_atoms_m3 = []
        self.tr_conc = []

        self.area_m2 = 0

    def post_process_face(self):
        """
        this function reads all the post process files and calculates the
        different scalar flows across the face
        """

        if self.face_type == "wall":
            return None

        self.post_process_flow, self.post_process_time, self.area_m2 = (
            get_post_process_list(
                self.post_process_path, "volFlow-", self.face_id, "sum(phi)"
            )
        )
        self.post_process_t_flow, _, _ = get_post_process_list(
            self.post_process_path, "volTFlow-", self.face_id, "sum(T)"
        )
        self.post_process_ta_flow, _, _ = get_post_process_list(
            self.post_process_path, "volTaFlow-", self.face_id, "sum(Ta)"
        )
        self.post_process_td_flow, _, _ = get_post_process_list(
            self.post_process_path, "volTdFlow-", self.face_id, "sum(Td)"
        )
        self.post_process_tr_flow, _, _ = get_post_process_list(
            self.post_process_path, "volTrFlow-", self.face_id, "sum(Tr)"
        )

        self.t_conc_atoms_m3 = [
            x / y if y != 0 else 0
            for x, y in zip(self.post_process_t_flow, self.post_process_flow)
        ]
        self.ta_conc_atoms_m3 = [
            x / y if y != 0 else 0
            for x, y in zip(self.post_process_ta_flow, self.post_process_flow)
        ]
        self.td_conc_atoms_m3 = [
            x / y if y != 0 else 0
            for x, y in zip(self.post_process_td_flow, self.post_process_flow)
        ]
        self.tr_conc = [
            x / y if y != 0 else 0
            for x, y in zip(self.post_process_tr_flow, self.post_process_flow)
        ]

        # check the sign of the flow
        if self.post_process_flow[-1] > 0 and self.face_type != "outlet":
            raise ValueError("Error with the flow sign")
        elif self.post_process_flow[-1] < 0 and self.face_type != "inlet":
            raise ValueError("Error with the flow sign")
        elif self.post_process_flow[-1] == 0 and self.face_type != "wall":
            raise ValueError("Error with the flow sign")
        if self.t_conc_atoms_m3[-1] < 0:
            raise ValueError("Error with the concentration sign")

        # print (self.post_process_t_flow)
        # print (self.t_conc_atoms_m3)

        return None
