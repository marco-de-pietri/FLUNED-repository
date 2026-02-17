import pickle

import h5py


def get_isotope_reactions_dicts(file_path):
    """
    this function  reads the depletion results and gets the dictionaries mapping the reactions and the parent isotopes
    """

    with h5py.File(file_path, "r") as f:
        nuc_idx = pickle.loads(f["index_nuc"][...].tobytes())
        rx_idx = pickle.loads(f["index_rx"][...].tobytes())

    return nuc_idx, rx_idx
