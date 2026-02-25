from pathlib import Path

from .ofoam_base import oFoamBase


class oFoamFluned(oFoamBase):
    def __init__(self, path: str | Path):
        super().__init__(path)
        self.reduction_rate = []
        self.normalized_average_td = []
        self.inlet_td_conc_atoms_m3 = []
        self.outlet_rr_conc_atoms_m3 = []
        self.average_rr_conc_atoms_m3 = []
        self.post_process_td_average = []
        self.average_ta = []
        self.volume_m3 = 0
        self.vtk_file_path = ""
        self.scaled_vtk_file_path = ""
        self.vtk_dimensions = []
