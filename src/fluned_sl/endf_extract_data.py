"""
code to read the ENDF data - author Patrick Sauvan - UNED
"""
import sys
import os

DECAY_LIB_NAME = 'eaf_dec_20070.total'
MODULE_PATH = os.path.dirname(os.path.abspath(__file__))
DECAY_LIB_PATH = os.path.join(MODULE_PATH,DECAY_LIB_NAME)

from .endf   import ENDF, MF8


class ENDF_library(object):
    """
    ENDF data class
    """

    def __init__(self):
        """
        initialization
        """
        self.endf = ENDF()
        self.endf.open_endf(DECAY_LIB_PATH)
        self.endf.generate_file_index()

    def neutron_spectra(self,nuclide):
        """
        get the neutron spectra
        """
        #try:
        #    nuclide = character_ENDFname(za)
        #except:
        #    nuclide = za
        #MF (first argument indicates data block)
        #MT (second data indicates kind of data inside of the data block)
        file8=MF8()
        self.endf.read_endf(nuclide)
        decay_data = self.endf.copy_MT(8,457)

        file8.insert_MT(457,decay_data)
        # spectra = file8.get_neutron_decay()
        spectra = file8.get_neutron_spectra()
        print(spectra)

        sys.exit()
        # spectra[0] => Energy    (eV)
        # spectra[1] => Intensity (gammas/bq)
        # spectra[2] => # MIO Programacion pone que es RTYP, pero yo creo que es STYP:
                            # 0 => Gamma
                            # 9 => x-ray


    #def heat_energy(self,za,styp=[0,9,1,2,4,5,6,7,8]):
    #    try:
    #        nuclide = character_ENDFname(za)
    #    except:
    #        nuclide = za

    #    #nuclide = character_ENDFname(za)
    #    #nuc_lines,Eavrg=endf.get_line_yield_heat(ZA,de)
    #    #MF (first argument indicates data block)
    #    #MT (second data indicates kind of data inside of the data block)
    #    file8=MF8()
    #    self.endf.read_endf(nuclide)
    #    decay_data = self.endf.copy_MT(8,457)

    #    file8.insert_MT(457,decay_data)
    #    # STYP
    #    # 0 gamma : 1 beta - : 2 E.C and beta+ : 4 alpha
    #    # 5 neutrons :  6 SF (Spontaneous fission): 7 protons
    #    # 8 Discrete electrons : 9X-ray
    #    #

    #    data = file8.get_avrE(styp=styp,detail=True)

    #    EM = data[0]
    #    EM = EM[1]   # Electro magnetic energy
    #    CP = data[1]
    #    CP = CP[1]   # Charged particle energy
    #    return EM,CP



    def gamma_spectra(self,nuclide):
        """
        function description
        """
        #try:
        #    nuclide = character_ENDFname(za)
        #except:
        #    nuclide = za
        #nuclide = character_ENDFname(za)

    #    nuc_lines,Eavrg=endf.get_line_yield_heat(ZA,de)
        #MF (first argument indicates data block)
        #MT (second data indicates kind of data inside of the data block)
        file8=MF8()
        self.endf.read_endf(nuclide)
        decay_data = self.endf.copy_MT(8,457)

        file8.insert_MT(457,decay_data)
        spectra = file8.get_gamma_spectra()

        # spectra[0] => Energy    (eV)
        # spectra[1] => Intensity (gammas/bq)
        # spectra[2] => # MIO Programacion pone que es RTYP, pero yo creo que es STYP:
                            # 0 => Gamma
                            # 9 => x-ray
        return spectra

    def gamma_yield(self,nuclide,delet=0.0):
        """
        function description
        """
        #try:
        #    nuclide = character_ENDFname(za)
        #except:
        #    nuclide = za

        #    nuc_lines,Eavrg=endf.get_line_yield_heat(ZA,de)
        file8=MF8()
        self.endf.read_endf(nuclide)
        decay_data = self.endf.copy_MT(8,457)

        file8.insert_MT(457,decay_data)
        spectra_tmp=file8.get_gamma_spectra()

        intensity = 0.0
        for line in spectra_tmp:
            if line[0]*line[1]*1e-6<delet:
                continue
            intensity = intensity + line[1]

        return intensity

    def  get_lambda(self,nuclide):
        """
        get the decay constant of a nuclide
        """
        #try:
        #    nuclide = character_ENDFname(za_val)
        #except:
        #    nuclide = za_val

        file8=MF8()
        self.endf.read_endf(nuclide)
        decay_data = self.endf.copy_MT(8,457)
        file8.insert_MT(457,decay_data)
        t_12=file8.get_halflife()
        if t_12 == 'stable':
            lmbda=0
        else:
            lmbda=0.69314718056/t_12
        return lmbda
