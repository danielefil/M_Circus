#classe Spettro


class Spectra:
    def __init__(self, Spectra):
        self.mz = Spectra[:, 0] 
        self.int = Spectra[:, 1]
        