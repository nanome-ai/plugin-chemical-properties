import rdkit.Chem.Descriptors as Desc
import rdkit.Chem.rdMolDescriptors as mDesc
from .ESOLCalculator import ESOLCalculator

class PropertyCalculator():
    def __init__(self):
        self.esol = ESOLCalculator()
        self._properties = [
            ("MW", "Molecular Weight", "%.3f", Desc.MolWt),
            ("logP", "Lipophilicity (logP)", "%.3f", lambda mol: mDesc.CalcCrippenDescriptors(mol)[0]),
            ("TPSA", "Total Polar Surface Area", "%.3f", mDesc.CalcTPSA),
            ("ESOL", "Estimated Solubility", "%.3f", self.esol.calc_esol),
            ("HBA", "# H-Bond Acceptors", "%d", mDesc.CalcNumHBA),
            ("HBD", "# H-Bond Donors", "%d", mDesc.CalcNumHBD),
            ("RB", "# Rotatable Bonds", "%d", mDesc.CalcNumRotatableBonds),
            ("AR", "# Aromatic Rings", "%d", mDesc.CalcNumAromaticRings)
        ]

    def calc_property(self, mol, index):
        (short_lbl, long_lbl, fmt, fn) = self._properties[index]
        return (short_lbl, long_lbl, fmt % fn(mol))

    @property
    def num_props(self):
        return len(self._properties)

    @property
    def short_labels(self):
        return list(list(zip(*self._properties))[0])

    @property
    def long_labels(self):
        return list(list(zip(*self._properties))[1])
