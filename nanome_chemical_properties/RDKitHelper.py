from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import rdkit.Chem.Descriptors as Desc
import rdkit.Chem.rdMolDescriptors as mDesc
from .ESOLCalculator import ESOLCalculator

import tempfile
import shutil
from cairosvg import svg2png

# mol 2d image drawing options
Draw.DrawingOptions.atomLabelFontSize = 40
Draw.DrawingOptions.dotsPerAngstrom = 100
Draw.DrawingOptions.bondLineWidth = 8

class RDKitHelper:
    def __init__(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)

        self.esol = ESOLCalculator()
        self._properties = [
            ('MW', 'Molecular Weight', '%.3f', Desc.MolWt),
            ('logP', 'Lipophilicity (logP)', '%.3f', lambda mol: mDesc.CalcCrippenDescriptors(mol)[0]),
            ('TPSA', 'Total Polar Surface Area', '%.3f', mDesc.CalcTPSA),
            ('ESOL', 'Estimated Solubility', '%.3f', self.esol.calc_esol),
            ('HBA', '# H-Bond Acceptors', '%d', mDesc.CalcNumHBA),
            ('HBD', '# H-Bond Donors', '%d', mDesc.CalcNumHBD),
            ('RB', '# Rotatable Bonds', '%d', mDesc.CalcNumRotatableBonds),
            ('AR', '# Aromatic Rings', '%d', mDesc.CalcNumAromaticRings)
        ]

    def __del__(self):
        shutil.rmtree(self.temp_dir.name)

    @property
    def num_props(self):
        return len(self._properties)

    @property
    def short_labels(self):
        return list(list(zip(*self._properties))[0])

    @property
    def long_labels(self):
        return list(list(zip(*self._properties))[1])

    # adds complex.rdmol
    def prepare_complex(self, complex):
        # get only the current conformer
        m = next(complex.molecules)
        m.move_conformer(m.current_conformer, 0)
        m.set_conformer_count(1)

        complex.io.to_sdf(self.temp_sdf.name)
        complex.rdmol = Chem.SDMolSupplier(self.temp_sdf.name)[0]
        if complex.rdmol is None:
            raise ValueError('rdkit was unable to parse the complex')

    # adds complex.properties
    def add_properties(self, complex):
        complex.properties = []
        for short_lbl, long_lbl, fmt, fn in self._properties:
            complex.properties.append((short_lbl, long_lbl, fmt % fn(complex.rdmol)))

    # adds complex.imag
    def add_image(self, complex):
        complex.image = tempfile.NamedTemporaryFile(delete=False, suffix='.png', dir=self.temp_dir.name).name

        mol = complex.rdmol
        Chem.AssignStereochemistryFrom3D(mol)
        AllChem.Compute2DCoords(mol)
        mol = Draw.rdMolDraw2D.PrepareMolForDrawing(mol)

        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(256, 192)
        drawer.drawOptions().additionalAtomLabelPadding = 0.3
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        svg = svg.replace('stroke-linecap:butt', 'stroke-linecap:round')

        svg2png(bytestring=svg, write_to=complex.image, output_width=1024, output_height=768)

    def add_smiles(self, complex):
        complex.smiles = Chem.MolToSmiles(complex.rdmol)
