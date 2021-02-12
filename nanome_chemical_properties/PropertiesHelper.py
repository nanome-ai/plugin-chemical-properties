from nanome.util import Logs

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import rdkit.Chem.Descriptors as Desc
import rdkit.Chem.rdMolDescriptors as mDesc
from .ESOLCalculator import ESOLCalculator

import json
import os
import re
import requests
import shutil
import tempfile
from cairosvg import svg2png
from datetime import datetime, timedelta
from functools import partial

# mol 2d image drawing options
Draw.DrawingOptions.atomLabelFontSize = 40
Draw.DrawingOptions.dotsPerAngstrom = 100
Draw.DrawingOptions.bondLineWidth = 8

API_CACHE_TIME = timedelta(seconds=1)
API_SETTINGS = os.path.join(os.path.dirname(__file__), '..', 'config.json')

class PropertiesHelper:
    def __init__(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)

        self.api_cache = {}
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

        if not os.path.exists(API_SETTINGS):
            return

        with open(API_SETTINGS, 'r') as f:
            self.api = json.load(f)

        if self.api.get('overwrite'):
            self._properties = []

        # validate config
        try:
            required_endpoint_keys = ['url', 'method', 'data']
            required_property_keys = ['description', 'format', 'path']

            for endpoint in self.api.get('endpoints'):
                if not endpoint.get('name'):
                    raise Exception('Invalid config: missing endpoint name')

                test = [endpoint[k] for k in required_endpoint_keys]
                for prop, info in endpoint.get('properties').items():
                    test = [info[k] for k in required_property_keys]

        except KeyError as key:
            raise Exception(f'Invalid config: missing {key} on {endpoint["name"]}')

        except TypeError:
            raise Exception(f'Invalid config: array where object should be')

        # register properties
        for endpoint in self.api.get('endpoints'):
            for prop, info in endpoint['properties'].items():
                fn = partial(self.fetch_property, endpoint, prop)
                p = (prop, info['description'], info['format'], fn)
                self._properties.append(p)

    @property
    def num_props(self):
        return len(self._properties)

    @property
    def short_labels(self):
        return list(list(zip(*self._properties))[0])

    @property
    def long_labels(self):
        return list(list(zip(*self._properties))[1])

    # query external property
    def fetch_property(self, endpoint, prop, rdmol):
        name = endpoint['name']
        cache = self.api_cache.get(name)

        if cache is None or datetime.now() - cache['time'] > API_CACHE_TIME:
            url = endpoint['url']
            method = endpoint['method']
            data = endpoint['data']

            try:
                if data == 'smiles' and method == 'GET':
                    smiles = Chem.MolToSmiles(rdmol)
                    url = url.replace(':smiles', smiles)
                    json = requests.get(url).json()

                elif data == 'sdf' and method == 'POST':
                    Chem.SDWriter(self.temp_sdf.name).write(rdmol)
                    files = {'file': open(self.temp_sdf.name, 'rb')}
                    json = requests.post(url, files=files).json()

                else:
                    Logs.error(f'Unsupported request type: {method} {data} on {name}')
                    return None

            except:
                Logs.error(f'Failed to fetch {name}')
                return None

            data = {}
            for item, info in endpoint['properties'].items():
                value = json
                data[item] = None

                try:
                    for path in info['path'].replace('][', '].[').split('.'):
                        match = re.search(r'([^\[]+)?(?:\[(\d+)\])?', path)
                        (key, index) = match.groups()

                        if key is not None:
                            value = value.get(key)
                        if index is not None:
                            value = value[int(index)]

                    if type(value) not in [str, int, float, bool]:
                        raise TypeError

                    data[item] = value

                except:
                    Logs.error(f'Invalid path for {item} on {name}')

            cache = { 'time': datetime.now(), 'data': data }
            self.api_cache[name] = cache

        return cache['data'][prop]

    # adds complex.rdmol
    def prepare_complex(self, complex):
        # get only the current conformer
        m = next(complex.molecules)
        m.move_conformer(m.current_conformer, 0)
        m.set_conformer_count(1)

        complex.io.to_sdf(self.temp_sdf.name)
        complex.rdmol = Chem.SDMolSupplier(self.temp_sdf.name)[0]
        return complex.rdmol is not None

    # adds complex.properties
    def add_properties(self, complex):
        complex.properties = []
        for short_lbl, long_lbl, fmt, fn in self._properties:
            value = fn(complex.rdmol)
            value = fmt % value if value is not None else 'ERR'
            complex.properties.append((short_lbl, long_lbl, value))

    # adds complex.image
    def add_image(self, complex):
        complex.thumbnail = tempfile.NamedTemporaryFile(delete=False, suffix='.png', dir=self.temp_dir.name).name
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

        svg2png(bytestring=svg, write_to=complex.thumbnail, output_width=256, output_height=192)
        svg2png(bytestring=svg, write_to=complex.image, output_width=1024, output_height=768)

    def add_smiles(self, complex):
        complex.smiles = Chem.MolToSmiles(complex.rdmol)
