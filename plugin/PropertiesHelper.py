from nanome.util import Color, Logs

from rdkit import Chem
from rdkit.Chem import Draw
import rdkit.Chem.Descriptors as Desc
import rdkit.Chem.rdMolDescriptors as mDesc
from .ESOLCalculator import ESOLCalculator

import json
import os
import re
import requests
import tempfile
from cairosvg import svg2png
from collections import namedtuple
from datetime import datetime, timedelta
from functools import partial
from math import inf
from urllib.parse import quote

# mol 2d image drawing options
Draw.DrawingOptions.atomLabelFontSize = 40
Draw.DrawingOptions.dotsPerAngstrom = 100
Draw.DrawingOptions.bondLineWidth = 8

API_CACHE_TIME = timedelta(seconds=1)
API_SETTINGS = os.path.join(os.path.dirname(__file__), '..', 'config.json')

Property = namedtuple(
    'Property',
    ['name', 'description', 'format', 'fn', 'color_fn'],
    defaults=['', '', '%s', lambda x: None, lambda x: Color.White()])

PropertyValue = namedtuple(
    'PropertyValue',
    ['name', 'description', 'value', 'color'],
    defaults=['', '', '', Color.White()])

def within(min=-inf, max=inf):
    return lambda x: Color.White() if min <= x <= max else Color.Red()

def gradient(colors):
    for item in colors:
        item[1] = Color.from_int(int(item[1], 16))

    def color(x):
        for (x1, c1), (x2, c2) in zip(colors[:-1], colors[1:]):
            if x <= x1: return c1
            if x > x2: continue
            t = (x - x1) / (x2 - x1)
            return c1 * (1 - t) + c2 * t
        return colors[-1][1]

    return color

class PropertiesHelper:
    def __init__(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)

        self.api_cache = {}
        self.esol = ESOLCalculator()
        self._properties = [
            Property('MW', 'Molecular Weight', '%.3f', Desc.MolWt, within(max=500)),
            Property('logP', 'Lipophilicity (logP)', '%.3f', lambda mol: mDesc.CalcCrippenDescriptors(mol)[0], within(max=5)),
            Property('TPSA', 'Total Polar Surface Area', '%.3f', mDesc.CalcTPSA),
            Property('ESOL', 'Estimated Solubility', '%.3f', self.esol.calc_esol),
            Property('HBA', '# H-Bond Acceptors', '%d', mDesc.CalcNumHBA, within(max=10)),
            Property('HBD', '# H-Bond Donors', '%d', mDesc.CalcNumHBD, within(max=5)),
            Property('RB', '# Rotatable Bonds', '%d', mDesc.CalcNumRotatableBonds),
            Property('AR', '# Aromatic Rings', '%d', mDesc.CalcNumAromaticRings)
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

                color_fn = lambda x: Color.White()
                color_settings = info.get('color')

                if color_settings:
                    if color_settings['type'] == 'within':
                        color_fn = within(**color_settings['args'])
                    elif color_settings['type'] == 'gradient':
                        color_fn = gradient(**color_settings['args'])

                p = Property(prop, info['description'], info['format'], fn, color_fn)
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
                    url = url.replace(':smiles', quote(smiles))
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
        for short_lbl, long_lbl, fmt, fn, color_fn in self._properties:
            value = fn(complex.rdmol)
            color = Color.Red() if value is None else color_fn(value)
            value = 'ERR' if  value is None else fmt % value
            prop = PropertyValue(short_lbl, long_lbl, value, color)
            complex.properties.append(prop)

    # adds complex.image
    def add_image(self, complex):
        complex.thumbnail = tempfile.NamedTemporaryFile(delete=False, suffix='.png', dir=self.temp_dir.name).name
        complex.image = tempfile.NamedTemporaryFile(delete=False, suffix='.png', dir=self.temp_dir.name).name

        mol = complex.rdmol
        Chem.AssignStereochemistryFrom3D(mol)
        Chem.rdCoordGen.AddCoords(mol)
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
