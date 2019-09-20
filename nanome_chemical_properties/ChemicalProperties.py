import nanome
from nanome.util.enums import HorizAlignOptions

import os
import tempfile

from rdkit import Chem
import rdkit.Chem.Descriptors as Desc
import rdkit.Chem.rdMolDescriptors as mDesc
from .ESOLCalculator import ESOLCalculator

MENU_PATH = os.path.join(os.path.dirname(__file__), 'menu.json')

class PropertyCalculator():
    def __init__(self):
        self.esol = ESOLCalculator()
        self._properties = [
            ("MW", "Molecular Weight", "%.3f", Desc.MolWt),
            ("EMW", "Exact Molecular Weight", "%.3f", mDesc.CalcExactMolWt),
            ("logP", "Lipophilicity (logP)", "%.3f", lambda mol: mDesc.CalcCrippenDescriptors(mol)[0]),
            ("TPSA", "Total Polar Surface Area", "%.3f", mDesc.CalcTPSA),
            ("ESOL", "Estimated Solubility", "%.3f", self.esol.calc_esol),
            ("# HBA", "# H-Bond Donors", "%d", mDesc.CalcNumHBA),
            ("# HBD", "# H-Bond Acceptors", "%d", mDesc.CalcNumHBD),
            ("# RB", "# Rotatable Bonds", "%d", mDesc.CalcNumRotatableBonds),
            ("# AR", "# Aromatic Rings", "%d", mDesc.CalcNumAromaticRings)
        ]

    def calc_property(self, mol, index):
        (short_lbl, long_lbl, fmt, fn) = self._properties[index]
        return (short_lbl, long_lbl, fmt % fn(mol))

    @property
    def short_labels(self):
        return list(list(zip(*self._properties))[0])

    @property
    def long_labels(self):
        return list(list(zip(*self._properties))[1])

class ChemicalProperties(nanome.PluginInstance):
    def start(self):
        self.calc = PropertyCalculator()

        self.temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self.selected_complex = None
        self.selected_properties = []

        self.create_menu()
        self.on_run()

    def on_run(self):
        # self.do_the_thing()
        self.display_menu()

    def do_the_thing(self):
        def got_deep(complexes):
            complex = complexes[0]
            complex.io.to_sdf(self.temp_sdf.name)
            complex = nanome.api.structure.Complex.io.from_sdf(path=self.temp_sdf.name)
            self.add_to_workspace([complex])

        def got_shallow(complexes):
            self.request_complexes([complexes[0].index], got_deep)

        self.request_complex_list(got_shallow)

    def create_menu(self):
        self.menu = nanome.ui.Menu.io.from_json(MENU_PATH)
        menu = self.menu

        self.pfb_complex = menu.root.find_node("Prefab Complex")
        self.pfb_property = menu.root.find_node("Prefab Property")
        self.pfb_result = menu.root.find_node("Prefab Result")

        self.ls_complexes = menu.root.find_node("Complex List").get_content()
        self.ls_properties = menu.root.find_node("Property List").get_content()
        self.ls_results = menu.root.find_node("Results List").get_content()

        self.btn_refresh = menu.root.find_node("Refresh Button").get_content()
        self.btn_refresh.register_pressed_callback(self.refresh_lists)

        self.display_menu()

    def display_menu(self):
        self.menu.enabled = True
        self.refresh_complexes()
        self.refresh_properties()
        self.update_menu(self.menu)

    def refresh_lists(self, button=None):
        self.refresh_complexes()
        self.refresh_results()

    def refresh_complexes(self, button=None):
        def select_complex(button):
            self.selected_complex = button.complex
            for item in self.ls_complexes.items:
                btn = item.get_content()
                if btn.selected:
                    btn.selected = False
                    self.update_content(btn)
                    break
            button.selected = True
            self.update_content(button)
            self.refresh_results()

        def display_complexes(complexes):
            self.ls_complexes.items.clear()
            for complex in complexes:
                item = self.pfb_complex.clone()
                btn = item.get_content()
                btn.set_all_text(complex.name)
                if self.selected_complex:
                    btn.selected = complex.index == self.selected_complex.index
                btn.complex = complex
                btn.register_pressed_callback(select_complex)
                self.ls_complexes.items.append(item)
            self.update_content(self.ls_complexes)

        self.request_complex_list(display_complexes)

    def refresh_properties(self, button=None):
        def property_pressed(button):
            button.selected = not button.selected
            self.update_content(button)

            self.selected_properties.clear()
            for item in self.ls_properties.items:
                btn = item.get_content()
                if btn.selected:
                    self.selected_properties.append(btn.property_index)
            self.refresh_results()

        self.ls_properties.items.clear()
        for i, label in enumerate(self.calc.long_labels):
            item = self.pfb_property.clone()
            btn = item.get_content()
            btn.property_index = i
            btn.selected = i in self.selected_properties
            btn.register_pressed_callback(property_pressed)
            lbl = item.find_node("Label").get_content()
            lbl.text_value = label
            self.ls_properties.items.append(item)
        self.update_content(self.ls_properties)

    def refresh_results(self, button=None):
        if self.selected_complex is None:
            return

        def full_complexes_received(complexes):
            complex = complexes[0]
            complex.io.to_sdf(self.temp_sdf.name)

            try:
                mol = Chem.SDMolSupplier(self.temp_sdf.name)[0]
            except:
                self.send_notification(NotificationTypes.error, "rdkit was unable to parse the molecule")
                return

            self.ls_results.items.clear()
            for i in self.selected_properties:
                prop = self.calc.calc_property(mol, i)
                item = self.pfb_result.clone()
                name = item.find_node("Name").get_content()
                name.text_value = prop[0]
                value = item.find_node("Value").get_content()
                value.text_value = prop[2]
                self.ls_results.items.append(item)
            self.update_content(self.ls_results)

        self.request_complexes([self.selected_complex.index], full_complexes_received)

def main():
    plugin = nanome.Plugin("Chemical Properties", "Calculates and displays different properties of chemicals", "", False)
    plugin.set_plugin_class(ChemicalProperties)
    plugin.run('192.168.1.49', 8888)

if __name__ == "__main__":
    main()
