import nanome
from nanome.util.enums import NotificationTypes

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
            ("logP", "Lipophilicity (logP)", "%.3f", lambda mol: mDesc.CalcCrippenDescriptors(mol)[0]),
            ("TPSA", "Total Polar Surface Area", "%.3f", mDesc.CalcTPSA),
            ("ESOL", "Estimated Solubility", "%.3f", self.esol.calc_esol),
            ("HBA", "# H-Bond Acceptors", "%d", mDesc.CalcNumHBA),
            ("HBD", "# H-Bond Donors", "%d", mDesc.CalcNumHBD),
            ("RB", "# Rotatable Bonds", "%d", mDesc.CalcNumRotatableBonds),
            ("AR", "# Aromatic Rings", "%d", mDesc.CalcNumAromaticRings),
            ("EMW", "Exact Molecular Weight", "%.3f", mDesc.CalcExactMolWt)
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

class ChemicalProperties(nanome.PluginInstance):
    def start(self):
        self.calc = PropertyCalculator()

        self.temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self.displayed_complexes = []
        self.selected_complex = None
        self.selected_properties = [0, 1, 2, 3, 4, 5, 6, 7]
        self.table_complexes = []
        self.table_sort = [None, 0]

        self.create_menu()
        self.on_run()

    def on_run(self):
        self.display_menu()

    def on_complex_added(self):
        self.refresh_complexes()

    def on_complex_removed(self):
        self.refresh_complexes()

    def create_menu(self):
        self.menu = nanome.ui.Menu.io.from_json(MENU_PATH)
        menu = self.menu

        self.pfb_heading = menu.root.find_node("Prefab Heading")
        self.pfb_value = menu.root.find_node("Prefab Value")
        self.pfb_complex = menu.root.find_node("Prefab Complex")
        self.pfb_property = menu.root.find_node("Prefab Property")
        self.pfb_result = menu.root.find_node("Prefab Result")

        self.pg_single = menu.root.find_node("Single View")
        self.pg_table = menu.root.find_node("Table View")

        # single view elements
        self.ls_complexes = menu.root.find_node("Complex List").get_content()
        self.ls_properties = menu.root.find_node("Property List").get_content()
        self.ls_results = menu.root.find_node("Results List").get_content()

        self.btn_view_table = menu.root.find_node("View Table Button").get_content()
        self.btn_view_table.register_pressed_callback(self.view_table)

        self.btn_refresh = menu.root.find_node("Refresh Button").get_content()
        self.btn_refresh.register_pressed_callback(self.refresh_lists)

        self.btn_add_to_table = menu.root.find_node("Add to Table Button").get_content()
        self.btn_add_to_table.register_pressed_callback(self.toggle_complex_in_table)

        # table view elements
        self.ln_header = menu.root.find_node("Header")
        self.ln_sidepanel = menu.root.find_node("Table Side Panel")
        self.ls_table = menu.root.find_node("Table List").get_content()
        self.ls_table_properties = menu.root.find_node("Table Property List").get_content()

        self.btn_exit_table = menu.root.find_node("Exit Table Button").get_content()
        self.btn_exit_table.register_pressed_callback(self.view_single)

        self.btn_select_properties = menu.root.find_node("Select Properties Button").get_content()
        self.btn_select_properties.register_pressed_callback(self.toggle_table_sidepanel)

        self.display_menu()

    def display_menu(self):
        self.menu.enabled = True
        self.refresh_complexes()
        self.refresh_add_button()
        self.populate_properties()
        self.update_menu(self.menu)

    def view_single(self, button=None):
        self.pg_single.enabled = True
        self.pg_table.enabled = False
        self.update_menu(self.menu)

    def view_table(self, button=None):
        self.pg_single.enabled = False
        self.pg_table.enabled = True
        self.toggle_table_sidepanel(enable=False)
        self.refresh_table()

    def refresh_add_button(self):
        self.btn_add_to_table.unusable = not self.selected_complex

        complex_indices = [complex.index for complex in self.table_complexes]
        complex_in_table = self.selected_complex and self.selected_complex.index in complex_indices
        text = "remove from table" if complex_in_table else "add to table"
        self.btn_add_to_table.set_all_text(text)

        self.update_content(self.btn_add_to_table)

    def refresh_lists(self, button=None):
        self.refresh_complexes()
        self.refresh_results()

    def refresh_complexes(self, button=None):
        def select_complex(button):
            self.selected_complex = None
            self.refresh_add_button()

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
            self.displayed_complexes = [complex.index for complex in complexes]

            if self.selected_complex and self.selected_complex.index not in self.displayed_complexes:
                self.selected_complex = None
                self.refresh_results()

            self.ls_complexes.items.clear()
            for complex in complexes:
                item = self.pfb_complex.clone()
                btn = item.get_content()
                btn.set_all_text(complex.full_name)
                if self.selected_complex:
                    btn.selected = complex.index == self.selected_complex.index
                btn.complex = complex
                btn.register_pressed_callback(select_complex)
                self.ls_complexes.items.append(item)
            self.update_content(self.ls_complexes)

        self.request_complex_list(display_complexes)

    def refresh_properties(self):
        self.refresh_property_list(self.ls_properties)
        self.refresh_property_list(self.ls_table_properties)

    def refresh_property_list(self, property_list):
        for i, item in enumerate(property_list.items):
            btn = item.get_content()
            btn.selected = i in self.selected_properties
        self.update_content(property_list)

    def populate_properties(self):
        self.populate_property_list(self.ls_properties)
        self.populate_property_list(self.ls_table_properties)

    def populate_property_list(self, property_list):
        def property_pressed(button):
            button.selected = not button.selected
            self.update_content(button)

            self.selected_properties.clear()
            for item in property_list.items:
                btn = item.get_content()
                if btn.selected:
                    self.selected_properties.append(btn.property_index)
            self.refresh_results()
            self.refresh_properties()

            if property_list == self.ls_table_properties:
                self.refresh_table()

        property_list.items.clear()
        for i, label in enumerate(self.calc.long_labels):
            item = self.pfb_property.clone()
            btn = item.get_content()
            btn.property_index = i
            btn.selected = i in self.selected_properties
            btn.register_pressed_callback(property_pressed)
            lbl = item.find_node("Label").get_content()
            lbl.text_value = label
            property_list.items.append(item)
        self.update_content(property_list)

    def refresh_results(self, button=None):
        if self.selected_complex is None:
            self.ls_results.items.clear()
            self.update_content(self.ls_results)
            return

        def full_complexes_received(complexes):
            properties = self.calculate_properties(complexes[0], self.selected_properties)
            if properties is None:
                self.selected_complex = None
                self.refresh_add_button()
                return
            else:
                self.refresh_add_button()

            self.ls_results.items.clear()
            for prop in properties:
                item = self.pfb_result.clone()
                name = item.find_node("Name").get_content()
                name.text_value = prop[0]
                value = item.find_node("Value").get_content()
                value.text_value = prop[2]
                self.ls_results.items.append(item)
            self.update_content(self.ls_results)

        self.request_complexes([self.selected_complex.index], full_complexes_received)

    def calculate_properties(self, complex, property_indices):
        complex.io.to_sdf(self.temp_sdf.name)

        try:
            mol = Chem.SDMolSupplier(self.temp_sdf.name)[0]
            return [self.calc.calc_property(mol, i) for i in property_indices]
        except:
            self.send_notification(NotificationTypes.error, "rdkit was unable to parse the molecule")

    def toggle_complex_in_table(self, button=None):
        # if in table, remove it
        for complex in self.table_complexes:
            if complex.index == self.selected_complex.index:
                self.table_complexes.remove(complex)
                self.refresh_add_button()
                return

        def full_complexes_received(complexes):
            complex = complexes[0]
            complex.properties = self.calculate_properties(complex, range(self.calc.num_props))
            self.table_complexes.append(complex)
            self.refresh_add_button()

        self.request_complexes([self.selected_complex.index], full_complexes_received)

    def toggle_table_sidepanel(self, button=None, enable=None):
        show_sidepanel = enable if enable is not None else not self.ln_sidepanel.enabled
        self.ln_sidepanel.enabled = show_sidepanel
        text = "done" if show_sidepanel else "select properties"
        self.btn_select_properties.set_all_text(text)
        if enable is None:
            self.update_menu(self.menu)

    def set_table_sort(self, button):
        prop_index = button.prop_index

        if self.table_sort[0] == prop_index: # if already sorting by prop
            if self.table_sort[1] == -1: # if reverse, set to no sort
                self.table_sort = [None, 0]
                prop_index = None
            else: # reverse sort order
                self.table_sort[1] = -self.table_sort[1]
        else: # set to new prop sort
            self.table_sort = [prop_index, 1]

        for heading in self.ln_header.get_children():
            btn = heading.get_content()
            old_selected = btn.selected
            btn.selected = btn.prop_index == prop_index

            if old_selected != btn.selected:
                self.update_content(btn)

        self.refresh_table_values()

    def refresh_table(self, button=None):
        self.refresh_table_header()
        self.refresh_table_values()

    def refresh_table_header(self, button=None):
        self.ln_header.clear_children()

        heading = self.pfb_heading.clone()
        btn = heading.get_content()
        btn.set_all_text("complex")
        btn.prop_index = -1
        btn.selected = self.table_sort[0] == -1 and self.table_sort[1] != 0
        btn.register_pressed_callback(self.set_table_sort)
        self.ln_header.add_child(heading)

        for prop_index in self.selected_properties:
            heading = self.pfb_heading.clone()
            btn = heading.get_content()
            btn.set_all_text(self.calc.short_labels[prop_index])
            btn.prop_index = prop_index
            btn.selected = self.table_sort[0] == prop_index
            btn.register_pressed_callback(self.set_table_sort)
            self.ln_header.add_child(heading)

    def refresh_table_values(self, button=None):
        sorted_complexes = self.table_complexes
        if self.table_sort[1] != 0:
            prop_index = self.table_sort[0]
            reverse = self.table_sort[1] == -1

            if self.table_sort[0] == -1: # sort by name
                sort_fn = lambda complex: complex.full_name
            else: # sort by property value
                sort_fn = lambda complex: float(complex.properties[prop_index][2])

            sorted_complexes = sorted(sorted_complexes, key=sort_fn, reverse=reverse)

        self.ls_table.items.clear()
        for complex in sorted_complexes:
            ln = nanome.ui.LayoutNode()
            ln.layout_orientation = nanome.util.enums.LayoutTypes.horizontal

            ln_lbl = self.pfb_value.clone()
            lbl = ln_lbl.get_content()
            lbl.text_value = complex.full_name
            ln.add_child(ln_lbl)

            for prop_index in self.selected_properties:
                ln_lbl = self.pfb_value.clone()
                lbl = ln_lbl.get_content()
                lbl.text_value = complex.properties[prop_index][2]
                ln.add_child(ln_lbl)

            self.ls_table.items.append(ln)

        self.update_menu(self.menu)

def main():
    plugin = nanome.Plugin("Chemical Properties", "Calculates and displays different properties of chemicals using the RDKit Python library", "", False)
    plugin.set_plugin_class(ChemicalProperties)
    plugin.run('127.0.0.1', 8888)

if __name__ == "__main__":
    main()
