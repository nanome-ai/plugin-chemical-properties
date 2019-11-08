import nanome
from nanome.util import Logs
from nanome.util.enums import NotificationTypes

import os
import tempfile
import shutil
from cairosvg import svg2png
from datetime import datetime
from functools import partial

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from .PropertyCalculator import PropertyCalculator

BASE_DIR = os.path.dirname(__file__)
MENU_PATH = os.path.join(BASE_DIR, 'menu.json')
EDIT_ICON = os.path.join(BASE_DIR, 'icons', 'edit.png')
CHECK_ICON = os.path.join(BASE_DIR, 'icons', 'check.png')

# mol 2d image drawing options
Draw.DrawingOptions.atomLabelFontSize = 40
Draw.DrawingOptions.dotsPerAngstrom = 100
Draw.DrawingOptions.bondLineWidth = 8

class ChemicalProperties(nanome.PluginInstance):
    def start(self):
        self.calc = PropertyCalculator()

        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf", dir=self.temp_dir.name)
        self.temp_svg = tempfile.NamedTemporaryFile(delete=False, suffix=".svg", dir=self.temp_dir.name)

        self.displayed_complexes = []
        self.selected_complex = None
        self.selected_properties = [0, 1, 2, 3, 4, 5, 6, 7]
        self.snapshots = []
        self.snapshots_sort = [None, 0]
        self.snapshot_index = 0

        self.create_menu()
        self.on_run()

    def on_run(self):
        self.display_menu()

    def on_stop(self):
        # clean up temp files
        shutil.rmtree(self.temp_dir.name)

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
        self.pg_snapshots = menu.root.find_node("Snapshot View")

        # single view elements
        self.ln_results = menu.root.find_node("Results List")
        self.ln_error = menu.root.find_node("Results Error")

        self.ls_complexes = menu.root.find_node("Complex List").get_content()
        self.ls_properties = menu.root.find_node("Property List").get_content()
        self.ls_results = self.ln_results.get_content()

        self.btn_view_snapshots = menu.root.find_node("View Snapshots Button").get_content()
        self.btn_view_snapshots.register_pressed_callback(self.view_snapshots)

        self.btn_refresh = menu.root.find_node("Refresh Button").get_content()
        self.btn_refresh.register_pressed_callback(self.refresh_lists)

        self.btn_snapshot = menu.root.find_node("Snapshot Button").get_content()
        self.btn_snapshot.register_pressed_callback(self.snapshot_complex)

        # snapshots view elements
        self.ln_header = menu.root.find_node("Header")
        self.ln_properties_panel = menu.root.find_node("Properties Panel")
        self.ln_complex_info = menu.root.find_node("Complex Info Panel")

        self.ls_snapshots = menu.root.find_node("Snapshot List").get_content()
        self.ls_snapshots_properties = menu.root.find_node("Snapshot Property List").get_content()

        self.btn_exit_snapshots = menu.root.find_node("Exit Snapshots Button").get_content()
        self.btn_exit_snapshots.register_pressed_callback(self.view_single)

        self.btn_export = menu.root.find_node("Export Button").get_content()
        self.btn_export.register_pressed_callback(self.export_snapshots)

        self.btn_select_properties = menu.root.find_node("Select Properties Button").get_content()
        self.btn_select_properties.register_pressed_callback(self.toggle_properties_panel)

        # complex info panel elements
        self.ln_complex_name = menu.root.find_node("Complex Name")
        self.ln_complex_input = menu.root.find_node("Complex Name Input")
        self.ln_complex_image = menu.root.find_node("Complex Image")

        self.lbl_complex_name = self.ln_complex_name.get_content()
        self.lbl_complex_info = menu.root.find_node("Complex Info Label").get_content()

        self.inp_complex_name = self.ln_complex_input.get_content()

        self.btn_complex_edit = menu.root.find_node("Complex Edit Button").get_content()
        self.btn_complex_edit.register_pressed_callback(self.toggle_complex_edit)

        self.btn_complex_info_close = menu.root.find_node("Complex Close Button").get_content()
        self.btn_complex_info_close.register_pressed_callback(self.toggle_complex_panel)

        self.btn_complex_swap = menu.root.find_node("Complex Swap Button").get_content()
        self.btn_complex_swap.register_pressed_callback(partial(self.load_complex, True))

        self.btn_complex_load = menu.root.find_node("Complex Load Button").get_content()
        self.btn_complex_load.register_pressed_callback(partial(self.load_complex, False))

        self.btn_complex_delete = menu.root.find_node("Complex Delete Button").get_content()
        self.btn_complex_delete.register_pressed_callback(self.delete_complex)

        self.update_menu(self.menu)

    def display_menu(self):
        self.menu.enabled = True
        self.refresh_complexes()
        self.refresh_snapshot_button()
        self.populate_properties()
        self.view_single()

    def set_error(self, error=None):
        self.ln_results.enabled = error is None
        self.ln_error.enabled = error is not None

        if error is not None:
            lbl = self.ln_error.get_content()
            lbl.text_value = error

        self.update_menu(self.menu)

    def view_single(self, button=None):
        self.pg_single.enabled = True
        self.pg_snapshots.enabled = False
        self.update_menu(self.menu)

    def view_snapshots(self, button=None):
        self.pg_single.enabled = False
        self.pg_snapshots.enabled = True
        self.toggle_properties_panel(enable=False)
        self.refresh_snapshots()

    def refresh_snapshot_button(self):
        self.btn_snapshot.unusable = not self.selected_complex
        self.update_content(self.btn_snapshot)

    def refresh_lists(self, button=None):
        self.refresh_complexes()
        self.refresh_results()
        self.refresh_snapshot_button()

    def refresh_complexes(self, button=None):
        def select_complex(button):
            self.selected_complex = None
            self.refresh_snapshot_button()

            self.selected_complex = button.complex
            for item in self.ls_complexes.items:
                btn = item.get_content()
                if btn.selected:
                    btn.selected = False
                    self.update_content(btn)
                    break
            button.selected = True
            self.update_content(button)

            self.ls_results.items.clear()
            self.update_content(self.ls_results)
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
                btn.complex = complex
                btn.register_pressed_callback(select_complex)

                if self.selected_complex:
                    btn.selected = complex.index == self.selected_complex.index
                if btn.selected:
                    self.selected_complex = complex

                self.ls_complexes.items.append(item)
            self.update_content(self.ls_complexes)

        self.request_complex_list(display_complexes)

    def refresh_properties(self):
        self.refresh_property_list(self.ls_properties)
        self.refresh_property_list(self.ls_snapshots_properties)

    def refresh_property_list(self, property_list):
        for i, item in enumerate(property_list.items):
            btn = item.get_content()
            btn.selected = i in self.selected_properties
        self.update_content(property_list)

    def populate_properties(self):
        self.populate_property_list(self.ls_properties)
        self.populate_property_list(self.ls_snapshots_properties)

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

            if property_list == self.ls_snapshots_properties:
                self.refresh_snapshots()

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
        self.set_error(None)

        if self.selected_complex is None:
            self.ls_results.items.clear()
            self.update_content(self.ls_results)
            return

        def full_complexes_received(complexes):
            properties = self.calculate_properties(complexes[0], self.selected_properties)
            if properties is None:
                self.selected_complex = None
                self.refresh_snapshot_button()
                return

            self.refresh_snapshot_button()

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
        try:
            # get only the current conformer
            m = next(complex.molecules)
            m.move_conformer(m.current_conformer, 0)
            m.set_conformer_count(1)

            complex.io.to_sdf(self.temp_sdf.name)
            complex.mol = Chem.SDMolSupplier(self.temp_sdf.name)[0]
            if complex.mol is None:
                raise

            return [self.calc.calc_property(complex.mol, i) for i in property_indices]
        except:
            error = "rdkit was unable to parse the complex"
            self.send_notification(NotificationTypes.error, error)
            self.set_error(error)
            return None

    def snapshot_complex(self, button=None):
        def generate_image(complex):
            complex.image = tempfile.NamedTemporaryFile(delete=False, suffix=".png", dir=self.temp_dir.name)

            mol = complex.mol
            Chem.AssignStereochemistryFrom3D(mol)
            AllChem.Compute2DCoords(mol)
            mol = Draw.rdMolDraw2D.PrepareMolForDrawing(mol)

            drawer = Draw.rdMolDraw2D.MolDraw2DSVG(256, 192)
            drawer.drawOptions().additionalAtomLabelPadding = 0.3
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            svg = svg.replace("stroke-linecap:butt", "stroke-linecap:round")

            svg2png(bytestring=svg, write_to=complex.image.name, output_width=1024, output_height=768)

        def full_complexes_received(complexes):
            complex = complexes[0]
            complex.snapshot_name = str(self.snapshot_index)
            complex.timestamp = datetime.now()
            complex.properties = self.calculate_properties(complex, range(self.calc.num_props))
            generate_image(complex)

            self.snapshots.append(complex)
            self.snapshot_index += 1
            self.send_notification(NotificationTypes.success, "snapshot %s created" % complex.snapshot_name)

        self.request_complexes([self.selected_complex.index], full_complexes_received)

    def toggle_properties_panel(self, button=None, enable=None):
        show_sidepanel = enable if enable is not None else not self.ln_properties_panel.enabled
        self.ln_properties_panel.enabled = show_sidepanel

        text = "done" if show_sidepanel else "select properties"
        self.btn_select_properties.set_all_text(text)

        if enable is None:
            self.update_menu(self.menu)

    def toggle_complex_panel(self, button=None):
        self.ln_complex_info.enabled = not self.ln_complex_info.enabled

        if self.ln_complex_info.enabled:
            complex = button.complex
            self.ln_complex_info.complex = complex
            self.lbl_complex_name.text_value = complex.snapshot_name
            self.lbl_complex_info.text_value = complex.timestamp.strftime("%Y-%m-%d %H:%M:%S")

            img = self.ln_complex_image.add_new_image(complex.image.name)
            img.scaling_option = nanome.util.enums.ScalingOptions.fit

            self.btn_complex_swap.complex = complex
            self.btn_complex_load.complex = complex
            self.btn_complex_delete.complex = complex

            self.toggle_complex_edit(edit=False)

        self.update_menu(self.menu)

    def toggle_complex_edit(self, button=None, edit=None):
        is_editing = edit if edit is not None else self.ln_complex_name.enabled
        self.ln_complex_name.enabled = not is_editing
        self.ln_complex_input.enabled = is_editing

        icon = CHECK_ICON if is_editing else EDIT_ICON
        self.btn_complex_edit.set_all_icon(icon)

        if is_editing:
            complex_name = self.ln_complex_info.complex.snapshot_name
            self.inp_complex_name.input_text = complex_name
        elif edit is None:
            complex_name = self.inp_complex_name.input_text
            self.ln_complex_info.complex.snapshot_name = complex_name
            self.lbl_complex_name.text_value = complex_name
            self.refresh_snapshots_values()

        if edit is None:
            self.update_menu(self.menu)

    def set_snapshots_sort(self, button):
        prop_index = button.prop_index

        if self.snapshots_sort[0] == prop_index: # if already sorting by prop
            if self.snapshots_sort[1] == -1: # if reverse, set to no sort
                self.snapshots_sort = [None, 0]
                prop_index = None
            else: # reverse sort order
                self.snapshots_sort[1] = -self.snapshots_sort[1]
        else: # set to new prop sort
            self.snapshots_sort = [prop_index, 1]

        for heading in self.ln_header.get_children():
            btn = heading.get_content()
            old_selected = btn.selected
            btn.selected = btn.prop_index == prop_index

            if old_selected != btn.selected:
                self.update_content(btn)

        self.refresh_snapshots_values()

    def refresh_snapshots(self, button=None):
        self.refresh_snapshots_header()
        self.refresh_snapshots_values()

    def refresh_snapshots_header(self, button=None):
        self.ln_header.clear_children()

        heading = self.pfb_heading.clone()
        btn = heading.get_content()
        btn.set_all_text("ID")
        btn.prop_index = -1
        btn.selected = self.snapshots_sort[0] == -1 and self.snapshots_sort[1] != 0
        btn.register_pressed_callback(self.set_snapshots_sort)
        self.ln_header.add_child(heading)

        for prop_index in self.selected_properties:
            heading = self.pfb_heading.clone()
            btn = heading.get_content()
            btn.set_all_text(self.calc.short_labels[prop_index])
            btn.prop_index = prop_index
            btn.selected = self.snapshots_sort[0] == prop_index
            btn.register_pressed_callback(self.set_snapshots_sort)
            self.ln_header.add_child(heading)

    def refresh_snapshots_values(self, button=None):
        sorted_complexes = self.snapshots
        if self.snapshots_sort[1] != 0:
            prop_index = self.snapshots_sort[0]
            reverse = self.snapshots_sort[1] == -1

            if self.snapshots_sort[0] == -1: # sort by name first, then timestamp
                sort_fn = lambda complex: (complex.snapshot_name, complex.timestamp)
            else: # sort by property value
                sort_fn = lambda complex: float(complex.properties[prop_index][2])

            sorted_complexes = sorted(sorted_complexes, key=sort_fn, reverse=reverse)

        self.btn_export.unusable = not self.snapshots

        self.ls_snapshots.items.clear()
        for complex in sorted_complexes:
            ln = nanome.ui.LayoutNode()
            ln.layout_orientation = nanome.util.enums.LayoutTypes.horizontal

            ln_lbl = self.pfb_heading.clone()
            btn = ln_lbl.get_content()
            btn.set_all_text(complex.snapshot_name)
            btn.complex = complex
            btn.register_pressed_callback(self.toggle_complex_panel)
            ln.add_child(ln_lbl)

            for prop_index in self.selected_properties:
                ln_lbl = self.pfb_value.clone()
                lbl = ln_lbl.get_content()
                lbl.text_value = complex.properties[prop_index][2]
                ln.add_child(ln_lbl)

            self.ls_snapshots.items.append(ln)

        self.update_menu(self.menu)

    def load_complex(self, swap, button):
        complex = button.complex

        # hack to fix bonds not being loaded
        old_complex = complex
        complex.io.to_sdf(self.temp_sdf.name)
        complex = nanome.structure.Complex.io.from_sdf(path=self.temp_sdf.name)
        complex.index = old_complex.index
        complex.snapshot_name = old_complex.snapshot_name
        complex.full_name = old_complex.full_name
        complex.position = old_complex.position
        complex.rotation = old_complex.rotation

        if not swap:
            complex.index = -1
            complex.full_name = complex.snapshot_name
        self.update_structures_deep([complex])

    def delete_complex(self, button):
        for complex in self.snapshots:
            if complex == button.complex:
                self.snapshots.remove(complex)
                self.refresh_snapshot_button()
                break

        self.refresh_snapshots_values()
        self.toggle_complex_panel()

    def export_snapshots(self, button=None):
        file = nanome.util.FileSaveData()
        filename = datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
        file.path = 'snapshots\\' + filename
        file.write_text(','.join(['ID', 'SMILES'] + self.calc.short_labels) + '\n')

        for complex in self.snapshots:
            smiles = Chem.MolToSmiles(complex.mol)
            values = list(list(zip(*complex.properties))[2])
            file.write_text(','.join([complex.snapshot_name, smiles] + values) + '\n')

        def on_save_files_result(result_list):
            result = result_list[0]

            if result.error_code == nanome.util.FileErrorCode.no_error:
                self.send_notification(NotificationTypes.success, "saved to Documents\\nanome\\snapshots")
            else:
                self.send_notification(NotificationTypes.error, "error exporting csv")

        self.save_files([file], on_save_files_result)

def main():
    plugin = nanome.Plugin("Chemical Properties", "Calculates and displays different properties of chemicals using the RDKit Python library", "", False)
    plugin.set_plugin_class(ChemicalProperties)
    plugin.run('127.0.0.1', 8888)

if __name__ == "__main__":
    main()
