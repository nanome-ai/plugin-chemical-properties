import nanome
from nanome.util import async_callback
from nanome.util.enums import HorizAlignOptions, NotificationTypes, ScalingOptions

import os
from datetime import datetime

MAX_ATOM_COUNT = 200

BASE_DIR = os.path.join(os.path.dirname(__file__))
MENU_PATH = os.path.join(BASE_DIR, 'json/main.json')
IMG_FRAME = os.path.join(BASE_DIR, '../icons/frame.png')
IMG_SELECTION = os.path.join(BASE_DIR, '../icons/selection.png')

NO_CALLBACK = lambda *_: None

class MainMenu:
    def __init__(self, plugin: nanome.AsyncPluginInstance):
        self.plugin = plugin
        self.selected_complex_index = None
        self.selected_ligand_index = None
        self.selected_complex = None
        self.selected_ligand = None
        self.create_menu()

    def create_menu(self):
        self.menu = nanome.ui.Menu.io.from_json(MENU_PATH)
        root = self.menu.root

        self.pfb_complex = root.find_node('Prefab Complex')
        self.pfb_result = root.find_node('Prefab Result')

        self.ln_panel_left = root.find_node('Panel Left')
        self.ln_panel_right = root.find_node('Panel Right')

        self.ln_preview = root.find_node('Complex Preview')
        self.ln_ligand_list = root.find_node('Ligand List Panel')
        self.ln_message = root.find_node('Preview Message')
        self.ln_frame = root.find_node('Preview Frame')
        self.ln_image = root.find_node('Preview Image')

        self.ln_results = root.find_node('Results')
        self.ln_no_selection = root.find_node('No Selection')

        self.btn_refresh = root.find_node('Button Refresh').get_content()
        self.btn_snapshots = root.find_node('Button Snapshots').get_content()
        self.btn_snapshot = root.find_node('Button Snapshot').get_content()

        self.lst_complexes = root.find_node('Structure List').get_content()
        self.lst_ligands = root.find_node('Ligand List').get_content()
        self.lst_results = root.find_node('Results List').get_content()

        self.lbl_message = self.ln_message.get_content()
        self.lbl_complex = root.find_node('Complex Title').get_content()

        # add images, callbacks, etc
        img = self.ln_frame.add_new_image(IMG_FRAME)
        img.scaling_option = ScalingOptions.fit
        img = root.find_node('Select Image').add_new_image(IMG_SELECTION)
        img.scaling_option = ScalingOptions.fit

        self.btn_refresh.register_pressed_callback(self.refresh_lists)
        self.btn_snapshot.register_pressed_callback(self.snapshot_complex)

        show_snapshots = lambda b: self.plugin.menu_snapshots.show_menu()
        self.btn_snapshots.register_pressed_callback(show_snapshots)

    def show_menu(self):
        self.menu.enabled = True
        self.refresh_lists()
        self.refresh_snapshot_button()
        self.plugin.update_menu(self.menu)

    def reset_selection(self):
        self.selected_complex = None
        self.selected_complex_index = None
        self.selected_ligand = None
        self.selected_ligand_index = None

    def refresh_lists(self, button=None):
        self.refresh_complexes()
        self.refresh_selection()

    def hide_ligand_list(self):
        self.lst_ligands.items.clear()
        self.ln_ligand_list.enabled = False
        self.plugin.update_node(self.ln_panel_left)

    @async_callback
    async def refresh_complexes(self):
        complexes = await self.plugin.request_complex_list()
        complex_indexes = [complex.index for complex in complexes]

        if self.selected_complex_index not in complex_indexes:
            self.reset_selection()
            self.refresh_snapshot_button()
            self.refresh_results()

        self.lst_complexes.items.clear()
        for complex in complexes:
            item = self.pfb_complex.clone()
            btn = item.get_content()
            btn.text.value.set_all(complex.full_name)
            btn.complex_index = complex.index
            btn.register_pressed_callback(self.select_complex)
            if self.selected_complex_index:
                btn.selected = complex.index == self.selected_complex_index
            self.lst_complexes.items.append(item)

        if not complexes:
            self.lst_complexes.items.append(nanome.ui.LayoutNode())
            ln = nanome.ui.LayoutNode()
            lbl = ln.add_new_label('no structures')
            lbl.text_max_size = 0.4
            lbl.text_horizontal_align = HorizAlignOptions.Middle
            self.lst_complexes.items.append(ln)

        self.plugin.update_content(self.lst_complexes)

    @async_callback
    async def refresh_ligands(self):
        if not self.selected_complex:
            self.hide_ligand_list()
            return

        molecule = next(
            mol for i, mol in enumerate(self.selected_complex.molecules)
            if i == self.selected_complex.current_frame)
        ligands = await molecule.get_ligands()

        if not ligands:
            self.hide_ligand_list()
            self.compute_results()
            return

        ligand_complexes = []
        for ligand in ligands:
            ligand_complex = nanome.structure.Complex()
            ligand_complex.name = ligand.name
            ligand_molecule = nanome.structure.Molecule()
            ligand_chain = nanome.structure.Chain()
            ligand_complex.add_molecule(ligand_molecule)
            ligand_molecule.add_chain(ligand_chain)
            for residue in ligand.residues:
                ligand_chain.add_residue(residue)
            ligand_complexes.append(ligand_complex)

        if len(ligand_complexes) == 1:
            mol_num_atoms = sum(1 for _ in molecule.atoms)
            lig_num_atoms = sum(1 for _ in ligand_complexes[0].atoms)
            if mol_num_atoms == lig_num_atoms:
                self.hide_ligand_list()
                self.compute_results()
                return

        self.lst_ligands.items.clear()
        for i, ligand in enumerate([None, *ligand_complexes]):
            item = self.pfb_complex.clone()
            btn = item.get_content()
            btn.text.value.set_all(ligand.name if ligand else 'none')
            btn.selected = not ligand
            btn.ligand = ligand
            btn.ligand_index = i
            btn.register_pressed_callback(self.select_ligand)
            self.lst_ligands.items.append(item)

        self.ln_ligand_list.enabled = True
        self.plugin.update_node(self.ln_panel_left)

        select_index = self.selected_ligand_index or 1
        self.select_ligand(self.lst_ligands.items[select_index].get_content())

    @async_callback
    async def refresh_selection(self, complex=None):
        if not self.selected_complex_index:
            return

        self.lbl_complex.text_value = 'loading...'
        self.lst_results.items.clear()
        self.plugin.update_node(self.ln_results)
        self.update_preview(text='loading...')

        [complex] = await self.plugin.request_complexes([self.selected_complex_index])
        self.selected_complex = complex
        self.refresh_ligands()

    def select_complex(self, button):
        if self.selected_complex_index == button.complex_index:
            return

        if self.selected_complex:
            self.selected_complex.register_complex_updated_callback(NO_CALLBACK)

        self.reset_selection()
        self.refresh_snapshot_button()

        self.selected_complex_index = button.complex_index
        for item in self.lst_complexes.items:
            btn = item.get_content()
            if btn.selected:
                btn.selected = False
                self.plugin.update_content(btn)
                break
        button.selected = True
        self.plugin.update_content(button)
        self.refresh_selection()

    def select_ligand(self, button):
        if self.selected_ligand == button.ligand:
            return

        self.selected_ligand_index = button.ligand_index
        self.selected_ligand = button.ligand

        for item in self.lst_ligands.items:
            btn = item.get_content()
            if btn.selected:
                btn.selected = False
                self.plugin.update_content(btn)
                break
        button.selected = True
        self.plugin.update_content(button)
        self.compute_results()

    @async_callback
    async def compute_results(self, button=None):
        complex = self.selected_ligand or self.selected_complex

        error = None
        if len(list(complex.atoms)) > MAX_ATOM_COUNT:
            success = False
            error = 'too large'
        else:
            success = self.plugin.helper.prepare_complex(complex)

        self.ln_no_selection.enabled = not success
        self.ln_results.enabled = success
        self.plugin.update_node(self.ln_panel_right)

        if not success:
            text = error or 'rdkit error'
            self.update_preview(text=text)
            return

        self.selected_complex.register_complex_updated_callback(self.refresh_selection)

        self.lbl_complex.text_value = complex.full_name
        self.plugin.update_content(self.lbl_complex)

        self.plugin.helper.add_image(complex)
        self.plugin.helper.add_properties(complex)
        self.plugin.helper.add_smiles(complex)

        self.refresh_results()

        self.update_preview(image=complex.image)

    def refresh_results(self):
        if not self.selected_complex:
            self.update_preview(text='preview')
            self.ln_no_selection.enabled = True
            self.ln_results.enabled = False
            self.plugin.update_node(self.ln_panel_right)
            return

        selected_complex = self.selected_ligand or self.selected_complex
        if not getattr(selected_complex, 'properties', None):
            return

        self.lst_results.items.clear()
        for index in self.plugin.selected_properties:
            prop = selected_complex.properties[index]
            item = self.pfb_result.clone()
            name = item.find_node('Name').get_content()
            name.text_value = prop[0]
            value = item.find_node('Value').get_content()
            value.text_value = prop[2]
            self.lst_results.items.append(item)
        self.plugin.update_content(self.lst_results)

        self.refresh_snapshot_button()

    def refresh_snapshot_button(self):
        snapshot_exists = self.selected_complex in self.plugin.snapshots
        is_ligand = self.selected_ligand is not None
        complex_selected = self.selected_complex_index is not None
        self.btn_snapshot.unusable = is_ligand or not complex_selected or snapshot_exists
        self.plugin.update_content(self.btn_snapshot)

    def update_preview(self, text=None, image=None):
        lbl_enabled = self.ln_message.enabled
        if text:
            self.lbl_message.text_value = text
            if lbl_enabled:
                self.plugin.update_content(self.lbl_message)
            else:
                self.ln_message.enabled = True
                self.ln_frame.enabled = True
                self.ln_image.enabled = False
                self.ln_image.remove_content()
                self.plugin.update_node(self.ln_preview)
        elif image:
            self.ln_message.enabled = False
            self.ln_frame.enabled = False
            self.ln_image.enabled = True
            img = self.ln_image.add_new_image(image)
            img.scaling_option = ScalingOptions.fit
            self.plugin.update_node(self.ln_preview)

    def snapshot_complex(self, button=None):
        complex = self.selected_complex
        complex.timestamp = datetime.now()

        self.plugin.snapshots.append(complex)
        self.refresh_snapshot_button()
        self.plugin.menu_snapshots.refresh_snapshots()
        msg = 'snapshot created for ' + complex.full_name
        self.plugin.send_notification(NotificationTypes.success, msg)
