import nanome
from nanome.util.enums import NotificationTypes

import os
from datetime import datetime
from functools import partial

BASE_DIR = os.path.join(os.path.dirname(__file__))
MENU_PATH = os.path.join(BASE_DIR, 'json/snapshot.json')
EDIT_ICON = os.path.join(BASE_DIR, '../icons/edit.png')
CHECK_ICON = os.path.join(BASE_DIR, '../icons/check.png')

class SnapshotMenu:
    def __init__(self, plugin, complex, index):
        self.plugin = plugin
        self.complex = complex
        self.index = index
        self.create_menu()

    def create_menu(self):
        self.menu = nanome.ui.Menu.io.from_json(MENU_PATH)
        title = 'Snapshot ' + self.complex.full_name
        title = title[:20] + (title[20:] and '...')
        self.menu.title = title
        self.menu.index = self.index
        root = self.menu.root

        self.pfb_result = root.find_node('Prefab Result')

        self.ln_panel_left = root.find_node('Panel Left')
        self.ln_panel_right = root.find_node('Panel Right')

        self.ln_preview = root.find_node('Complex Preview')
        self.ln_image = root.find_node('Preview Image')
        self.lbl_timestamp = root.find_node('Snapshot Timestamp').get_content()
        self.lbl_smiles = root.find_node('SMILES').get_content()

        self.btn_edit = root.find_node('Button Edit').get_content()
        self.btn_load = root.find_node('Button Load').get_content()
        self.btn_swap = root.find_node('Button Swap').get_content()
        self.btn_delete = root.find_node('Button Delete').get_content()

        self.lst_results = root.find_node('Results List').get_content()

        self.ln_header = root.find_node('Complex Header')
        self.ln_title = root.find_node('Complex Title')
        self.lbl_title = self.ln_title.get_content()

        self.ln_title_input = root.find_node('Input Title')
        self.inp_title = self.ln_title_input.get_content()
        self.inp_title.register_submitted_callback(self.toggle_edit)

        self.btn_edit.register_pressed_callback(self.toggle_edit)
        self.btn_load.register_pressed_callback(partial(self.load_snapshot, swap=False))
        self.btn_swap.register_pressed_callback(partial(self.load_snapshot, swap=True))

        self.btn_edit.icon.active = True
        self.btn_edit.icon.value.set_all(EDIT_ICON)

        self.btn_delete.register_pressed_callback(self.delete_snapshot)
        self.btn_delete.mesh.active = True
        self.btn_delete.mesh.enabled.set_all(True)
        self.btn_delete.mesh.color.set_all(nanome.util.Color(r=220))
        self.btn_delete.mesh.color.highlighted = nanome.util.Color(r=200)
        self.btn_delete.outline.active = False

        self.btn_load.tooltip.title = 'Load as new entry'
        self.btn_load.tooltip.positioning_target = self.btn_load.ToolTipPositioning.bottom
        self.btn_load.tooltip.positioning_origin = self.btn_load.ToolTipPositioning.top
        self.btn_load.tooltip.bounds.y = 0.25
        self.btn_swap.tooltip.title = 'Replace original'
        self.btn_swap.tooltip.positioning_target = self.btn_swap.ToolTipPositioning.bottom
        self.btn_swap.tooltip.positioning_origin = self.btn_swap.ToolTipPositioning.top
        self.btn_swap.tooltip.bounds.y = 0.25

        self.populate_menu()
        self.plugin.update_menu(self.menu)

    def populate_menu(self):
        complex = self.complex
        self.lbl_title.text_value = complex.full_name

        self.lst_results.items.clear()
        for index in self.plugin.selected_properties:
            prop = complex.properties[index]
            item = self.pfb_result.clone()
            name = item.find_node('Name').get_content()
            name.text_value = prop[0]
            value = item.find_node('Value').get_content()
            value.text_value = prop[2]
            self.lst_results.items.append(item)
        self.plugin.update_content(self.lst_results)

        img = self.ln_image.add_new_image(complex.image)
        img.scaling_option = nanome.util.enums.ScalingOptions.fit

        self.lbl_timestamp.text_value = complex.timestamp.strftime('%Y-%m-%d %H:%M:%S')
        self.lbl_smiles.text_value = complex.smiles

    def toggle_edit(self, button=None):
        is_editing = self.ln_title.enabled

        icon = CHECK_ICON if is_editing else EDIT_ICON
        self.btn_edit.icon.value.set_all(icon)

        if is_editing:
            self.ln_title.enabled = False
            self.ln_title_input.enabled = True
            self.inp_title.input_text = self.complex.full_name
        else:
            self.ln_title.enabled = True
            self.ln_title_input.enabled = False
            title = self.inp_title.input_text
            self.complex.full_name = title
            self.lbl_title.text_value = title

            title = 'Snapshot ' + title
            title = title[:20] + (title[20:] and '...')
            self.menu.title = title

            self.plugin.refresh()

        self.plugin.update_menu(self.menu)

    def delete_snapshot(self, button=None):
        self.menu.enabled = False
        self.plugin.update_menu(self.menu)
        self.plugin.snapshots.remove(self.complex)
        self.plugin.refresh()

    def load_snapshot(self, button=None, swap=False):
        def on_complexes(complex, complex_list):
            if complex_list[0]:
                complex.position = complex_list[0].position
                complex.rotation = complex_list[0].rotation
            self.plugin.update_structures_deep([complex])

        complex = self.complex

        if not swap:
            index = complex.index
            complex.index = -1
            self.plugin.update_structures_deep([complex])
            complex.index = index
        else:
            self.plugin.request_complexes([complex.index], partial(on_complexes, complex))
