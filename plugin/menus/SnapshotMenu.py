import nanome
from nanome import ui
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

        self.pfb_result: ui.LayoutNode = root.find_node('Prefab Result')

        self.ln_panel_left: ui.LayoutNode = root.find_node('Panel Left')
        self.ln_panel_right: ui.LayoutNode = root.find_node('Panel Right')

        self.ln_preview: ui.LayoutNode = root.find_node('Complex Preview')
        self.ln_image: ui.LayoutNode = root.find_node('Preview Image')
        self.lbl_timestamp: ui.Label = root.find_node('Snapshot Timestamp').get_content()
        self.lbl_smiles: ui.Label = root.find_node('SMILES').get_content()

        self.btn_preview: ui.Button = self.ln_preview.get_content()
        self.btn_edit: ui.Button = root.find_node('Button Edit').get_content()
        self.btn_load: ui.Button = root.find_node('Button Load').get_content()
        self.btn_swap: ui.Button = root.find_node('Button Swap').get_content()
        self.btn_delete: ui.Button = root.find_node('Button Delete').get_content()

        self.lst_results: ui.UIList = root.find_node('Results List').get_content()

        self.ln_header: ui.LayoutNode = root.find_node('Complex Header')
        self.ln_title: ui.LayoutNode = root.find_node('Complex Title')
        self.lbl_title: ui.Label = self.ln_title.get_content()

        self.ln_title_input: ui.LayoutNode = root.find_node('Input Title')
        self.inp_title: ui.TextInput = self.ln_title_input.get_content()
        self.inp_title.register_submitted_callback(self.toggle_edit)

        self.btn_preview.register_pressed_callback(self.open_preview)
        self.btn_edit.register_pressed_callback(self.toggle_edit)
        self.btn_load.register_pressed_callback(partial(self.load_snapshot, swap=False))
        self.btn_swap.register_pressed_callback(partial(self.load_snapshot, swap=True))
        self.btn_delete.register_pressed_callback(self.delete_snapshot)

        self.btn_edit.icon.active = True
        self.btn_edit.icon.value.set_all(EDIT_ICON)

        self.populate_menu()
        self.plugin.update_menu(self.menu)

    def populate_menu(self):
        complex = self.complex
        self.lbl_title.text_value = complex.full_name

        self.lst_results.items.clear()
        for index in sorted(self.plugin.selected_properties):
            prop = complex.properties[index]
            item: ui.LayoutNode = self.pfb_result.clone()
            item.enabled = True
            item.find_node('Name').get_content().text_value = prop.name
            item.get_content().tooltip.content = prop.description
            value = item.find_node('Value').get_content()
            value.text_value = prop.value
            value.text_color = prop.color
            self.lst_results.items.append(item)
        self.plugin.update_content(self.lst_results)

        img = self.ln_image.add_new_image(complex.image)
        img.scaling_option = nanome.util.enums.ScalingOptions.fit

        self.lbl_timestamp.text_value = complex.timestamp.strftime('%Y-%m-%d %H:%M:%S')
        self.lbl_smiles.text_value = complex.smiles

    def open_preview(self, button=None):
        complex = self.complex
        self.plugin.send_files_to_load((complex.image, complex.name))

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

    @nanome.util.async_callback
    async def load_snapshot(self, button=None, swap=False):
        complex = self.complex

        if not swap:
            index = complex.index
            complex.index = -1
            self.plugin.update_structures_deep([complex])
            complex.index = index
        else:
            complex_list = await self.plugin.request_complexes([complex.index])
            if complex_list[0]:
                complex.position = complex_list[0].position
                complex.rotation = complex_list[0].rotation
            self.plugin.update_structures_deep([complex])
