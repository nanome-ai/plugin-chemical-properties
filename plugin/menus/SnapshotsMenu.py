import nanome
from nanome import ui
from nanome.util.enums import NotificationTypes

import os
from datetime import datetime

BASE_DIR = os.path.join(os.path.dirname(__file__))
MENU_PATH = os.path.join(BASE_DIR, 'json/snapshots.json')

class SnapshotsMenu:
    def __init__(self, plugin):
        self.plugin = plugin
        self.snapshots_sort = [None, 0]
        self.create_menu()

    def create_menu(self):
        self.menu = nanome.ui.Menu.io.from_json(MENU_PATH)
        self.menu.index = 2
        self.menu.enabled = False
        root = self.menu.root

        self.pfb_heading: ui.LayoutNode = root.find_node('Prefab Heading')
        self.pfb_value: ui.LayoutNode = root.find_node('Prefab Value')

        self.ln_header: ui.LayoutNode = root.find_node('Header')
        self.lst_snapshots: ui.UIList = root.find_node('Snapshot List').get_content()

        self.btn_export: ui.Button = root.find_node('Export Button').get_content()
        self.btn_export.register_pressed_callback(self.export_snapshots)

    def show_menu(self):
        self.menu.enabled = True
        self.refresh_snapshots()

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

        for heading in self.ln_header.get_children()[1:]:
            btn = heading.get_content()
            old_selected = btn.selected
            btn.selected = btn.prop_index == prop_index

            if old_selected != btn.selected:
                self.plugin.update_content(btn)

        self.refresh_snapshots_values()

    def refresh_snapshots(self, button=None):
        self.refresh_snapshots_header()
        self.refresh_snapshots_values()

    def refresh_snapshots_header(self, button=None):
        self.ln_header.clear_children()

        blank: ui.LayoutNode = self.pfb_value.clone()
        blank.enabled = True
        blank.get_content().text_value = ''
        self.ln_header.add_child(blank)

        for prop_index in sorted(self.plugin.selected_properties):
            heading: ui.LayoutNode = self.pfb_heading.clone()
            heading.enabled = True
            btn: ui.Button = heading.get_content()
            btn.text.value.set_all(self.plugin.helper.short_labels[prop_index])
            btn.prop_index = prop_index
            btn.selected = self.snapshots_sort[0] == prop_index
            btn.register_pressed_callback(self.set_snapshots_sort)
            self.ln_header.add_child(heading)

    def refresh_snapshots_values(self, button=None):
        complexes = self.plugin.snapshots
        if self.snapshots_sort[1] != 0:
            prop_index = self.snapshots_sort[0]
            reverse = self.snapshots_sort[1] == -1

            if self.snapshots_sort[0] == -1: # sort by name first, then timestamp
                sort_fn = lambda complex: (complex.full_name, complex.timestamp)
            else: # sort by property value
                sort_fn = lambda complex: float(complex.properties[prop_index][2])

            complexes = sorted(complexes, key=sort_fn, reverse=reverse)

        self.btn_export.unusable = not self.plugin.snapshots

        self.lst_snapshots.items.clear()
        for complex in complexes:
            ln = nanome.ui.LayoutNode()
            ln.layout_orientation = nanome.util.enums.LayoutTypes.horizontal
            ln.forward_dist = 0.002

            btn = ln.add_new_button('')
            btn.outline.active = False
            btn.mesh.active = True
            btn.mesh.enabled.idle = False
            btn.mesh.color.set_all(nanome.util.Color.Grey())
            btn.tooltip.title = complex.full_name
            # btn.tooltip.content = complex.timestamp.strftime('%Y-%m-%d %H:%M:%S')
            btn.tooltip.bounds.x = 1
            btn.tooltip.bounds.y = 0.25
            btn.tooltip.positioning_target = btn.ToolTipPositioning.bottom_left
            btn.tooltip.positioning_origin = btn.ToolTipPositioning.top
            btn.complex = complex
            btn.register_pressed_callback(lambda b: self.plugin.view_snapshot(b.complex))

            ln_img: ui.LayoutNode = self.pfb_value.clone()
            ln_img.enabled = True
            img = ln_img.add_new_image(complex.thumbnail)
            img.scaling_option = nanome.util.enums.ScalingOptions.fit
            ln.add_child(ln_img)

            for prop_index in sorted(self.plugin.selected_properties):
                prop = complex.properties[prop_index]
                ln_lbl: ui.LayoutNode = self.pfb_value.clone()
                ln_lbl.enabled = True
                lbl: ui.Label = ln_lbl.get_content()
                lbl.text_value = prop.value
                lbl.text_color = prop.color
                ln.add_child(ln_lbl)

            self.lst_snapshots.items.append(ln)

        self.plugin.update_menu(self.menu)

    @nanome.util.async_callback
    async def export_snapshots(self, button=None):
        file = nanome.util.FileSaveData()
        filename = datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '.csv'
        file.path = 'snapshots\\' + filename
        file.write_text(','.join(['NAME', 'SMILES'] + self.plugin.helper.short_labels) + '\n')

        for complex in self.plugin.snapshots:
            values = list(list(zip(*complex.properties))[2])
            file.write_text(','.join([complex.full_name, complex.smiles] + values) + '\n')

        [result] = await self.plugin.save_files([file])
        if result.error_code == nanome.util.FileErrorCode.no_error:
            self.plugin.send_notification(NotificationTypes.success, 'saved to Documents\\nanome\\snapshots')
        else:
            self.plugin.send_notification(NotificationTypes.error, 'error exporting csv')
