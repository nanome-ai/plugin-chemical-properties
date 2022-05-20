import nanome
import os

from nanome import ui

BASE_DIR = os.path.join(os.path.dirname(__file__))
MENU_PATH = os.path.join(BASE_DIR, 'json/settings.json')

class SettingsMenu:
    def __init__(self, plugin: nanome.PluginInstance):
        self.plugin = plugin
        self.create_menu()

    def create_menu(self):
        self.menu = ui.Menu.io.from_json(MENU_PATH)
        self.menu.index = 1
        root: ui.LayoutNode = self.menu.root

        self.dd_properties: ui.Dropdown = root.find_node('Dropdown Properties').get_content()
        self.dd_properties.register_item_clicked_callback(self.toggle_property)

        ln_toggle_colors: ui.LayoutNode = root.find_node('Toggle Colors')
        btn_toggle_colors = ln_toggle_colors.add_new_toggle_switch('Lipinski\'s Rule of 5')
        btn_toggle_colors.register_pressed_callback(self.toggle_colors)
        btn_toggle_colors.selected = self.plugin.enable_colors

        self.populate_properties()

    def show_menu(self):
        self.menu.enabled = True
        self.refresh_properties()
        self.plugin.update_menu(self.menu)

    def toggle_colors(self, btn: ui.Button):
        self.plugin.enable_colors = btn.selected
        self.plugin.refresh()

    def refresh_properties(self):
        for i, ddi in enumerate(self.dd_properties.items):
            ddi.selected = i in self.plugin.selected_properties
        num_selected = len(self.plugin.selected_properties)
        title = f'{num_selected} propert{"ies" if num_selected != 1 else "y"} selected'
        self.dd_properties.permanent_title = title
        self.plugin.update_content(self.dd_properties)

    def populate_properties(self):
        labels = self.plugin.helper.long_labels
        self.dd_properties.items.clear()
        for i, label in enumerate(labels):
            ddi = ui.DropdownItem(label)
            ddi.property_index = i
            ddi.selected = i in self.plugin.selected_properties
            ddi.close_on_selected = False
            self.dd_properties.items.append(ddi)
        self.plugin.update_content(self.dd_properties)

    def toggle_property(self, dd: ui.Dropdown, ddi: ui.DropdownItem):
        if ddi.property_index in self.plugin.selected_properties:
            self.plugin.selected_properties.remove(ddi.property_index)
        else:
            self.plugin.selected_properties.add(ddi.property_index)
        self.refresh_properties()
        self.plugin.refresh()
