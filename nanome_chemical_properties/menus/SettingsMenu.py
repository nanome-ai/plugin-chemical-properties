import nanome
import os

BASE_DIR = os.path.join(os.path.dirname(__file__))
MENU_PATH = os.path.join(BASE_DIR, 'json/settings.json')

class SettingsMenu:
    def __init__(self, plugin):
        self.plugin = plugin
        self.create_menu()

    def create_menu(self):
        self.menu = nanome.ui.Menu.io.from_json(MENU_PATH)
        self.menu.index = 1
        root = self.menu.root

        self.pfb_property = root.find_node('Prefab Property')
        self.lst_properties = root.find_node('Property List').get_content()

        self.populate_properties()

    def show_menu(self):
        self.menu.enabled = True
        self.refresh_properties()
        self.plugin.update_menu(self.menu)

    def refresh_properties(self):
        for i, item in enumerate(self.lst_properties.items):
            btn = item.get_content()
            btn.selected = i in self.plugin.selected_properties
        self.plugin.update_content(self.lst_properties)

    def populate_properties(self):
        def property_pressed(button):
            button.selected = not button.selected
            self.plugin.update_content(button)

            self.plugin.selected_properties.clear()
            for item in self.lst_properties.items:
                btn = item.get_content()
                if btn.selected:
                    self.plugin.selected_properties.append(btn.property_index)
            self.refresh_properties()
            self.plugin.refresh()

        labels = self.plugin.rdk.long_labels

        self.lst_properties.items.clear()
        for i, label in enumerate(labels):
            item = self.pfb_property.clone()
            btn = item.get_content()
            btn.property_index = i
            btn.selected = i in self.plugin.selected_properties
            btn.register_pressed_callback(property_pressed)
            lbl = item.find_node('Label').get_content()
            lbl.text_value = label
            self.lst_properties.items.append(item)
        self.plugin.update_content(self.lst_properties)
