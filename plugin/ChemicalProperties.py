import nanome
from nanome.util.enums import Integrations

from .menus import MainMenu, SettingsMenu, SnapshotMenu, SnapshotsMenu
from .PropertiesHelper import PropertiesHelper

class ChemicalProperties(nanome.AsyncPluginInstance):
    def start(self):
        self.set_plugin_list_button(self.PluginListButtonType.run, 'Open')
        self.set_plugin_list_button(self.PluginListButtonType.advanced_settings, 'Select Properties')

        self.helper = PropertiesHelper()
        self.selected_properties = list(range(self.helper.num_props))

        self.menu_main = MainMenu(self)
        self.menu_settings = SettingsMenu(self)
        self.menu_snapshots = SnapshotsMenu(self)

        self.snapshots = []
        self.snapshot_menu_index = 3
        self.snapshot_menus = []

        self.on_run()

    def on_stop(self):
        del self.helper

    def on_run(self):
        self.menu_main.show_menu()

    def on_complex_added(self):
        self.menu_main.refresh_complexes()

    def on_complex_removed(self):
        self.menu_main.refresh_complexes()

    def on_advanced_settings(self):
        self.menu_settings.show_menu()

    def refresh(self):
        self.menu_main.refresh_results()
        self.menu_snapshots.refresh_snapshots()

    def view_snapshot(self, complex):
        # reopen menu if already exists
        for menu in self.snapshot_menus:
            if complex == menu.complex:
                menu.menu.enabled = True
                self.update_menu(menu.menu)
                return

        # create new menu
        menu_snapshot = SnapshotMenu(self, complex, self.snapshot_menu_index)
        self.snapshot_menus.append(menu_snapshot)
        self.snapshot_menu_index += 1

def main():
    plugin = nanome.Plugin("Chemical Properties", "Calculates and displays different properties of chemicals using the RDKit Python library", "Analysis", True, integrations=[Integrations.analysis])
    plugin.set_plugin_class(ChemicalProperties)
    plugin.run()

if __name__ == "__main__":
    main()
