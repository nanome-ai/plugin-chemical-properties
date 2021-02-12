# Nanome - Chemical Properties

A Nanome plugin to display chemical properties for a complex using rdkit

### Installation

```sh
$ pip install nanome-chemical-properties
```

### Dependencies

This plugin requires `rdkit` to be installed and in the `PATH` variable.

Installation instructions for `rdkit` can be found [here](http://www.rdkit.org/docs/Install.html)

### Usage

To start the plugin:

```sh
$ nanome-chemical-properties -a <plugin_server_address>
```

### Docker Usage

To run in a Docker container:

```sh
$ cd docker
$ ./build.sh
$ ./deploy.sh -a <plugin_server_address>
```

#### Basic properties

In Nanome:

- Activate Plugin
- Select Properties
- Select Complex
- View Results

#### External Properties

This plugin supports configuration to fetch chemical properties from external sources. To get started, rename the `config.example.json` to `config.json`. In that file:

`overwrite` - `true` will replace the default properties with the external ones. `false` will append external properties to the end of the list.\
`endpoints` - list of url endpoints to fetch properties from.

##### `endpoints`

An endpoint is a url that can be queried for chemical properties. Each endpoint can be used to populate multiple properties. The currently supported modes are `GET smiles` and `POST sdf`, with a response type of `json`.

In `GET smiles` mode, the endpoint `url` is expected to contain the string `:smiles`. When this endpoint is queried, the smiles for the chemical in question will replace `:smiles`. Example: `example.com/chem/:smiles` or `example.com/chem?q=:smiles`.

`name` - unique name for the endpoint\
`url` - endpoint url\
`method` - either `GET` or `POST`\
`data` - either `sdf` or `smiles`\
`properties` - a mapping of property names to config

###### `properties`

The properties configuration allows one endpoint to populate multiple properties. Each property has a `path` to find the property in the response body. For example, if the endpoint response body is `{"properties":{"prop1":1}}`, the `path` would be `properties.prop1`.

`description` - description to show up in the "Select Properties" menu\
`format` - format string to format the unit. example `"%.3f"`\
`path` - path to the property in the response body

#### Snapshots

`nanome-chemical-properties` has a "snapshots" feature where you can take a snapshot of a complex to compare properties against other snapshots. Snapshots only persist while the plugin is active, so deactivating the plugin will lose the current snapshots.

To create a snapshot, select a complex and press the "snapshot" button.

To view a table comparison of your snapshots, press the "view snapshots" button.

In the snapshots view:

- Pressing a column header will switch through the sorting options for that column. Pressing it once will sort ascending, a second time will sort descending, and a third time will remove sorting on that column.
- Pressing on a snapshot row will bring up a window containing a 2D rendering of the snapshot, as well as the snapshot timestamp and options to rename, load, or delete the snapshot.
- Pressing "export csv" will save the snapshot names, SMILES string, and properties to ~\Documents\nanome\snaphots on the computer running Nanome.

### License

MIT
