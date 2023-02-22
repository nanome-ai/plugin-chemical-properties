# Nanome - Chemical Properties

A Nanome Plugin to display chemical properties for a complex using RDKit.

## Dependencies

[Docker](https://docs.docker.com/get-docker/)

When running outside of Docker:

This plugin requires `rdkit` to be installed and in the `PATH` variable.

Installation instructions for `rdkit` can be found [here](http://www.rdkit.org/docs/Install.html)

## Usage

To run in a Docker container:

```sh
$ cd docker
$ ./build.sh
$ ./deploy.sh -a <plugin_server_address>
```

## Basic Properties

In Nanome:

- Activate Plugin
- Select Properties
- Select Complex
- View Results

## External Properties

This plugin supports configuration to fetch chemical properties from external sources. To get started, rename the `config.example.json` to `config.json`. In that file:

`overwrite` - `true` will replace the default properties with the external ones. `false` will append external properties to the end of the list.\
`endpoints` - list of url endpoints to fetch properties from.

### `endpoints`

An endpoint is a url that can be queried for chemical properties. Each endpoint can be used to populate multiple properties. The currently supported modes are `GET smiles`, `POST smiles`, and `POST sdf`, with a response type of `json`.

In `GET smiles` mode, the endpoint `url` is expected to contain the string `:smiles`. When this endpoint is queried, the smiles for the chemical in question will replace `:smiles`. Example: `example.com/chem/:smiles` or `example.com/chem?q=:smiles`.

In `POST smiles` mode, the `payload` is expected to contain the string `:smiles`. When this endpoint is queried, the smiles for the chemical in question will replace `:smiles`. Example: `{"smiles":":smiles"}`.

`name` - unique name for the endpoint\
`url` - endpoint url\
`method` - either `GET` or `POST`\
`data` - either `sdf` or `smiles`\
`cache_time` - time in seconds to cache result for same SMILES (default 30).\
`headers` - optional headers object to send with the request\
`payload` - for `POST smiles` requests, the payload to send containing `:smiles`\
`properties` - a mapping of property names to config

#### `properties`

The properties configuration allows one endpoint to populate multiple properties. Each property has a `path` to find the property in the response body. For example, if the endpoint response body is `{"properties":{"prop1":1}}`, the `path` would be `properties.prop1`.

`description` - description to show up in the "Select Properties" menu\
`format` - format string to format the unit. example `"%.3f"`\
`path` - path to the property in the response body
`color` - optional color scheme to apply to the property, read below

##### `path` syntax

`path` supports dot notation, array indexing, and array searching.

`path.to.prop` - dot notation\
`path[0]` - array indexing\
`path[key=value]` - array searching

Dot notation: if the endpoint response body is `{"properties":{"prop1":1}}`, the `path` would be `properties.prop1`.

Array indexing: if the endpoint response body is `{"properties":[1,2]}`, the `path` would be `properties[0]`.

Array searching: if the endpoint response body is `{"properties":[{"key":"prop1","value":1},{"key":"prop2","value":2}]}`, the `path` would be `properties[key=prop1].value`.

##### `property color`

To colorize external properties, you can add the `color` config to a property. The two color `type`s are `within` and `gradient`, where `within` colors a property white if it is within the provided range and red otherwise, and `gradient` colors a property according to a gradient with the provided color stops. Examples are provided in the `config.example.json` file.

`within` `args`:\
  `min` - optional min bound\
  `max` - optional max bound

`gradient` `args`:\
  `colors` - list of `[value, color]` pairs where `color` is a hexadecimal RGBA value (e.g. `"8000FFFF"` for purple). must be sorted ascending `value`

## Snapshots

Chemical Properties has a "snapshots" feature where you can take a snapshot of a complex to compare properties against other snapshots. Snapshots only persist while the plugin is active, so deactivating the plugin will lose the current snapshots.

To create a snapshot, select a complex and press the "snapshot" button.

To view a table comparison of your snapshots, press the "view snapshots" button.

In the snapshots view:

- Pressing a column header will switch through the sorting options for that column. Pressing it once will sort ascending, a second time will sort descending, and a third time will remove sorting on that column.
- Pressing on a snapshot row will bring up a window containing a 2D rendering of the snapshot, as well as the snapshot timestamp and options to rename, load, or delete the snapshot.
- Pressing "export csv" will save the snapshot names, SMILES string, and properties to ~\Documents\nanome\snaphots on the computer running Nanome.

## Development

To run Chemical Properties with autoreload:

```sh
$ python3 -m pip install -r requirements.txt
$ python3 run.py -r -a <plugin_server_address> [optional args]
```

## License

MIT
