const topojson = require('topojson-client')
const topology = require('../data/great-lakes.json')

// var topology = {
//   type: 'Topology',
//   transform: { scale: [1, 1], translate: [0, 0] },
//   objects: { foo: { type: 'Polygon', arcs: [[0]] }, bar: { type: 'Polygon', arcs: [[0, 1]] } },
//   arcs: [[[0, 0], [1, 1]], [[1, 1], [-1, -1]]]
// }
let geojson = topojson.feature(topology, topology.objects.ne_50m_lakes)
console.log(JSON.stringify(geojson))
