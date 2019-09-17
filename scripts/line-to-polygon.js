var lineToPolygon = require('turf-line-to-polygon')
// var lineFeature = {
//   type: 'Feature',
//   properties: {},
//   geometry: {
//     type: 'LineString',
//     coordinates: [[0, 0], [0, 1], [1, 1], [1, 0]]
//   }
// }
const lineFeature = require('../data/thames-river.json').features[0]
// console.log(line)

var polyFeature = lineToPolygon(lineFeature)
console.log(JSON.stringify(polyFeature, null, 2))
