// const OSMBuildings = require('osmbuildings')
// const OSMBuildings = require('./assets/osmBuildings')

// console.log(OSMBuildings)
// // map.addMapTiles('https://{s}.tiles.mapbox.com/v3/[YOUR_MAPBOX_KEY]/{z}/{x}/{y}.png')

// map.addGeoJSONTiles('https://{s}.data.osmbuildings.org/0.2/anonymous/tile/{z}/{x}/{y}.json')

// https://www.openstreetmap.org/way/15560300#map=17/43.65441/-79.43705
// https://d.data.osmbuildings.org/0.2/anonymous/tile/15/9153/11958.json
var map = new OSMBuildings({
  container: 'map',
  backgroundColor: '#fdfdfd',
  position: { latitude: 43.65453, longitude: -79.43413 },
  zoom: 17.5,

  minZoom: 15,

  maxZoom: 20,

  tilt: 25
})
console.log(map)
const dufferinMall = require('./data/dufferinMall.json')
map.addGeoJSON(dufferinMall)
const cycling = require('./data/subway.json')
map.addGeoJSON(cycling)

// map.addMapTiles(
//   'https://{s}.tiles.mapbox.com/v3/pk.eyJ1Ijoic3BlbmNlcm1vdW50YWluIiwiYSI6Inp5UVZEY3cifQ.dh-_SvkPgv9YOQZLG5ZHKg/{z}/{x}/{y}.png'
// )

// map.addGeoJSONTiles('https://{s}.data.osmbuildings.org/0.2/anonymous/tile/{z}/{x}/{y}.json')
