const Deck = require('@deck.gl/core').Deck
const GeoJsonLayer = require('@deck.gl/layers').GeoJsonLayer
const scaleLinear = require('./scale')
const properties = require('./properties')

const latitude = 43.6528
const longitude = -79.37

let scale = scaleLinear({
  world: [0, 30000],
  minmax: [0, 180]
})

const INITIAL_VIEW_STATE = {
  latitude: latitude,
  longitude: longitude,
  zoom: 13.8,
  // bearing: -45,
  pitch: 50
}

let layers = [
  { id: 'buildings', path: './data/buildings.json', elevation: 0.2 },
  {
    id: 'greatLakes',
    path: './data/lake-ontario-partial.json',
    elevation: 0.01,
    fill: [91, 131, 186]
  }
]

const color = function() {
  let r = Math.random() * 10
  return [163 + r, 155, 149]
}

layers = layers.map(o => {
  // let prop = properties[o.id] || {}
  return new GeoJsonLayer({
    id: o.id,
    data: o.path,
    stroked: false,
    filled: true,
    extruded: true,
    opacity: 1,
    getElevation: scale(o.elevation || 0.2),
    getFillColor: () => {
      return o.fill || color()
    },
    pickable: true,
    getStrokeColor: [70, 130, 180]
    // onClick: ({ object, x, y }) => {
    //   console.log(object)
    // }
  })
})

layers.push(
  new GeoJsonLayer({
    id: 'subway',
    data: './data/subway.json',
    pickable: false,
    stroked: false,
    filled: true,
    extruded: true,
    lineWidthScale: 10,
    lineWidthMinPixels: 2,
    getFillColor: [160, 160, 180, 200],
    getLineColor: [194, 22, 30],
    getRadius: 100,
    getLineWidth: 4
  })
)
layers.push(
  new GeoJsonLayer({
    id: 'ontario-east',
    data: './data/ontario-line-east.json',
    pickable: false,
    stroked: true,
    filled: true,
    extruded: true,
    lineWidthScale: 10,
    lineWidthMinPixels: 2,
    getLineColor: [97, 34, 156],
    getRadius: 100,
    getLineWidth: 4
  })
)
layers.push(
  new GeoJsonLayer({
    id: 'ontario-west',
    data: './data/ontario-west.json',
    pickable: false,
    stroked: true,
    filled: true,
    extruded: true,
    lineWidthScale: 10,
    lineWidthMinPixels: 2,
    getLineColor: [97, 34, 156],
    getRadius: 100,
    getLineWidth: 4
  })
)

let deck = new Deck({
  initialViewState: INITIAL_VIEW_STATE,
  controller: true,
  layers: layers
})
