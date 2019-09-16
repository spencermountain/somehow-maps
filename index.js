const Deck = require('@deck.gl/core').Deck
const layers = require('@deck.gl/layers')
const GeoJsonLayer = layers.GeoJsonLayer
const PolygonLayer = layers.PolygonLayer

// source: Natural Earth http://www.naturalearthdata.com/ via geojson.xyz
const lakeOntario = './data/lake-ontario.json'
const lakeHuron = './data/lake-huron.json'
const lakeErie = './data/lake-erie.json'
const lakeStClair = './data/lake-st-clair.json'

const INITIAL_VIEW_STATE = {
  latitude: 43.6542,
  longitude: -79.5074,
  zoom: 6.8,
  bearing: -20,
  pitch: 50
}

new Deck({
  initialViewState: INITIAL_VIEW_STATE,
  controller: true,
  layers: [
    new GeoJsonLayer({
      id: 'lakeOntario',
      data: lakeOntario,
      stroked: false,
      filled: true,
      extruded: true,
      lineWidthMinPixels: 2,
      opacity: 1,
      getElevation: 10000,
      getFillColor: [70, 130, 180]
    }),
    new GeoJsonLayer({
      id: 'lakeHuron',
      data: lakeHuron,
      stroked: false,
      filled: true,
      extruded: true,
      lineWidthMinPixels: 2,
      opacity: 1,
      getElevation: 30000,
      getFillColor: [70, 130, 180]
    }),
    new GeoJsonLayer({
      id: 'lakeErie',
      data: lakeErie,
      stroked: false,
      filled: true,
      extruded: true,
      lineWidthMinPixels: 2,
      opacity: 1,
      getElevation: 30000,
      getFillColor: [70, 130, 180]
    }),
    new GeoJsonLayer({
      id: 'lakeStClair',
      data: lakeStClair,
      stroked: false,
      filled: true,
      extruded: true,
      lineWidthMinPixels: 2,
      opacity: 1,
      getElevation: 30000,
      getFillColor: [70, 130, 180]
    })
  ]
})

// document.body.style.margin = '0px'
