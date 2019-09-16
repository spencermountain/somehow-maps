const Deck = require('@deck.gl/core').Deck
const layers = require('@deck.gl/layers')
const GeoJsonLayer = layers.GeoJsonLayer
const PolygonLayer = layers.PolygonLayer

// source: Natural Earth http://www.naturalearthdata.com/ via geojson.xyz
const COUNTRIES =
  'https://d2ad6b4ur7yvpq.cloudfront.net/naturalearth-3.3.0/ne_50m_admin_0_scale_rank.geojson' //eslint-disable-line
const subway = './data/subway.json'

const INITIAL_VIEW_STATE = {
  // latitude: 49.2573,
  // longitude: -123.1241,
  latitude: 43.6542,
  longitude: -79.5074,
  zoom: 10,
  bearing: 0,
  pitch: 30
}

const landCover = [[[-79.0, 43.7], [-74.02, 40.7], [-74.02, 40.72], [-74.0, 40.72]]]
// 'https://raw.githubusercontent.com/uber-common/deck.gl-data/master/examples/geojson/vancouver-blocks.json'
const deck = new Deck({
  initialViewState: INITIAL_VIEW_STATE,
  controller: true,
  layers: [
    new PolygonLayer({
      id: 'ground',
      data: landCover,
      stroked: false,
      getPolygon: f => f,
      getFillColor: [0, 0, 0, 0]
    }),
    new GeoJsonLayer({
      id: 'geojson',
      data: './data/dufferinMall.json',
      opacity: 1,
      stroked: true,
      filled: true,
      extruded: true,
      wireframe: false,
      getElevation: f => {
        console.log(f)
        return 100
      },
      getFillColor: f => [120, 120, 120],
      getLineColor: [0, 0, 0],
      pickable: true,
      onHover: this._onHover
    }),

    new GeoJsonLayer({
      id: 'base-map',
      data: COUNTRIES,
      // Styles
      stroked: true,
      filled: true,
      lineWidthMinPixels: 2,
      opacity: 0.4,
      // getLineDashArray: [3, 3],
      getLineColor: [60, 60, 60],
      getFillColor: [200, 200, 200]
    }),
    new GeoJsonLayer({
      id: 'subway',
      data: subway,
      // Styles
      stroked: true,
      filled: true,
      lineWidthMinPixels: 2,
      opacity: 0.4,
      // getLineDashArray: [3, 3],
      getLineColor: [60, 60, 60],
      getFillColor: [200, 200, 200]
    })
  ]
})

// For automated test cases
/* global document */
document.body.style.margin = '0px'
