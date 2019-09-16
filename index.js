const Deck = require('@deck.gl/core').Deck
const GeoJsonLayer = require('@deck.gl/layers').GeoJsonLayer
const scaleLinear = require('./scale')

const data = {
  lakeOntario: './data/lake-ontario.json',
  lakeHuron: './data/lake-huron.json',
  lakeErie: './data/lake-erie.json',
  lakeStClair: './data/lake-st-clair.json',
  stClairNorth: './data/st-clair-north.json',
  stClairSouth: './data/st-clair-south.json',
  niagaraSouth: './data/niagara-south.json',
  niagaraNorth: './data/niagara-north.json'
}

let scale = scaleLinear({
  world: [0, 30000],
  minmax: [0, 180]
})

// cn tower - 553 m

const INITIAL_VIEW_STATE = {
  latitude: 43.6542,
  longitude: -79.5074,
  zoom: 6.8,
  bearing: 70,
  pitch: 50
}

let layers = [
  { id: 'lakeHuron', elevation: 176 },
  { id: 'stClairNorth', elevation: 175 },
  { id: 'lakeStClair', elevation: 175 },
  { id: 'stClairSouth', elevation: 174 },
  { id: 'lakeErie', elevation: 173 },
  { id: 'niagara-south', elevation: 173 },
  { id: 'niagaraSouth', elevation: 170 },
  { id: 'niagaraNorth', elevation: 76 },
  { id: 'lakeOntario', elevation: 74 }
]
layers = layers.map(o => {
  return new GeoJsonLayer({
    id: o.id,
    data: data[o.id],
    stroked: false,
    filled: true,
    extruded: true,
    lineWidthMinPixels: 2,
    opacity: 1,
    getElevation: scale(o.elevation),
    getFillColor: [70, 130, 180]
  })
})

let deck = new Deck({
  initialViewState: INITIAL_VIEW_STATE,
  controller: true,
  layers: layers
})

setTimeout(() => {
  console.log('hi')
  deck.viewState.transitionDuration = 8000
  deck.viewState.bearing = 20
}, 4000)
// document.body.style.margin = '0px'
