const Deck = require('@deck.gl/core').Deck
const fitAspect = require('fit-aspect-ratio')
const Shape = require('./shapes/Shape')
const Line = require('./shapes/Line')

const INITIAL_VIEW_STATE = {
  zoom: 12.0,
  altitude: 1.5,
  bearing: -66.66,
  height: 1010,
  latitude: 43.66,
  longitude: -79.36,
  maxPitch: 60,
  maxZoom: 20,
  minPitch: 0,
  minZoom: 0,
  pitch: 50
}

class World {
  constructor(obj = {}) {
    this.width = 100
    this.height = 100
    this._el = obj._el || obj.id || 'body'
    this.state = Object.assign({}, INITIAL_VIEW_STATE, obj)
    obj.aspect = obj.aspect || 'widescreen'
    if (obj.aspect) {
      this.aspect = obj.aspect
      let res = fitAspect({
        aspect: obj.aspect,
        width: 100
      })
      this.width = res.width || 100
      this.height = res.height || 100
    }
    this._controller = true
    this.shapes = []
    if (typeof this.el === 'string') {
      this.el = document.querySelector(this.el)
    }
  }
  shape(obj) {
    let shape = new Shape(obj)
    this.shapes.push(shape)
    return shape
  }
  line(obj) {
    let shape = new Line(obj)
    this.shapes.push(shape)
    return shape
  }
  controller(obj) {
    this._controller = obj
    return this
  }
  build() {
    let el = document.querySelector(this._el)
    var canvas = document.createElement('canvas')
    el.appendChild(canvas)
    let ctx = canvas.getContext('webgl')
    // build each deckgl layer
    let layers = this.shapes.map(s => s.build())
    new Deck({
      gl: ctx,
      width: el.offsetWidth,
      height: el.offsetHeight,
      initialViewState: this.state,
      controller: this._controller,
      layers: layers
    })
  }
}

module.exports = World
