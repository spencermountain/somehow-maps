// const Deck = require('@deck.gl/core').Deck
const fitAspect = require('fit-aspect-ratio')
const Shape = require('./shapes/Shape')
const Line = require('./shapes/Line')
const { Deck } = require('@deck.gl/core')

const INITIAL_VIEW_STATE = {
  zoom: 7.0,
  altitude: 1.5,
  height: 1010,
  bearing: 0,
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
    this._el = obj.el || obj.id || 'body'
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
    this._controller = obj.contoller || { scrollZoom: false }
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
    window.el = el
    var canvas = document.createElement('canvas')
    el.appendChild(canvas)
    let gl = canvas.getContext('webgl')

    gl.viewport(0, 0, canvas.width, canvas.height)
    gl.clearColor(1.0, 1.0, 1.0, 1.0)
    gl.clear(gl.COLOR_BUFFER_BIT)
    // ctx.clearColor(255, 255, 255, 1.0)
    // ctx.clear(ctx.COLOR_BUFFER_BIT)

    // build each deckgl layer
    let layers = this.shapes.map(s => s.build())
    console.log({
      width: el.offsetWidth,
      height: el.offsetHeight
    })
    console.log(this._controller)
    new Deck({
      gl: gl,
      width: el.offsetWidth,
      height: el.offsetHeight,
      initialViewState: this.state,
      controller: true,
      // views: new OrthographicView(),
      layers: layers
    })
  }
}

module.exports = World
