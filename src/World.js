const fitAspect = require('fit-aspect-ratio')
const htm = require('htm')
const vhtml = require('vhtml')
const Dot = require('./shapes/Dot')
const Text = require('./shapes/Text')
const Shape = require('./shapes/Shape')
const Line = require('./shapes/Line')
const d3Geo = require('d3-geo')

class World {
  constructor(obj = {}) {
    this.width = obj.width || 600
    this.height = obj.height || 400
    if (obj.aspect) {
      this.aspect = obj.aspect
      let res = fitAspect(obj)
      this.width = res.width || 600
      this.height = res.height || 400
    }
    this.shapes = []
    this.html = htm.bind(vhtml)

    this.projection = d3Geo
      .geoMercator()
      .scale(1050)
      .center([-79.3961, 43.6601])
  }
  bind(fn) {
    this.html = htm.bind(fn)
  }
  dot(obj) {
    let dot = new Dot(obj, this)
    this.shapes.push(dot)
    return dot
  }
  line(obj) {
    let dot = new Line(obj, this)
    this.shapes.push(dot)
    return dot
  }
  text(obj) {
    let dot = new Text(obj, this)
    this.shapes.push(dot)
    return dot
  }
  shape(obj) {
    let dot = new Shape(obj, this)
    this.shapes.push(dot)
    return dot
  }
  build() {
    let h = this.html
    let shapes = this.shapes.sort((a, b) => (a._order > b._order ? 1 : -1))
    let elements = []
    elements = elements.concat(shapes.map(shape => shape.build()))
    let attrs = {
      width: this.width,
      height: this.height,
      viewBox: `0,0,${this.width},${this.height}`,
      preserveAspectRatio: 'xMidYMid meet',
      style: 'overflow:hidden; margin: 10px 20px 25px 25px;' // border:1px solid lightgrey;
    }
    return h`<svg ...${attrs} class="outline">
      ${elements}
    </svg>`
  }
}
module.exports = World
