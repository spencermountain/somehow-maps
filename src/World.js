const fitAspect = require('fit-aspect-ratio')
const htm = require('htm')
const vhtml = require('vhtml')
const d3Geo = require('d3-geo')
const Dot = require('./shapes/Dot')
const Text = require('./shapes/Text')
const Shape = require('./shapes/Shape')
const Line = require('./shapes/Line')
const Graticule = require('./shapes/Graticule')
const data = require('./data')
const Background = require('./Background')

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
    this.back = []
    this.html = htm.bind(vhtml)
    this._clip = true
    this.projection = d3Geo.geoMercator().scale(258)
  }
  mercator() {
    this.projection = d3Geo.geoMercator().scale(450)
  }
  globe() {
    this.projection = d3Geo
      .geoOrthographic()
      .scale(958)
      .translate([-190, -590])
      .rotate([77, -51, 0])
  }
  background(str) {
    let shape = new Background(str, this)
    this.back.push(shape)
    return shape
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
  graticule(obj) {
    let dot = new Graticule(obj, this)
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
  clip(bool) {
    this._clip = bool
    return this
  }
  center(point) {
    if (typeof point === 'string') {
      point = data.points[point]
    }
    this.projection.center(point)
  }
  fit() {}
  build() {
    let h = this.html
    let shapes = this.shapes.sort((a, b) => (a._order > b._order ? 1 : -1))
    let elements = []
    elements = elements.concat(shapes.map(shape => shape.build()))
    elements = elements.concat(this.back.map(sh => sh.build()))
    let attrs = {
      width: this.width,
      height: this.height,
      viewBox: `0,0,${this.width},${this.height}`,
      preserveAspectRatio: 'xMidYMid meet',
      style: 'margin: 10px 20px 25px 25px;' // border:1px solid lightgrey;
    }
    if (this._clip) {
      attrs.style += 'overflow:hidden; border:1px solid #a3a5a5;'
    } else {
      attrs.style += 'overflow:visible;'
    }
    return h`<svg ...${attrs}>
      ${elements}
    </svg>`
  }
}
module.exports = World
