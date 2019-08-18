const fitAspect = require('fit-aspect-ratio')
const htm = require('htm')
const vhtml = require('vhtml')
const d3Geo = require('d3-geo')
const Dot = require('./shapes/Dot')
const Text = require('./shapes/Text')
const Shape = require('./shapes/Shape')
const Line = require('./shapes/Line')
const Latitude = require('./shapes/Latitude')
const data = require('./data')
const fns = require('./_fns')

const Graticule = require('./background/Graticule')
const Background = require('./background/Background')

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
    this.projection = d3Geo.geoMercator().scale(550)
  }
  fit() {
    let ranges = this.shapes.map(sh => {
      sh.bounds()
    })
    ranges = ranges.filter(o => o)
    let left = fns.bounds(ranges.map(o => o.left)).min
    let top = fns.bounds(ranges.map(o => o.top)).max
    // console.log(min, max)
    return this
  }
  globe() {
    this.projection = d3Geo
      .geoOrthographic()
      .scale(158)
      .translate([90, 90])
      // .scale(958)
      // .translate([-190, -590])
      .rotate([77, -10, 0])
  }
  background(str) {
    let shape = new Background(str, this)
    this.back.push(shape)
    return shape
  }
  graticule(obj) {
    let dot = new Graticule(obj, this)
    this.back.push(dot)
    return dot
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
  latitude(obj) {
    let dot = new Latitude(obj, this)
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
