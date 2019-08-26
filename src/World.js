const fitAspect = require('fit-aspect-ratio')
const htm = require('htm')
const vhtml = require('vhtml')
const d3Geo = require('d3-geo')
const Dot = require('./shapes/Dot')
const Text = require('./shapes/Text')
const Shape = require('./shapes/Shape')
const Line = require('./shapes/Line')
const Latitude = require('./shapes/Latitude')
const Longitude = require('./shapes/Longitude')
const data = require('./data')
const fns = require('./_fns')

const Graticule = require('./background/Graticule')
const Background = require('./background/Background')

class World {
  constructor(obj = {}) {
    this.width = obj.width || 500
    this.height = obj.height || 500
    if (obj.aspect) {
      this.aspect = obj.aspect
      let res = fitAspect(obj)
      this.width = res.width || 500
      this.height = res.height || 500
    }
    this.shapes = []
    this.back = []
    this.html = htm.bind(vhtml)
    this._clip = true
    this.projection = d3Geo.geoMercator() //.fitSize([500, 500], feature) //.scale(350)
  }
  mercator() {
    this.projection = d3Geo.geoMercator().scale(200)
  }
  globe() {
    this.projection = d3Geo.geoOrthographic().scale(288)
    // .translate([90, 90])
    // .scale(958)
    // .translate([-190, -590])
    // .rotate([29, -10, 0])
    //       .rotate([29, -20, 0])
  }
  fit(input) {
    let box = []
    if (input) {
      box = fns.parseBounds(input)
    } else {
      let ranges = this.shapes.map(sh => sh.bounds())
      ranges = ranges.filter(o => o)
      let east = fns.bounds(ranges.map(o => o.east)).min
      let west = fns.bounds(ranges.map(o => o.west)).max
      let top = fns.bounds(ranges.map(o => o.top)).min
      let bottom = fns.bounds(ranges.map(o => o.bottom)).max
      // give it margins
      // east = east + 13
      // console.log(east)
      // top *= 1.01
      // bottom *= 1.1
      // east += 5
      box = [[west, top], [east, bottom]]
    }
    // box = [[69, -122], [-71, 43]]
    // box = [[80, 80], [-170, 170]]

    // console.log(box)
    // box = [[-9.0882278, 72.2460938].reverse(), [-55.3228175, 168.2249543].reverse()]
    // box = fns.parseBounds(input)
    let shape = {
      type: 'Feature',
      geometry: {
        type: 'LineString',
        coordinates: [box[0].reverse(), box[1].reverse()]
      }
    }
    // this.projection.fitSize([this.width - 10, this.height - 10], shape)
    let margin = 25
    let extent = [[margin, 0], [500 - margin, 500]]
    this.projection.fitExtent(extent, shape)

    return this
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

  clip(bool) {
    this._clip = bool
    return this
  }
  center(input) {
    let point = fns.parsePoint(input)
    this.projection.center(point)
  }
  zoom(mult) {
    // let scale = fns.parseZoom(input)
    // this.projection.scale(scale)
    let current = this.projection.scale()
    this.projection.scale(current * mult)
    return this
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
  longitude(obj) {
    let dot = new Longitude(obj, this)
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
    elements = elements.concat(this.back.map(sh => sh.build()))
    let attrs = {
      // width: this.width,
      // height: this.height,
      viewBox: `0,0,${this.width},${this.height}`,
      preserveAspectRatio: 'xMidYMid meet',
      style: 'margin: 10px 20px 25px 25px;' // border:1px solid lightgrey;
    }
    if (this._clip) {
      attrs.style += 'overflow:hidden; border:1px solid #a3a5a5;'
    } else {
      attrs.style += 'overflow:visible; border:1px solid #a3a5a5;'
    }
    return h`<svg ...${attrs}>
      ${elements}
    </svg>`
  }
}
module.exports = World
