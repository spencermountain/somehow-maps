const data = require('../data')
const fns = require('../_fns')
const d3Geo = require('d3-geo')
const topojson = require('topojson-client')
const colors = require('spencer-color').colors

const defaults = {
  fill: 'none',
  stroke: 'grey'
}
class Shape {
  constructor(obj = {}, world) {
    this.world = world
    this._data = obj.data || []
    this._id = obj.id //|| fns.uid('input')
    this.attrs = Object.assign({}, defaults, obj)
    this.shape = data.shapes[obj.shape] || data.points[obj.point] || obj.shape
    this.point = []
    this.style = {}
    this._type = 'Shape'
  }
  at(str) {
    if (typeof str === 'string') {
      str = str.toLowerCase().trim()
      this.point = fns.parsePoint(str)
    } else {
      this.point = str
    }
    return this
  }
  set(input) {
    this._data = input
    return this
  }
  bounds() {
    let x = fns.bounds(this._data.map(a => a[0]))
    let y = fns.bounds(this._data.map(a => a[1]))
    return {
      top: x.min,
      bottom: x.max,
      east: y.min,
      west: y.max
    }
  }
  color(c) {
    this.attrs.fill = colors[c] || c
    return this
  }
  stroke(c) {
    this.attrs.stroke = colors[c] || c
    return this
  }
  fill(c) {
    this.attrs.fill = colors[c] || c
    return this
  }
  opacity(n) {
    this.attrs.opacity = n
    return this
  }
  drawSyle() {
    return Object.keys(this.style)
      .map(k => {
        return `${k}:${this.style[k]};`
      })
      .join(' ')
  }
  geoJSON() {
    if (typeof this.shape === 'object') {
      let key = Object.keys(this.shape.objects)[0]
      return topojson.feature(this.shape, this.shape.objects[key])
    }
    return []
  }
  build() {
    let h = this.world.html
    let projection = this.world.projection
    const toPath = d3Geo.geoPath().projection(projection)
    let geojson = this.geoJSON()
    let d = toPath(geojson)
    let attrs = Object.assign({}, this.attrs, {
      id: this._id,
      d: d
    })
    return h`<path ...${attrs}></path>`
  }
}
module.exports = Shape
