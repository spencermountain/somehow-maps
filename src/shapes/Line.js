const Shape = require('./Shape')
const colors = require('spencer-color').colors
const d3Geo = require('d3-geo')
const fns = require('../_fns')

const defaults = {
  fill: 'none',
  stroke: colors.blue,
  'stroke-width': 4
}

class Line extends Shape {
  constructor(obj = {}, world) {
    obj = Object.assign({}, defaults, obj)
    super(obj, world)
    this._type = 'Line'
    this._radius = obj.radius || 5
    this._data = []
    this._showPoints = true
  }
  showPoints(bool) {
    this._showPoints = bool
    return this
  }
  color(c) {
    this.attrs.stroke = colors[c] || c
    return this
  }
  // bounds() {
  //   console.log(this._data)
  //   return {}
  // }
  toData() {
    return {
      type: 'Feature',
      geometry: {
        type: 'LineString',
        coordinates: this._data
      }
    }
  }
  from(input) {
    this._data[0] = fns.parsePoint(input)
    return this
  }
  to(input) {
    this._data[1] = fns.parsePoint(input)
    return this
  }
  makePoint(arr) {
    let h = this.world.html
    let point = this.world.projection([arr[0], arr[1]])
    return h`<circle cx="${point[0]}" cy="${point[1]}" r="5" fill="${this.attrs.stroke}"></circle>`
  }
  build() {
    let h = this.world.html
    let projection = this.world.projection
    const toPath = d3Geo.geoPath().projection(projection)
    let geoJSON = this.toData()
    let d = toPath(geoJSON)
    let attrs = Object.assign({}, this.attrs, {
      id: this._id,
      d: d
    })
    let pointA = null
    let pointB = null
    if (this._showPoints) {
      pointA = this.makePoint(this._data[0])
      pointB = this.makePoint(this._data[1])
    }
    return h`<g>
      ${pointA}
      <path ...${attrs}></path>
      ${pointB}
    </g>`
  }
}
module.exports = Line
