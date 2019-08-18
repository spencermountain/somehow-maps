const Shape = require('./Shape')
const colors = require('spencer-color').colors
const d3Geo = require('d3-geo')
const data = require('../data')
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
        coordinates: [this._data[0].reverse(), this._data[1].reverse()]
      }
    }
  }
  from(str) {
    this._data[0] = data.points[str] || str
    return this
  }
  to(str) {
    this._data[1] = data.points[str] || str
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
    return h`<g>
      ${this.makePoint(this._data[0])}
      <path ...${attrs}></path>
      ${this.makePoint(this._data[1])}
    </g>`
  }
}
module.exports = Line
