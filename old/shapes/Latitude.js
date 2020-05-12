const Shape = require('./Shape')
const colors = require('spencer-color').colors
const d3Geo = require('d3-geo')

const defaults = {
  fill: 'none',
  stroke: colors.blue,
  'stroke-width': 4
}
class Latitude extends Shape {
  constructor(obj = {}, world) {
    obj = Object.assign({}, defaults, obj)
    super(obj, world)
    this._type = 'Dot'
    this._lng = obj.lng || 5
    this.fill = obj.fill || defaults.fill
  }
  color(c) {
    this.attrs.stroke = colors[c] || c
    return this
  }
  at(lng) {
    this._lng = lng
    return this
  }
  toData() {
    let arr = []
    for (let lat = -180; lat <= 180; lat += 10) {
      arr.push([lat, this._lng]) //create a ton of small line-segments
    }
    return {
      type: 'Feature',
      geometry: {
        type: 'LineString',
        coordinates: arr
      }
    }
  }
  build() {
    let h = this.world.html
    let projection = this.world.projection
    const toPath = d3Geo.geoPath().projection(projection)
    let geoJSON = this.toData()
    let d = toPath(geoJSON)
    let attrs = Object.assign({}, this.attrs, {
      id: this._id,
      d: d,
      fill: this.fill
    })
    return h`<g>
      <path ...${attrs}></path>
    </g>`
  }
}
module.exports = Latitude
