// const Shape = require('./Shape')
const colors = require('spencer-color').colors
const d3Geo = require('d3-geo')
const Shape = require('../shapes/Shape')

const defaults = {
  fill: 'none',
  stroke: colors.lighter
}

class Graticule extends Shape {
  constructor(obj = {}, world) {
    obj = Object.assign({}, defaults, obj)
    super(obj, world)
    this.world = world
    this._type = 'Graticule'
    this.lat = 0
    this.lon = 4
  }
  build() {
    let h = this.world.html
    let projection = this.world.projection
    const toPath = d3Geo.geoPath().projection(projection)
    let grat = d3Geo.geoGraticule().step([this.lon, this.lat])
    // grat.extent([0, 90, 0, 90])
    let lines = grat.lines()
    lines = lines.map(line => {
      let d = toPath(line)
      let attrs = Object.assign({}, this.attrs, {
        id: this._id,
        d: d,
        stroke: colors.lighter
      })
      return h`<path ...${attrs}></path>`
    })
    return h`<g>${lines}</g>`
  }
}
module.exports = Graticule
