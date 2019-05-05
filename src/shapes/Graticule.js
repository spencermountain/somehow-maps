const Shape = require('./Shape')
const colors = require('spencer-color').colors
const d3Geo = require('d3-geo')

const defaults = {
  fill: 'none',
  stroke: colors.lighter
}

class Graticule extends Shape {
  constructor(obj = {}, world) {
    obj = Object.assign({}, defaults, obj)
    super(obj, world)
  }
  build() {
    let h = this.world.html
    let projection = this.world.projection
    const toPath = d3Geo.geoPath().projection(projection)
    let grat = d3Geo.geoGraticule().step([4, 4])
    let lines = grat.lines()
    lines = lines.map(line => {
      let d = toPath(line)
      let attrs = Object.assign({}, this.attrs, {
        id: this._id,
        d: d
      })
      return h`<path ...${attrs}></path>`
    })
    return h`<g>${lines}</g>`
  }
}
module.exports = Graticule
