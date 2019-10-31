const Shape = require('./Shape')
const GeoJsonLayer = require('@deck.gl/layers').GeoJsonLayer

class Line extends Shape {
  constructor(obj, world) {
    super(obj, world)
  }
  build() {
    return new GeoJsonLayer({
      id: this._id,
      data: this._data,
      pickable: false,
      stroked: false,
      filled: true,
      extruded: true,
      lineWidthScale: 10,
      lineWidthMinPixels: 2,
      getFillColor: [160, 160, 180, 200],
      getLineColor: [194, 22, 30],
      getRadius: 100,
      getLineWidth: 4
    })
  }
}
module.exports = Line
