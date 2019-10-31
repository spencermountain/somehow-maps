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
      getLineColor: this._color,
      getRadius: 100,
      getLineWidth: 4
    })
  }
}
module.exports = Line
