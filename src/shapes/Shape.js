const GeoJsonLayer = require('@deck.gl/layers').GeoJsonLayer
const colors = require('spencer-color').colors
const scale = require('./fns/_scale')
const uuid = require('./fns/_id')
const toRgb = require('./fns/_rgb')

class Shape {
  constructor(obj, world) {
    this._data = obj.data
    this._color = obj.color || [20, 20, 90]
    this._opacity = obj.opacity || 1
    this._id = obj.id || uuid()
  }
  data(obj) {
    this._data = obj
    return this
  }
  color(c) {
    this._color = toRgb(colors[c]) || toRgb(c)
    return this
  }
  build() {
    console.log(this._color)
    return new GeoJsonLayer({
      id: this._id,
      data: this._data,
      pickable: false,
      stroked: false,
      filled: true,
      extruded: true,
      opacity: this._opacity,
      getFillColor: this._color,
      getElevation: scale(0.01)
    })
  }
}
module.exports = Shape
