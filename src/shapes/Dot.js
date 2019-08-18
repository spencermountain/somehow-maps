const Shape = require('./Shape')
const colors = require('spencer-color').colors

const defaults = {
  fill: colors.blue,
  stroke: 'none'
}
class Dot extends Shape {
  constructor(obj = {}, world) {
    obj = Object.assign({}, defaults, obj)
    super(obj, world)
    this._type = 'Dot'
    this._radius = obj.radius || 5
  }
  build() {
    let h = this.world.html
    let point = this.world.projection([this.point[1], this.point[0]])
    let attrs = Object.assign({}, this.attrs, {
      id: this._id,
      cx: point[0],
      cy: point[1],
      r: this._radius
    })
    return h`<circle ...${attrs} ><title>${this._title}</title></circle>`
  }
}
module.exports = Dot
