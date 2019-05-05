const Shape = require('./Shape')
const defaults = {
  fill: 'grey',
  stroke: 'none',
  'stroke-width': 1,
  'stroke-linecap': 'round'
}

class Text extends Shape {
  constructor(obj, world) {
    let text = null
    let textFn = null
    if (typeof obj === 'string') {
      text = [obj]
      obj = {}
    } else if (typeof obj === 'function') {
      textFn = obj
      obj = {}
    } else if (Array.isArray(obj)) {
      text = obj
      obj = []
    }
    obj = Object.assign({}, defaults, obj)
    super(obj, world)
    this.textFn = textFn
    this.textLines = text || obj.text || []
  }
  build() {
    let h = this.world.html
    let textArr = this.textLines
    if (this.textFn !== null) {
      textArr = this.textFn(this.world)
      textArr = typeof textArr === 'string' ? [textArr] : textArr
    }
    let point = this.world.projection([this.point[1], this.point[0]])
    let inside = textArr.map(str => h`<tspan x="0" dy="1.2em">${String(str)}</tspan>`)
    let transform = `translate(${point[0]} ${point[1]})`
    return h`<g transform="${transform}" style="${this.drawSyle()}">
    <text ...${this.attrs}>
      ${inside}
    </text>
  </g>`
  }
}
module.exports = Text
