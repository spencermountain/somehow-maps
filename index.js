const somehowMaps = require('./src')

let w = somehowMaps({
  height: 400,
  aspect: 'widescreen'
})

// w.background('states')
w.background('world')
// w.background('rivers')

// w.line()
//   .from('toronto')
//   .to('winnipeg')
//   .color('red')
let lng = 0
w.line()
  .set([[-180, lng], [-90, lng], [0, lng], [90, lng], [180, lng]])
  .color('red')
  .showPoints(false)

// tropic of cancer
w.latitude()
  .at(23)
  .color('sky')

// tropic of capricorn
w.latitude()
  .at(-23)
  .color('sky')

// w.line()
//   .from([170, 1])
//   .to([20, 1])
//   .color('red')

// w.dot()
//   .at('barrie')
//   .color('blue')

// w.text('Toronto')
//   .at('toronto')
//   .color('red')

w.clip(true)
w.graticule()
w.globe()

w.center([-72.3961, 43.6601])
w.fit()

document.querySelector('#stage').innerHTML = w.build()
