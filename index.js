const somehowMaps = require('./src')
const cities = require('/Users/spencer/mountain/somehow-geo/src/data/points/cities.js')
let w = somehowMaps({
  height: 500,
  aspect: 'square'
})

// w.background('states')
w.background('world')
// w.background('rivers')

Object.keys(cities).forEach((k) => {
  if (cities[k][0] > 40 || cities[k][0] < -40) {
    let dot = w.dot().at(cities[k])
    dot.color('blue').radius(4)
  } else {
    // dot.color('lightblue').radius(4)
  }
})
//   ALWAYS USE
//  [ lat, lng ]
//  [ Y, X ]
//  [90->90,   -180, 180]

w.latitude().at(40)
w.latitude().at(-40)

w.clip(false)
w.graticule()
w.globe()

// w.center([-72, 43])
w.rotate(70)
// w.fit()
// w.zoom(2)

let el = document.querySelector('#stage')
el.innerHTML = w.build()

el.addEventListener('resize', () => {
  console.log('SVG resized.')
})
