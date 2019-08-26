const somehowMaps = require('./src')

let w = somehowMaps({
  // height: 200,
  // aspect: 'widescreen'
})

// w.background('states')
w.background('world')
// w.background('rivers')

//   ALWAYS USE
//  [ lat, lng ]
//  [ Y, X ]
//  [90->90,   -180, 180]

// toronto
// [43, -79]
// [north +,  west -]

w.line()
  .from('toronto')
  // .to('cape town')
  .to('montreal')
  .color('red')
// w.line()
//   .from('iran')
//   .to('france')
//   .color('blue')
// w.line()
//   .set([[69, -122], [-71, 163]])

// .color('green')
// w.line()
//   .from([-58.3961, -50.6601])
//   .to([-22.3961, -43.6601])
//   .color('orange')

// w.line()
//   .from('vancouver')
//   .to('ghana')
// w.longitude()
//   .at('toronto')
//   .color('blue')

// w.clip(true)
// w.graticule()
// w.globe()

// w.center([-72, 43])
// w.center([43, -79])
// w.center('toronto')
// w.fit('germany')
w.fit()
w.zoom(1.1)

let el = document.querySelector('#stage')
el.innerHTML = w.build()

el.addEventListener('resize', () => {
  console.log('SVG resized.')
})
