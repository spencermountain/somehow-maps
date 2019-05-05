const somehowMaps = require('./src')

let w = somehowMaps({
  height: 300,
  aspect: 'widescreen'
})

w.shape({
  shape: 'great-lakes'
}).fill('lightblue')

w.shape({
  shape: 'provinces'
})
//.color('lightgrey')

w.shape({
  shape: 'north-america'
})

w.line()
  .from('toronto')
  .to('winnipeg')
  .color('red')

w.dot()
  .at('barrie')
  .color('blue')

w.text('Toronto')
  .at('toronto')
  .color('red')

w.clip(true)
w.graticule()
w.globe()

w.center([-72.3961, 43.6601])

document.querySelector('#stage').innerHTML = w.build()
