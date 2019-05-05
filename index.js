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
}).color('lightgrey')

w.shape({
  shape: 'north-america'
})
// w.shape({
//   shape: 'rivers'
// }).stroke('blue')

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

w.clip(false)

document.querySelector('#stage').innerHTML = w.build()
