const somehowMaps = require('./src')

let w = somehowMaps({
  height: 300,
  aspect: 'widescreen'
})

w.shape({
  shape: 'great-lakes'
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

document.querySelector('#stage').innerHTML = w.build()
