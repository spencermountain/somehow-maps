const somehowMaps = require('./src')

let w = somehowMaps({
  height: 300,
  aspect: 'widescreen'
})

w.dot()
  .at('toronto')
  .color('red')

w.text('toronto')
  .at('toronto')
  .color('red')

w.shape({
  shape: 'great-lakes'
})

w.line()
  .from('toronto')
  .to('winnipeg')
  .color('red')

document.querySelector('#stage').innerHTML = w.build()
