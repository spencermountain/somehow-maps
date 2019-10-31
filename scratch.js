const somehowGeo = require('./src')

let w = somehowGeo({
  el: '#stage'
})
w.line({
  data: './data/subway.json'
}).color('red')

w.shape({
  data: './data/great-lakes.json'
}).color('blue')

let timeout = false
// window.resize event listener
window.addEventListener('resize', function() {
  clearTimeout(timeout)
  timeout = setTimeout(() => {
    w.build()
  }, 250)
})

w.build()
