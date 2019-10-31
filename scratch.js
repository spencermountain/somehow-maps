const somehowGeo = require('./src')

let w = somehowGeo({
  idel: '#stage',
  width: 500,
  height: 500
})
w.line({
  data: './data/subway.json'
}).color('red')

w.build()
