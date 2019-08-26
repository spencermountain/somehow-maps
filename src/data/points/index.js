module.exports = Object.assign(
  {},
  require('./cities'),
  require('./ontario'),
  require('./manitoba'),
  require('./countries')
)

console.log(module.exports['france'])
