const cities = require('./cities')
const ontario = require('./ontario')
const manitoba = require('./north-america')
const countries = require('./countries')
module.exports = Object.assign({}, cities, ontario, manitoba, countries)
