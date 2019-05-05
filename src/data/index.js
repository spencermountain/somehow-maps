module.exports = {
  shapes: {
    lakes: require('./shapes/lakes'),
    rivers: require('./shapes/rivers'),
    'great-lakes': require('./shapes/great-lakes'),

    world: require('./shapes/world'),
    'north-america': require('./shapes/north-america'),
    countries: require('./shapes/countries'),
    provinces: require('./shapes/provinces')
  },
  points: require('./points/index')
}
