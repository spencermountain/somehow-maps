module.exports = {
  shapes: {
    // lakes: require('./shapes/lakes'),
    // rivers: require('./shapes/rivers'),
    'great-lakes': require('./shapes/great-lakes'),

    land: require('./shapes/land'),
    'north-america': require('./shapes/north-america'),
    world: require('./shapes/countries'),
    states: require('./shapes/states')
  },
  points: require('./points/index')
}
