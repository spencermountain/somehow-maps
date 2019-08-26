module.exports = {
  shapes: {
    // lakes: require('./shapes/lakes'),
    // rivers: require('./shapes/rivers'),
    land: require('./shapes/land'),
    'north-america': require('./shapes/north-america'),
    world: require('./shapes/countries'),
    states: require('./shapes/states'),
    'great-lakes': require('./shapes/great-lakes')
  },
  points: require('./points/index'),
  bounds: require('./bounds')
}
