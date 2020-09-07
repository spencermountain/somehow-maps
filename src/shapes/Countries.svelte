<script>
  import Shape from './Shape.svelte'
  import * as topojson from 'topojson-client'
  export let color = 'grey'
  // import countries from '../data/countries'
  import countries from '../data/more.geo.json'
  // let key = Object.keys(countries.objects)[0]
  // let geoJSON = topojson.feature(countries, countries.objects[key])
  // console.log(geoJSON)
  // console.log(countries.features[2])
  // console.log(countries.features.findIndex(d => d.properties.name === 'Canada'))
</script>

<Shape shape={countries} stroke={color} />
