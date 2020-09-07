<script>
  import Shape from './Shape.svelte'
  import * as topojson from 'topojson-client'
  export let stroke = 'grey'
  export let fill = 'none'
  import countries from '../data/more.geo.json'
</script>

<Shape shape={countries} {stroke} {fill} />
