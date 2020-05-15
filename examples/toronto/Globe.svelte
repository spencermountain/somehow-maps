<script>
  import { tweened } from 'svelte/motion'
  import { cubicOut } from 'svelte/easing'
  import { Map, Globe, Countries, Line, Dot, Shape, Intersections } from '../../src'
  import roads from './roads.geo.json'

  // import countries from '../../src/data/more.geo.json'
  // let focus = countries.features[3]
</script>

<!-- 
// https://overpass-turbo.eu/
[out:json][timeout:25];
(
  way["highway"="secondary"]({{bbox}});
);
out body;
>;
out skel qt;
 -->

<div>
  <Map focus={roads}>
    <Countries />
    <Shape shape={roads} fill="blue" />
    <Intersections shape={roads} />
  </Map>
</div>
