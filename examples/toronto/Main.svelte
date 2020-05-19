<script>
  import { tweened } from 'svelte/motion'
  import { cubicOut } from 'svelte/easing'
  import Road from './Road.svelte'
  import { Map, Globe, Countries, Line, Dot, Shape, Intersections, Graticule } from '../../src'
  import roads from './roads.geo.json'
  import { isDown } from './stores.js'

  let eachRoad = roads.features.map(f => {
    return {
      type: 'FeatureCollection',
      features: [f]
    }
  })
  function onMouseDown() {
    isDown.set(true)
  }
  function onMouseUp() {
    isDown.set(false)
  }
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

<div on:mousedown={onMouseDown} on:mouseup={onMouseUp}>
  <Map focus={roads} tilt={50}>
    <Graticule lat={0.5} lon={0.5} />
    <Countries />
    {#each eachRoad as road}
      <Road shape={road} />
    {/each}
    <!-- <Shape shape={roads} fill="blue" /> -->
    <!-- <Intersections shape={roads} /> -->
  </Map>
</div>
