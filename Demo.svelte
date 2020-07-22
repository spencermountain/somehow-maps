<script>
  import { Globe, Line, Graticule, Dot, Latitude, Countries, Longitude, Map } from './src'
  const toronto = [43.65, -79.43]
</script>

<style>
  .w16 {
    width: 16rem;
  }
</style>

<div class="main">
  <h2>somehow-maps</h2>
  <a href="https://github.com/spencermountain/somehow-maps">somehow-maps</a>
  <div class="w16">
    <Globe tilt={-10} rotate="30" focusOn={'toronto'}>
      <Graticule />
      <Countries color="lightgrey" />
      <Dot at={toronto} radius="2" />
    </Globe>
  </div>
  <hr class="mt5" />
  <div class="w16">
    <!-- <Map tilt={-10}>
      <Graticule />
      <Countries color="lightgrey" />
      <Dot at={toronto} radius="2" />
    </Map> -->
  </div>
</div>
