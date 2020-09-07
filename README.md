<div align="center">
  <div><b>somehow-maps</b></div>
  <img src="https://user-images.githubusercontent.com/399657/68222691-6597f180-ffb9-11e9-8a32-a7f38aa8bded.png"/>
  <div>— part of <a href="https://github.com/spencermountain/somehow">somehow</a> —</div>
  <div>WIP svelte infographics</div>
  <div align="center">
    <sub>
      by
      <a href="https://spencermounta.in/">Spencer Kelly</a> 
    </sub>
  </div>
</div>
<div align="right">
  <a href="https://npmjs.org/package/somehow-maps">
    <img src="https://img.shields.io/npm/v/somehow-maps.svg?style=flat-square" />
  </a>
</div>
<img height="25px" src="https://user-images.githubusercontent.com/399657/68221862-17ceb980-ffb8-11e9-87d4-7b30b6488f16.png"/>

<div align="center">
  <code>npm install somehow-maps</code>
</div>

this library just hoodwinks some d3 projection methods into Svelte components -

- so that basic SVG maps can be made in a leisurely way.
- **work-in-progress!**

```js
<script>
  import { Globe, Line, Graticule, Dot, Latitude } from 'somehow-map'
</script>

<Globe rotate={30} tilt={-10}>
  <Graticule />
  <Countries stroke="grey" />
  <Latitude at={40} />
  <Line from="toronto" to='jamaica' />
  <Dot at={[90, 0]} color="lightblue" radius={50} />
</Globe>
```

![image](https://user-images.githubusercontent.com/399657/92409352-494f6e80-f10e-11ea-82cf-51c8bdc8e4c0.png)

### Components

these two are the containers, and should be mostly interchangeable

- `<Globe/>`
- `<Map/>`
  both of these accept a 'focus' param, what zooms the map projection to fit a geojson object.

these are the basic primitives

- `<Line from={[lat, lng]} to={'city-name'}/>`
- `<Shape shape={geojson}/>`
- `<Dot radius={4} at={'chicago'}/>`

  these are continual lines around the whole world:

- `<Latitude/>`
- `<Longitude/>`

  these are some helpers

- `<Countries/>`
- `<Graticule/>`
- `<Intersections/>`

### Misc

to generate ad-hoc topojson shapes (without using QGIS):

- grab shapefiles from [natural earth](https://www.naturalearthdata.com/downloads/)
- draw a new shape in [geojson.io](http://geojson.io)
- add them both to [mapshaper.org](https://mapshaper.org/)
- run `mapshaper -clip myshape.geojson`
- export as **topojson**

---

to generate ad-hoc topojson shapes from openstreetmap:

- download some pbf data from [geofabrik.de](http://download.geofabrik.de/)

- use James Halliday's nodejs [pbf parser](https://github.com/substack/osm-pbf-parser) to filter down [specific tags](https://wiki.openstreetmap.org/wiki/Key:amenity)

### See also

- [somehow](https://github.com/spencermountain/somehow)
- [thensome.how](https://thensome.how)

MIT
