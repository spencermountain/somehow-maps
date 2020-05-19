<div align="center">
  <img src="https://cloud.githubusercontent.com/assets/399657/23590290/ede73772-01aa-11e7-8915-181ef21027bc.png" />
  <div>make a map in svelte, with minimal thinking</div>

  <a href="https://npmjs.org/package/somehow-maps">
    <img src="https://img.shields.io/npm/v/somehow-maps.svg?style=flat-square" />
  </a>
  <a href="https://unpkg.com/somehow-maps">
    <img src="https://badge-size.herokuapp.com/spencermountain/somehow/master/builds/somehow-maps.min.js" />
  </a>
</div>

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
  <Latitude at={40} />
  <Line from="toronto" to='jamaica' />
  <Dot at={[90, 0]} color="lightblue" radius={50} />
</Globe>
```

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
