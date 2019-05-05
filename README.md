# somehow-geo

make a map without thinking

---

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
