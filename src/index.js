const THREE = require('../assets/three.js')
require('../assets/controls.js')
const drawThreeGeo = require('./three-geojson.js')

const width = 500
const height = 400

//New scene and camera
var scene = new THREE.Scene()
var camera = new THREE.PerspectiveCamera(75, width / height, 0.5, 1000)
//New Renderer
var renderer = new THREE.WebGLRenderer()
renderer.setSize(width, height)
document.getElementById('stage').appendChild(renderer.domElement)

var planet = new THREE.Object3D()
//Create a sphere to make visualization easier.
var geometry = new THREE.SphereGeometry(10, 32, 32)
var material = new THREE.MeshBasicMaterial({
  color: 0x333333,
  wireframe: true,
  transparent: true
})
var sphere = new THREE.Mesh(geometry, material)
planet.add(sphere)
//Draw the GeoJSON
$.getJSON('src/countries.json', function(data) {
  drawThreeGeo(
    data,
    10,
    'sphere',
    {
      color: 0x80ff80
    },
    planet
  )
})
scene.add(planet)
//Set the camera position
camera.position.z = 20
//Enable controls
var controls = new THREE.TrackballControls(camera)
//Render the image
function render() {
  controls.update()
  requestAnimationFrame(render)
  renderer.render(scene, camera)
}
render()
