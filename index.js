/* global requestAnimationFrame */
const THREE = require('three')
require('./src/TerrainLoader')
require('./src/TrackballControls')

var width = 800,
  height = 600

var scene = new THREE.Scene()
scene.add(new THREE.AmbientLight(0xeeeeee))
// scene.background = 0xfdfdfd

scene.background = new THREE.Color(0xfdfdfd)

var camera = new THREE.PerspectiveCamera(85, width / height, 0.1, 1000)
camera.position.set(0, -30, 30)

var renderer = new THREE.WebGLRenderer()
renderer.setSize(width, height)

var terrainLoader = new THREE.TerrainLoader()
terrainLoader.load('./assets/jotunheimen.bin', function(data) {
  var geometry = new THREE.PlaneGeometry(60, 60, 199, 199)

  for (var i = 0, l = geometry.vertices.length; i < l; i++) {
    geometry.vertices[i].z = (data[i] / 65535) * 5
  }

  // var material = new THREE.MeshPhongMaterial({
  //   map: THREE.ImageUtils.loadTexture('../assets/jotunheimen-texture.jpg')
  // })
  var material = new THREE.MeshPhongMaterial({
    color: 0xdddddd,
    wireframe: true
  })

  var plane = new THREE.Mesh(geometry, material)
  scene.add(plane)
})

var controls = new THREE.TrackballControls(camera)

document.getElementById('webgl').appendChild(renderer.domElement)

function render() {
  controls.update()
  requestAnimationFrame(render)
  renderer.render(scene, camera)
}

render()
