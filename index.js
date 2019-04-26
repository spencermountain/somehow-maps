const htm = require('htm')
const h = htm.bind(require('vhtml'))

document.querySelector('#stage').innerHTML = h`<div class="">hello</div>`
