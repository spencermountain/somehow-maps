import svelte from 'rollup-plugin-svelte'
import resolve from '@rollup/plugin-node-resolve'
import commonjs from '@rollup/plugin-commonjs'
import livereload from 'rollup-plugin-livereload'
import { terser } from 'rollup-plugin-terser'
import json from 'rollup-plugin-json'

const example = 'geo'
const production = false

function serve() {
  let started = false
  return {
    writeBundle() {
      if (!started) {
        started = true
        require('child_process').spawn('npm', ['run', 'serve'], {
          stdio: ['ignore', 'inherit', 'inherit'],
          shell: true
        })
      }
    }
  }
}

export default {
  input: `./examples/${example}/app.js`,
  output: {
    sourcemap: false,
    format: 'iife',
    name: 'app',
    file: 'build/bundle.js'
  },
  plugins: [
    json({}),
    svelte({
      dev: !production,
      css: (css) => {
        css.write('build/bundle.css', false)
      }
    }),
    resolve({
      browser: true,
      dedupe: ['svelte']
    }),
    commonjs(),
    !production && serve(),
    !production && livereload('.'),
    production && terser()
  ],
  watch: {
    clearScreen: false
  }
}
