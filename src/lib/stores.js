import { writable } from 'svelte/store'

export let r = writable(0)
export let tilt = writable(0)

// // projection.rotate([rotate, tilt, 3])

// export let proj = derived(height, ($height) => $height * 2)
