import * as THREE from "three";
import "./style.css";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

// SCENE SETUP // SCENE SETUP // SCENE SETUP // SCENE SETUP // SCENE SETUP // SCENE SETUP // SCENE SETUP 

// First we need the basic 3js scene object
const scene = new THREE.Scene();

// Scene background slightly darker than white
scene.background = new THREE.Color(0xf0f0f0);

// Setting up the camera & camerea settings
const camera = new THREE.PerspectiveCamera(
  70, // FOV
  window.innerWidth / window.innerHeight, // Aspect ratio
  1, // Near
  100000 // Far
);

// Set the position & add it to the scene
// Found this poisition is best just through trial & error
camera.position.set(0, 250, 1000);
scene.add(camera);

// LIGHTING // LIGHTING // LIGHTING // LIGHTING // LIGHTING // LIGHTING // LIGHTING // LIGHTING // LIGHTING //

// Will need to tweak this once polishing off the final product

// Ambient + a single shadow-casting Directional is plenty
scene.add(new THREE.AmbientLight(0xf0f0f0, 2.2));

const sun = new THREE.DirectionalLight(0xffffff, 1.6);
sun.position.set(200, 300, 150);
sun.castShadow = true;
sun.shadow.mapSize.set(1024, 1024);
sun.shadow.camera.near = 50;
sun.shadow.camera.far = 2000;
sun.shadow.camera.left = -1200;
sun.shadow.camera.right = 1200;
sun.shadow.camera.top = 1200;
sun.shadow.camera.bottom = -1200;
scene.add(sun);

// FLOOR GRID // FLOOR GRID // FLOOR GRID // FLOOR GRID // FLOOR GRID // FLOOR GRID // FLOOR GRID // FLOOR GRID 

// So the floor plane is necessary as we want the point sphere to cast shadows
const planeWidth = 2000;
const floorY = -200;

const planeGeo = new THREE.PlaneGeometry(planeWidth, planeWidth);
planeGeo.rotateX(-Math.PI / 2);
const planeMat = new THREE.ShadowMaterial({ color: 0x000000, opacity: 0.2 });
// invisible-ish plane to catch shadows
const floor = new THREE.Mesh(planeGeo, planeMat);
floor.position.y = floorY;
floor.receiveShadow = true;
scene.add(floor);

// However I think the floor grid adds clutter
// May ammend this later
// Grid helper on top of the shadow plane
const gridRes = 50;
// const grid = new THREE.GridHelper(planeWidth, gridRes);
// grid.position.y = floorY + 1;
// grid.material.opacity = 0.95;
// grid.material.transparent = true;
// scene.add(grid);

// Renderer / Controls // Renderer / Controls // Renderer / Controls // Renderer / Controls // Renderer / Controls 

// Usual boiler plate, just tripped this out from example

const canvas = document.querySelector(".webgl");
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(window.devicePixelRatio);
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.shadowMap.enabled = true;
renderer.shadowMap.enabled = true;
renderer.shadowMap.type = THREE.PCFSoftShadowMap; 
renderer.render(scene, camera);

const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true;
controls.dampingFactor = 0.2;

// POINT GENERATION // POINT GENERATION // POINT GENERATION // POINT GENERATION // POINT GENERATION // POINT GENERATION 

// Should be independent of anything process specific
// We just want to add sphere at specified location
// Making sure they have properties which allow them to fade after certain amount of time

// reuse geometry/material for perf
const sphereRadius = 10;
const sphereGeo = new THREE.SphereGeometry(sphereRadius, 16, 16);
const sphereMat = new THREE.MeshLambertMaterial({
  color: 0xff0000,
  transparent: true,
  opacity: 1.0
});

function createPoint(now) {
  const x = (Math.random() - 0.5) * planeWidth;
  const z = (Math.random() - 0.5) * planeWidth;
  const mesh = new THREE.Mesh(sphereGeo, sphereMat.clone());
  mesh.position.set(x, floorY + sphereRadius, z);
  mesh.castShadow = true;      // spheres cast shadow
  scene.add(mesh);
  return { mesh, birthTime: now, deathTime: now + lifetime };
}

// POISSON PROCESS // POISSON PROCESS // POISSON PROCESS // POISSON PROCESS // POISSON PROCESS // POISSON PROCESS 

let lambda = 2; // points/sec       
const lifetime = 3; // seconds   
let points = []; // { mesh, birthTime, deathTime }     

function randomPoisson(mean) {
  const L = Math.exp(-mean);
  let p = 1.0, k = 0;
  while (p > L) { k++; p *= Math.random(); }
  return k - 1;
}

function updatePoints(delta, now) {
  const keep = [];
  for (const p of points) {
    const age = now - p.birthTime;
    const lifeFrac = age / lifetime;
    if (now < p.deathTime) {
      if (lifeFrac > 0.8) {              // fade final 20%
        const fade = (lifeFrac - 0.8) / 0.2;
        p.mesh.material.opacity = Math.max(0, 1 - fade);
      }
      keep.push(p);
    } else {
      scene.remove(p.mesh);
      p.mesh.geometry.dispose();
      p.mesh.material.dispose();
    }
  }
  points = keep;
}

// INTENSITY MESH // INTENSITY MESH // INTENSITY MESH // INTENSITY MESH // INTENSITY MESH // INTENSITY MESH 

// geometry lives in XZ-plane, we manipulate vertex Y (height)
const intensityGeo = new THREE.PlaneGeometry(
  planeWidth, planeWidth, gridRes, gridRes
);
intensityGeo.rotateX(-Math.PI / 2);

// red, semi-transparent sheet; no shadow casting/receiving
const intensityMat = new THREE.MeshStandardMaterial({
  color: 0xFFCCCB, // red tint
  emissive: 0xaa0000, // slight self-glow so it's visible
  emissiveIntensity: 0.35,
  roughness: 1,
  metalness: 0,
  side: THREE.DoubleSide,
  transparent: true,
  opacity: 0.45,
  depthWrite: false // DO NOT write depth -> see points beneath
});

const intensityMesh = new THREE.Mesh(intensityGeo, intensityMat);
intensityMesh.position.y = floorY + 10;
intensityMesh.castShadow = false; // do NOT cast shadows
intensityMesh.receiveShadow = true; // and do NOT receive either
intensityMesh.renderOrder = 1;

// If you want to visualize the intensity plane, rather than jsutt he wire frame,
// uncomment this
// Will add back as an option later
// scene.add(intensityMesh);

// a vertex cloud riding on the same geometry (nice for “detail”)
const dotMat = new THREE.PointsMaterial({
  size: 3,
  sizeAttenuation: false,
  color: 0x000000,
  depthWrite: false // keep dots visible through the sheet
});
const dots = new THREE.Points(intensityGeo, dotMat);
dots.position.copy(intensityMesh.position);
dots.rotation.copy(intensityMesh.rotation);
dots.scale.copy(intensityMesh.scale);
dots.castShadow = false; // do NOT cast shadows
dots.renderOrder = 2; // draw after the sheet

// Again, will allow this as an option later

// scene.add(dots);

// Height buffer we animate
// This will act as the local "intensity"
const vertexCount = intensityGeo.attributes.position.count;
const heights = new Float32Array(vertexCount).fill(0);
const posAttr = intensityGeo.attributes.position;

// WIRE FRAME // WIRE FRAME // WIRE FRAME // WIRE FRAME // WIRE FRAME // WIRE FRAME // WIRE FRAME 

const cols = gridRes + 1; // number of vertices per row
const rows = gridRes + 1; // number of rows

// total segments: horizontal + vertical
const totalSegments = rows * (cols - 1) + (rows - 1) * cols;

// positions for 2 endpoints per segment, 3 floats per endpoint
const wirePositions = new Float32Array(totalSegments * 2 * 3);
const wireGeom = new THREE.BufferGeometry();
wireGeom.setAttribute('position', new THREE.BufferAttribute(wirePositions, 3));

const wireMat = new THREE.LineBasicMaterial({
  color: 0x770000,
  transparent: true,
  opacity: 0.55,
  depthWrite: false // draw on top without writing depth
  // (lineWidth is ignored by most browsers/WebGL, so don't rely on it)
});

// We will work based off the intensity mesh
// This might be super inefficient, as we are copying the entire mesh, so essentially 
// doing 2x work for some things. Might be better to change the intensityMesh Material 
// or figure out a better way?

// But for right now, just assume the intensityMesh is up to date with the correct heights?
const wire = new THREE.LineSegments(wireGeom, wireMat);
wire.position.copy(intensityMesh.position);
wire.rotation.copy(intensityMesh.rotation);
wire.scale.copy(intensityMesh.scale);
wire.renderOrder = 2.2; // above dots/sheet
wire.castShadow = false;  // no shadows from the wire
wire.receiveShadow = false;
scene.add(wire);

// Fast refs
const pa = posAttr.array; // Float32Array of [x,y,z,...]
const wp = wirePositions;

// Fill wire endpoints from the surface vertices each frame
function updateWireFromSurface() {
  let ptr = 0;

  // helper which copies vertex i from posAttr into wirePositions
  function pushV(i) {
    const off = i * 3;
    wp[ptr++] = pa[off + 0];
    wp[ptr++] = pa[off + 1];
    wp[ptr++] = pa[off + 2];
  }

  // Horizontal segments (i,j) -> (i+1,j)
  for (let j = 0; j < rows; j++) {
    for (let i = 0; i < cols - 1; i++) {
      const a = j * cols + i;
      const b = a + 1;
      pushV(a);
      pushV(b);
    }
  }

  // Vertical segments (i,j) -> (i,j+1)
  for (let j = 0; j < rows - 1; j++) {
    for (let i = 0; i < cols; i++) {
      const a = j * cols + i;
      const b = a + cols;
      pushV(a);
      pushV(b);
    }
  }

  wireGeom.attributes.position.needsUpdate = true;
  // Don't know wtf this does, still works without it,
  // But it was suggested so keeping it in
  wireGeom.computeBoundingSphere();
}


// cache vertex x,z
const vx = new Float32Array(vertexCount);
const vz = new Float32Array(vertexCount);
for (let i = 0; i < vertexCount; i++) {
  vx[i] = posAttr.getX(i);
  vz[i] = posAttr.getZ(i);
}

function addSpike(x0, z0, strength = 100, radius = 100) {
  const r2 = radius * radius;
  const sigma2 = (radius * radius) / 4;    // (radius/2)^2
  for (let i = 0; i < vertexCount; i++) {
    const dx = vx[i] - x0;
    const dz = vz[i] - z0;
    const d2 = dx * dx + dz * dz;
    if (d2 < r2) {
      const falloff = Math.exp(-d2 / (2 * sigma2));
      heights[i] += strength * falloff;
    }
  }
}

// HAWKES SECTION // HAWKES SECTION // HAWKES SECTION // HAWKES SECTION // HAWKES SECTION // HAWKES SECTION 

// Need to add explanation section in html, maybe a modal, with multiple tabs explaining
// the process simulation in increasing levels of technicality?

// Parameters
let alpha = 0.6;
let beta = 1;
// **baseline TOTAL rate** mu_W (events/sec)
// Notice that I'm using lambda in the code, need to change to mu, the slider has "mu" correctly
let muW = lambda;

// Intensity field components on the grid
// F: history field (sum of k(·-x_i) with time decay, notice no alpha here
const F = new Float32Array(vertexCount).fill(0);

// Baseline surface mu(x) on the grid (here: uniform). 
// Need to change this later when we add non-homogeneous baseline (space)
const AREA = planeWidth * planeWidth; // world-units^2
function muPerArea() { return muW / AREA; } // mu(x) = mu_W / |W| (uniform)

// Pure display scaling
// Actual intensity heights are tiny in magnitude due to large sample area?
// Need to scale them up, however think there
const HEIGHT_SCALE = 5e6;

// Normalized Gaussian stamp for the grid, need to have integrate to 1 or I think 
// potential for things to go wrong, need to investigate why this is
const dx = planeWidth / gridRes;
const dz = planeWidth / gridRes;
const sigmaK = 60; // spatial kernel width in world units (adjust as you like)
const stamp = (function buildGaussianStamp(sigma=80, rad=3) {
  const rx = Math.max(1, Math.ceil((rad * sigma) / dx));
  const rz = Math.max(1, Math.ceil((rad * sigma) / dz));
  const W = 2*rx + 1, H = 2*rz + 1;
  const buf = new Float32Array(W * H);
  let sum = 0;
  for (let di = -rx; di <= rx; di++) {
    for (let dj = -rz; dj <= rz; dj++) {
      const val = Math.exp(-0.5 * ((di*dx)*(di*dx) + (dj*dz)*(dj*dz)) / (sigma*sigma));
      buf[(di+rx)*H + (dj+rz)] = val;
      sum += val;
    }
  }
  const scale = 1 / (sum * dx * dz); // ensures sum buf * dx * dz approx 1
  return { buf, rx, rz, W, H, scale };
})(sigmaK);

function addStampToF(x0, z0) {
  // map world (x0,z0) to grid indices
  const xi = Math.round((x0 + planeWidth/2) / dx);
  const zj = Math.round((z0 + planeWidth/2) / dz);
  for (let di = -stamp.rx; di <= stamp.rx; di++) {
    for (let dj = -stamp.rz; dj <= stamp.rz; dj++) {
      const i = xi + di, j = zj + dj;
      if (i < 0 || i > gridRes || j < 0 || j > gridRes) continue; // clamp edges (no torus)
      const vi = j * (gridRes + 1) + i;
      F[vi] += stamp.buf[(di + stamp.rx) * stamp.H + (dj + stamp.rz)] * stamp.scale;
    }
  }
}

// Mixture sampler pieces
function randn() { // Box–Muller
  const u = 1 - Math.random(), v = Math.random();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}
function sampleFromBaseline() {
  const x = (Math.random() - 0.5) * planeWidth;
  const z = (Math.random() - 0.5) * planeWidth;
  return [x, z];
}
function sampleFromKernel(cx, cz) {
  let x = cx + sigmaK * randn();
  let z = cz + sigmaK * randn();
  // clamp to window (you can switch to wrap-around if you prefer a torus)
  x = Math.max(-planeWidth/2, Math.min(planeWidth/2, x));
  z = Math.max(-planeWidth/2, Math.min(planeWidth/2, z));
  return [x, z];
}
function createPointAt(x, z, now) {
  const mesh = new THREE.Mesh(sphereGeo, sphereMat.clone());
  mesh.position.set(x, floorY + sphereRadius, z);
  mesh.castShadow = true;
  scene.add(mesh);
  return { mesh, birthTime: now, deathTime: now + lifetime };
}

// Time process (Ogata) state
let tSim  = performance.now() * 0.001;  // simulation time we’ve advanced the model to
let Ssum  = 0; // S(t) = sum alpha * e^{-beta(t - T_i)} (sum of parent weights)
let Lambda = muW + Ssum; // total rate big LAMBDA(t−)
let tNext = null; // next accepted event time (scheduled), or null if none
const parents = []; // { x, z, w } with w = current decayed weight alpha * e^{-beta delta_t} (w/ delta_t is change in time )

// decay model state (F, Ssum, parent weights) to a new time
function decayTo(tNew) {
  if (tNew <= tSim) return;
  const dt = tNew - tSim;
  const decay = Math.exp(-beta * dt);
  // spatial field decay
  for (let i = 0; i < vertexCount; i++) F[i] *= decay;
  // weights
  for (let p of parents) p.w *= decay;
  Ssum *= decay;
  Lambda = muW + Ssum;
  tSim = tNew;
}

// schedule the next accepted event time (does not modify state; only tNext)
// This is like the usual hawkes 1D simulation, see R sim comments if you forget my guy
function scheduleNext() {
  if (Lambda <= 0) { tNext = Infinity; return; }
  // We work from the current left-limit (predictable) bigLAMBDA(tSim−) = Lambda
  // Generate candidates until one is accepted
  let M = Lambda;
  let t = tSim;
  for (;;) {
    if (M <= 0) { tNext = Infinity; return; }
    const E = -Math.log(Math.random()) / M; // Exp(M)
    t += E;
    // candidate drift for total rate (no state mutation here)
    const decay = Math.exp(-beta * E);
    const LambdaPrime = muW + (M - muW) * decay; // what bigLAMBDA would be at t (since it only decays)
    // accept?
    if (Math.random() <= LambdaPrime / M) {
      tNext = t; // accepted time
      return;
    } else {
      // rejected, restart from the decayed left-limit at the candidate
      M = LambdaPrime;
      // loop continues
    }
  }
}

// sample event location from the mixture at current tSim (after decaying to tSim)
function sampleLocationNow() {
  const w0 = muW; // baseline weight
  const wS = Ssum; // offspring total weight
  if (w0 + wS <= 0) return sampleFromBaseline();
  if (Math.random() < w0 / (w0 + wS)) return sampleFromBaseline();
  // pick parent proportional to its current weight
  const total = parents.reduce((s, p) => s + p.w, 0);
  let r = Math.random() * total;
  for (const p of parents) { r -= p.w; if (r <= 0) return sampleFromKernel(p.x, p.z); }
  // fallback (numerical)
  const p = parents[parents.length - 1];
  return sampleFromKernel(p.x, p.z);
}

// render/update loop
let oldTime = performance.now() * 0.001;
function loop() {
  const now = performance.now() * 0.001;
  const deltaRender = now - oldTime;
  oldTime = now;

  // If we don't have a scheduled event yet (or after param changes), schedule one
  if (tNext == null) {
    muW = lambda; // keep mu_W synced w/ lambda slider
    scheduleNext();
  }

  // Process any Hawkes events that are due up to 'now'
  while (tNext <= now) {
    // advance model exactly to the event time
    decayTo(tNext);
    // sample location & create visual
    const [x, z] = sampleLocationNow();
    const evt = createPointAt(x, z, tNext);
    points.push(evt);
    // update spatial field and time-process jump
    addStampToF(x, z); // adds one unit of kernel mass to F
    const newParent = { x, z, w: alpha }; // new parent starts with weight alph a
    parents.push(newParent);
    Ssum += alpha;
    Lambda = muW + Ssum;

    // schedule the next event from the new left-limit
    scheduleNext();
  }

  // Finally decay model forward from last processed time to 'now'
  decayTo(now);

  // Push accurate heights to the mesh, height = mu(x) + alpha*F(x)
  const muHere = muPerArea();
  for (let i = 0; i < vertexCount; i++) {
    const h = (muHere + alpha * F[i]) * HEIGHT_SCALE;
    heights[i] = h;
    posAttr.setY(i, h);
  }

  console.log(heights);
  posAttr.needsUpdate = true;
  intensityGeo.computeVertexNormals();
  intensityGeo.attributes.normal.needsUpdate = true;
  updateWireFromSurface();

  // Fade spheres as you already do
  updatePoints(deltaRender, now);

  controls.update();
  renderer.render(scene, camera);
  requestAnimationFrame(loop);
}
loop();

// keep sliders in sync with Hawkes state 
// If baseline slider changes, update mu_W and reschedule from current left-limit
const sliderLambda = document.getElementById("lambdaSlider");
const lambdaValue = document.getElementById("lambdaValue");
if (sliderLambda && lambdaValue) {
  sliderLambda.addEventListener("input", (e) => {
    lambda = parseFloat(e.target.value);
    lambdaValue.textContent = lambda.toFixed(1);
    muW = lambda;
    // Invalidate schedule (optional but safer when parameters jump)
    tNext = null;
  });
}
const sliderBeta = document.getElementById("betaSlider");
const betaValue = document.getElementById("betaValue");
if (sliderBeta && betaValue) {
  sliderBeta.addEventListener("input", (e) => {
    beta = parseFloat(e.target.value);
    betaValue.textContent = beta.toFixed(1);
    tNext = null; // re-schedule with new decay rate
  });
}
const sliderAlpha = document.getElementById("alphaSlider");
const alphaValue = document.getElementById("alphaValue");
if (sliderAlpha && alphaValue) {
  sliderAlpha.addEventListener("input", (e) => {
    alpha = parseFloat(e.target.value);
    alphaValue.textContent = alpha.toFixed(1);
    // no need to reschedule; new events will use the new alpha
    // This needs to be fixed, dont want to spike the old events higher
    // but in the decay portion, I think they are being re-multiplied by alpha, 
    // so the new alpha value spikes them higher than the alpha they
    // were generated with 
  });
}


// =================== Loop ===================
// let oldTime = performance.now() * 0.001;
// const decayRate = 1.0; // s^-1

// function loop() {
//   const now = performance.now() * 0.001;
//   const delta = now - oldTime;
//   oldTime = now;

//   // arrivals
//   const n = randomPoisson(lambda * delta);
//   for (let i = 0; i < n; i++) {
//     const p = createPoint(now);
//     points.push(p);
//     addSpike(p.mesh.position.x, p.mesh.position.z);
//   }

//   // decay & push heights into geometry
//   const decay = Math.exp(-decayRate * delta);
//   for (let i = 0; i < vertexCount; i++) {
//     heights[i] *= decay;
//     posAttr.setY(i, heights[i]);
//   }
//   posAttr.needsUpdate = true;
//   intensityGeo.computeVertexNormals();
//   intensityGeo.attributes.normal.needsUpdate = true;

//   updateWireFromSurface();

//   updatePoints(delta, now);

//   controls.update();
//   renderer.render(scene, camera);
//   requestAnimationFrame(loop);
// }
// loop();


// =================== Resize ===================
window.addEventListener("resize", () => {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
});


