// Display options:

const CANVAS_WIDTH  = 900//1920;
const CANVAS_HEIGHT = 600//1080;
const FRAME_RATE    = 20;

let settings = {
	size:          1024,
	energy:        3E+4,
	median:        0.5,
	sigma:         0.01,
	timeStep:      1E-6,
	stepsPerFrame: 20,
	maxFrames:     1000,
	potential:     x => 2E+4*Math.pow((4*x - 1)*(4*x - 3),2),
	label:         'Double Well',
	momentumZoom:  4,
	scaleFactor:   1,
	underlay:      null,
	dataFile:      'doubleWell',
	imageFile:     null
};

let quantumParticle;

function setup() {
	frameRate(FRAME_RATE);
	createCanvas(CANVAS_WIDTH, CANVAS_HEIGHT);
	settings.underlay = createGraphics(CANVAS_WIDTH, CANVAS_HEIGHT);
	background(0);
	
	quantumParticle = new Schroedinger(settings);
}

// Draw loop:

function draw() {
	quantumParticle.simulationStep();
}
