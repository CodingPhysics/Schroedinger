/********************* SIMULATION OF THE TIME DEPENDENT SCHROEDINGER EQUATION IN 1D *************************
From: 	Coding Phyisics
email:	CodingPhysicSimulation@gmail.com
github:	https://github.com/CodingPhysics
*************************************************************************************************************/

/******************************************** PHYSICAL MODEL ************************************************

We cosider the motion of a single quantum particle in one dimension for a given potential V(x).
The state to this particle is discribed by a complex wave function Ψ(x,t) (probability amplitude).
The 'motion' of the particle is determined by the time development operator U(Δt):

	Ψ(x,t + Δt) = U(Δt) Ψ(x,t),
	U(Δt) = exp(- i H Δt)

where H is the Hamilton operator of the system:

	H = -(1/2) * d^2/dx^2 + V(x)

For a small time step Δt the time development can be approximated by

	(1 + i H Δt/2) Ψ(x,t + Δt) ≈ (1 - i H Δt/2) Ψ(x,t)	(*)
	
*************************************************************************************************************/

/********************************************** ALGORITHM ***************************************************

The continuous wave function is represented by an array of length N:

	Ψ[j] = Ψ(j Δx),	j = 0, ..., N - 1

with Ψ[0] = Ψ[N-1] = 0 for fixed boundary conditions.
The spatial derivative in the Hamilton operator is approximated by

	d^2Ψ/dx^2 (x) ≈ (Ψ[j - 1] - 2 Ψ[j] + Ψ[j])/ Δx^2

Therefore, equation (*) can be rewritten as a martix equation for the vectors Ψ' = Ψ(t + Δt) and Ψ = Ψ(t):

	T Ψ' = - T* Ψ	(**)

where T is an tridiagonal matrix with

	T[j][j-1] = T[j][j-1] = 1
	T[j][j]   = 4 i Δx^2/Δt - 2 Δx^2 V(j Δx) - 2
	T[j][k]   = 0	elsewhere

The matrix equation (**) is solved directly by

	a[j]  = - 1/(T[j][j] + a[j-1]),      a[0]    = 0
	b[j]  = (-T* Ψ)[j] + a[j-1] b[j-1],  b[0]    = 0
	Ψ'[j] = a[j] (Ψ'[j+1] + b[j]),       Ψ'[N-1] = 0

*************************************************************************************************************/

/************************************************ USAGE *****************************************************

Schroedinger is called with a settings object as a single argument, that summarizes all parameters of the simulation:

var mySettings = { 
	potential:     <Function>,     // potential energy V(x)
	energy:        <Number>,       // energy E of the initial wave packet (E = k^2/2)
	median:        <Number>,       // center of the initial wave packet within the interval (0,1)
	sigma:         <Number>,       // width of the initial wave packet
	size:          <Number>        // size N of the wave function array
	timeStep:      <Number>,       // time step Δt of the simulation
	stepsPerFrame: <Number>,       // number of interation steps per frame
	maxFrames:     <Number>,       // total number of simulation steps
	momentumZoom:  <Number>,       // zoom factor for plot of wave function in momentum space 
	scaleFactor:   <Number>,       // scaling factor for the plot of wave function in position space 
	label:         <String>,       // name of the simulation
	underlay:      <p5.Graphics>,  // p5.js graphics buffer for the underlay of the canvas
	imageFile:     <String>|null,  // file name for capturing of the animation frames
	dataFile:      <String>|null   // file name for the simulation data
};

var quantumParticle = new Schroedinger(mySettings);

In order to excute and animate the simulation, the function simulationStep() must be called in the draw function:

function draw() {
	quantumParticle.simulationStep();
}

Image and data files are automatically saved to the local download directory.

*************************************************************************************************************/

function Schroedinger(settings){
	
	this.size          = settings.size;
	this.energy        = settings.energy;
	this.median        = settings.median;
	this.sigma         = settings.sigma;
	this.timeStep      = settings.timeStep;
	this.stepsPerFrame = settings.stepsPerFrame;
	this.maxFrames     = settings.maxFrames;
	this.momentumZoom  = settings.momentumZoom;
	this.scaleFactor   = settings.scaleFactor;
	this.potential     = settings.potential;
	this.label         = settings.label;
	this.underlay      = settings.underlay;
	this.imageFile     = settings.imageFile;
	this.dataFile      = settings.dataFile;
	
	this.xStep         = 1/(this.size - 1);
	this.pStep         = 2*Math.PI;
	this.psiX          = new Array(this.size).fill(new Complex());
	this.psiP          = new Array(this.size).fill(new Complex());
	this.potEnergy     = new Array(this.size).fill(0.);
	this.diagT         = new Array(this.size).fill(0.);
	this.auxA          = new Array(this.size).fill(0.);
	this.rhoX          = new Array(this.size).fill(0.);
	this.rhoP          = new Array(this.size).fill(0.);
	this.plotPot       = new Array(width);
	this.plotX         = new Array(width);
	this.plotP         = new Array(width);
	this.maxPot        = -Infinity;
	this.minPot        = Infinity;
	this.maxFourierAmp = 0.;
	this.frameCount    = 0;
	this.statistics    = {time: 0, normX: 0, normP: 0, meanX: 0, meanP: 0, rmsX: 0, rmsP: 0, leftX: 0, leftP: 0};
	this.dataTable     = null;
	
	/*************************************************************************************************/
	
	this.initialize = function(){
		let start = Math.floor(this.median*(this.size - 1));
		this.computePsiAndV(start, +1);
		this.computePsiAndV(start, -1);
		this.updateMomenta();
		
		this.computeTridiagonalMatrix();
		this.initializePlotParameters();
		
		this.updateStatistics();
		this.show();
		this.frameCount++;
		this.statistics.time += this.stepsPerFrame*this.timeStep;
	} 
	
	/*************************************************************************************************/
	
	this.computePsiAndV = function(startIndex, sign){
		let phase  = 0;
		let decay  = 0;
		
		for(let i = startIndex; (i > 0) && (i < this.size - 1); i += sign){
			let envelope      = Math.exp(-Math.pow((this.xStep*i - this.median)/(2*this.sigma), 2) - decay);
			let phaseFactor   = new Complex(Math.cos(phase), Math.sin(phase)); 
			this.psiX[i]      = phaseFactor.mul(envelope);
			this.rhoX[i]      = this.psiX[i].sqr();
			this.potEnergy[i] = this.potential(this.xStep*i);
			let increment     = sign*Math.sqrt(2*Math.abs(this.energy - this.potEnergy[i]))*this.xStep;
			phase            += (this.potEnergy[i] < this.energy)? increment : 0;
			decay            += (this.potEnergy[i] < this.energy)? 0 : increment;
			this.maxPot       = Math.max(this.maxPot, this.potEnergy[i]);
			this.minPot       = Math.min(this.minPot, this.potEnergy[i]);
		}			
	}
	
	/*************************************************************************************************/
	
	this.updateMomenta = function(){
		for(let i = 0; i < this.size; i++){
			this.psiP[i] = (new Complex()).set(this.psiX[i]);
		}
		fft(this.psiP);		
		for(let i = 0; i < this.size/2; i++){
			[this.psiP[i], this.psiP[i + this.size/2]] = [this.psiP[i + this.size/2], this.psiP[i]];
		}
		
		this.maxFourierAmp = 0.;
		for(let i = 0; i < this.size; i++){
			this.maxFourierAmp = Math.max(this.maxFourierAmp, this.psiP[i].sqr());
			this.rhoP[i] += this.psiP[i].sqr();
		}
		this.maxFourierAmp = Math.sqrt(this.maxFourierAmp);
	}
	
	/*************************************************************************************************/
	
	this.computeTridiagonalMatrix = function(){
		let imaginary = new Complex(0, 4*this.xStep*this.xStep/this.timeStep);
		let scalar    = this.xStep*this.xStep*2;
		
		for(let i = 1; i < this.size-1; i++){
			this.diagT[i] = (new Complex()).add(imaginary,scalar*this.potEnergy[i], 2);
			this.auxA[i] = (i > 1) ? (new Complex(-1)).mul(this.auxA[i-1]) : new Complex();
			this.auxA[i].sub(imaginary).add(scalar*this.potEnergy[i], 2).invert();
		}
	}
	
	/*************************************************************************************************/
	
	this.initializePlotParameters = function(){
		for(let j = 1; j < width-1; j++){
			let x       = j/(width - 1);
			let inArray = { x: this.size*x,           p: this.size*(0.5 + (x - 0.5)/this.momentumZoom) };
			let index   = { x: Math.floor(inArray.x), p: Math.floor(inArray.p)                         };
			let weight  = { x: inArray.x - index.x,   p: inArray.p - index.p                           };
		
			this.plotX[j] = {i: index.x, w: weight.x};
			this.plotP[j] = {i: index.p, w: weight.p};
			this.plotPot[j] = this.potential(x);
		}	
	}
	
	/*************************************************************************************************/
	
	this.updateStatistics = function(){
		
		this.statistics.normX = 0;
		this.statistics.normP = 0;
		
		for(let i = 0; i < this.size/2; i++){
			this.statistics.normX += this.psiX[i].sqr();
			this.statistics.normP += this.psiP[i].sqr();
		}
		
		this.statistics.leftX = this.statistics.normX;
		this.statistics.leftP = this.statistics.normP;
		
		for(let i = this.size/2; i < this.size; i++){
			this.statistics.normX += this.psiX[i].sqr();
			this.statistics.normP += this.psiP[i].sqr();
		}
		
		this.statistics.leftX /= this.statistics.normX;
		this.statistics.leftP /= this.statistics.normP;
		
		this.statistics.meanX = 0;
		this.statistics.meanP = 0;
		
		for(let i = 0; i < this.size; i++){
			let x = i*this.xStep;
			let p = (i - this.size/2)*this.pStep;
			this.statistics.meanX += x*this.psiX[i].sqr();
			this.statistics.meanP += p*this.psiP[i].sqr();
		}
		
		this.statistics.meanX /= this.statistics.normX;
		this.statistics.meanP /= this.statistics.normP;
		
		this.statistics.rmsX = 0;
		this.statistics.rmsP = 0;
		
		for(let i = 0; i < this.size; i++){
			let deltaX = i*this.xStep - this.statistics.meanX;
			let deltaP = (i - this.size/2)*this.pStep - this.statistics.meanP;
			this.statistics.rmsX += deltaX*deltaX*this.psiX[i].sqr();
			this.statistics.rmsP += deltaP*deltaP*this.psiP[i].sqr();
		}
		
		this.statistics.rmsX = Math.sqrt(this.statistics.rmsX/this.statistics.normX);
		this.statistics.rmsP = Math.sqrt(this.statistics.rmsP/this.statistics.normP);
		
		this.appendDataTable();
	}
	
	/*************************************************************************************************/
	
	this.appendDataTable = function(){
		if(this.dataTable == null){
			this.dataTable = new p5.Table();
			
			for(let property in this.statistics){
				this.dataTable.addColumn(property);
			}
		}
		
		let newRow = this.dataTable.addRow();
		for(let property in this.statistics){
			newRow.setString(property, this.statistics[property].toExponential(6)); 
			console.log(property + ': ' + this.statistics[property].toExponential(6)); 
		}
		
		if(this.frameCount === this.maxFrames && this.dataFile !== null){
			saveTable(this.dataTable, this.dataFile + 'Statictics.csv');
			console.log('-> Statictics data saved as ' + this.dataFile + 'Statictics.csv');
		}
	}
	
	/*************************************************************************************************/
	
	this.show = function(){
		background(30);
		
		let textHeight = height/30;
		let topMargin  = 1.5*textHeight;
		let heightX    = 0.7*height;
		let seperator  = heightX + topMargin + 2;
		let heightP    = height - 2*topMargin - heightX - 2;
		
		this.drawUnderlay(seperator,heightX - 2*topMargin, 2*topMargin);
		this.plotComplexData(this.psiX, this.plotX, heightX, topMargin, this.scaleFactor);
		this.plotComplexData(this.psiP, this.plotP, heightP, height - heightP, 1/this.maxFourierAmp);
		this.drawOverlay(seperator, textHeight);
		
		if(this.imageFile !== null){
			saveCanvas(this.imageFile + this.frameCount + '.png');
			console.log('-> Frame saved as ' + this.imageFile + this.frameCount + '.png');
		}
	}
	
	/*************************************************************************************************/
	
	this.drawUnderlay = function(seperatorHeight, plotHeight, topMargin){
		if(this.frameCount == 0){
			let maxPotEnergyWaveVector = Math.sqrt(2*this.maxPot);
			let maxWaveVector          = Math.PI/this.xStep;
			let delta                  = 0.5*width*maxPotEnergyWaveVector/maxWaveVector*this.momentumZoom;
			
			this.underlay.stroke(235);
			this.underlay.strokeWeight(3);
			this.underlay.line(0,seperatorHeight,width,seperatorHeight);
			this.underlay.strokeWeight(1);
			this.underlay.fill(100);
			this.underlay.rect(0.5*width - delta,seperatorHeight, 2*delta, height-seperatorHeight);
			this.underlay.line(0.5*width, seperatorHeight, 0.5*width, height); 
						
			this.underlay.fill(100);
			
			this.underlay.beginShape();
			this.underlay.vertex(0,seperatorHeight-1);
			this.underlay.vertex(0,plotHeight + topMargin);
			let y = plotHeight + topMargin;
			if(this.maxPot > this.minPot){
				for(let j = 1; j < width-1; j++){
					y = plotHeight*(this.maxPot - this.plotPot[j])/(this.maxPot - this.minPot) + topMargin;
					this.underlay.vertex(j, y);
				}
			}
			this.underlay.vertex(width-1,y);
			this.underlay.vertex(width-1,seperatorHeight-1);
			this.underlay.endShape(CLOSE);
		}
		
		image(this.underlay,0,0);
	}
	
	/*************************************************************************************************/
	
	this.plotComplexData = function(data, plotMap, plotHeight, topMargin, scaleFactor){
		strokeWeight(1);
		colorMode(HSB, 255);
		
		for(let j = 1; j < width-1; j++){
			let a   = new Complex(1 - plotMap[j].w);
			let b   = new Complex(plotMap[j].w);
			let i   = plotMap[j].i;
			let val = (new Complex()).add(a.mul(data[i]),b.mul(data[i+1]));
			let y   = plotHeight*(1 - val.abs()*scaleFactor) + topMargin;
			let cl  = Math.floor((val.arg()/Math.PI + 1)*128);
			
			stroke(cl,225,235);
			line(j,plotHeight + topMargin,j,y);
		}
	}
	
	/*************************************************************************************************/
	
	this.drawOverlay = function(seperatorHeight, textHeight){
		noStroke();
		fill(235);
		textSize(textHeight);
		textAlign(LEFT,TOP);
		text('Position Space', 0.25*textHeight, 0.25*textHeight);
		text('Momentum Space (x' + this.momentumZoom + ')', 0.25*textHeight, seperatorHeight + 0.25*textHeight);
		textAlign(RIGHT,TOP);
		text('P(x<.5) = ' + this.statistics.leftX.toFixed(3), 0.41*width, 0.25*textHeight);
		text('P(p<0) = ' + this.statistics.leftP.toFixed(3), 0.41*width, seperatorHeight + 0.25*textHeight);
		textAlign(CENTER,TOP);
		text('<x> = ' + this.statistics.meanX.toPrecision(3), 0.5*width, 0.25*textHeight);
		text('<p> = ' + this.statistics.meanP.toPrecision(3), 0.5*width, seperatorHeight + 0.25*textHeight);
		textAlign(LEFT,TOP);
		text('RMS(x) = ' + this.statistics.rmsX.toPrecision(3), 0.59*width, 0.25*textHeight);
		text('RMS(p) = ' + this.statistics.rmsP.toPrecision(3), 0.59*width, seperatorHeight + 0.25*textHeight);
		textAlign(RIGHT,TOP);
		text('t = ' + this.statistics.time.toFixed(5), width - 0.25*textHeight, 0.25*textHeight);
		text(this.label, width - 0.25*textHeight, seperatorHeight + 0.25*textHeight);
	}
	
	/*************************************************************************************************/
	
	this.propagate = function(){
		let auxB = new Array(this.size);
		
		for(let step = 0; step < this.stepsPerFrame; step++){
			for(let i = 1; i < this.size - 1; i++){
				auxB[i] = (i > 1) ? (new Complex(auxB[i-1])).mul(this.auxA[i-1]) : new Complex();
				auxB[i].add((new Complex(this.diagT[i])).mul(this.psiX[i]).sub(this.psiX[i-1],this.psiX[i+1]));
			}
			for(let i = this.size - 2; i > 0; i--){
				this.psiX[i].set(this.psiX[i+1]).sub(auxB[i]).mul(this.auxA[i]);
				this.rhoX[i] += this.psiX[i].sqr();
			}
		}
	}
	
	/*************************************************************************************************/
	
	this.simulationStep = function(){
		console.log('--- ITERATION ' + this.frameCount + ' ---');
		
		this.propagate();
		this.updateMomenta();
		this.updateStatistics();
		this.show();
		
		this.frameCount++;
		this.statistics.time += this.stepsPerFrame*this.timeStep;
		if(this.frameCount > this.maxFrames){ 
			noLoop();
			if(this.dataFile !== null){
				this.saveAverageDensity();
			}
		}
	}
	
	/*************************************************************************************************/
	
	this.saveAverageDensity = function(){
		
		let densityData = new p5.Table();       
		densityData.addColumn('Position x');
		densityData.addColumn('Density rhoX');
		densityData.addColumn('Momentum p');
		densityData.addColumn('Density rhoP');
		
		let intRhoX = 0;
		let intRhoP = 0;
		for(let i = 0; i < this.size; i++){
			intRhoX += this.rhoX[i]*this.xStep;
			intRhoP += this.rhoP[i]*this.pStep;
		}
		
		for(let i = 0; i < this.size; i++){
			let newRow = densityData.addRow();
			newRow.setString('Position x', (i*this.xStep).toExponential(8));
			newRow.setString('Density rhoX', (this.rhoX[i]/intRhoX).toExponential(8));
			newRow.setString('Momentum p', ((i - this.size/2)*this.pStep).toExponential(8));
			newRow.setString('Density rhoP', (this.rhoP[i]/intRhoP).toExponential(8));
		}
		
		saveTable(densityData, this.dataFile + 'Density.csv');
		console.log('-> Average density saved as ' + this.dataFile + 'Density.csv');
	}
	
	/*************************************************************************************************/
	
	this.initialize();	
}

