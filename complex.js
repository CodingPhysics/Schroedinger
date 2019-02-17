/********************************************* COMPLEX NUMBERS **********************************************/

/*************************************************************************************************************
From: 	Coding Phyisics
email:	CodingPhysicSimulation@gmail.com
github:	https://github.com/CodingPhysics
*************************************************************************************************************/

/*************************************************************************************************************
  - cartesian representation (re, im) and polar coordinates (abs, arg)
  - basic aritmetric operations (add, sub, mul, div) with variable arguments
  - operations can be chained
*************************************************************************************************************/


function Complex(a, b) {
	
	this.re = (a instanceof Complex)? a.re : (a ? a : 0.0);
	this.im = (a instanceof Complex)? a.im : (b ? b : 0.0);
	
	this.absVal = undefined;
	this.sqrVal = undefined;
	this.argVal = undefined;
  
	this.sqr = function(){ 
		if(this.sqrVal === undefined){ 
			this.sqrVal = this.re*this.re + this.im*this.im;
		}
		return this.sqrVal;
	}
	
	this.abs = function(){ 
		if(this.absVal === undefined){ 
			this.absVal = Math.sqrt(this.sqr());
		}
		return this.absVal;
	}
	
	this.arg = function(){ 
		if(this.argVal === undefined){
			this.argVal = Math.atan2(this.im, this.re)      
		}
		return this.argVal;
	}
	
	this.clearBuffer = function(){
		this.absVal = undefined;
		this.sqrVal = undefined;
		this.argVal = undefined;
	}
	
	this.toString  = function(){ 
		return this.re.toString()+(this.im < 0 ? '-':'+')+Math.abs(this.im).toString()+'i'
	};
	
	
	this.set = function(x){
		if (typeof x == "number"){
			this.re = x;
			this.im = 0.0;
		} else if(x instanceof Complex) {
			this.re = x.re;
			this.im = x.im;
		}
		this.clearBuffer();
		return this;
	}
	
	this.conj = function(){ 
		this.im = - this.im;
		this.clearBuffer();
		return this;
	}
	
	this.add = function(...args) {
		for(let arg of args){
			if (typeof arg === "number") {
				this.re += arg;
			} else if (arg instanceof Complex) {
				this.re += arg.re;
				this.im += arg.im;
			}
		}
		this.clearBuffer();
		return this;
	}

	this.sub = function(...args) {
		for(let arg of args){
			if (typeof arg === "number") {
				this.re -= arg;
			} else if (arg instanceof Complex){
				this.re -= arg.re;
				this.im -= arg.im;
			}
		}
		this.clearBuffer();
		return this;
	}

	this.mul = function(...args) {
		for(let arg of args){
			if (typeof arg === "number") {
				this.re *= arg;
				this.im *= arg;
			} else if (arg instanceof Complex){
				let a = this.re * arg.re - this.im * arg.im;
				let b = this.re * arg.im + this.im * arg.re;
				this.re = a;
				this.im = b;
			}
		}
		this.clearBuffer();
		return this;
	}

	this.div = function(...args) {
		for(let arg of args){
			if (typeof arg === "number") {
				this.re /= arg;
				this.im /= arg;
			} else if (arg instanceof Complex){
				let a = this.re * arg.re + this.im * arg.im;
				let b = this.im * arg.re - this.re * arg.im;
				let c = arg.re * arg.re + arg.im * arg.im;
				this.re = a / c;
				this.im = b / c;
			}
		}
		this.clearBuffer();
		return this;
	}
	
	this.invert = function(){
		this.re = this.re/this.sqr();
		this.im = -this.im/this.sqr();
		this.clearBuffer();
		return this;
	}
}

function conjugate(x){
	if(x instanceof Complex){
		return x.conj();
	} else if (x instanceof Array){
		let len = x.length;
		for(let i = 0; i < len; i++){
			conjugate(x[i]);
		}
		return x;
	}
}

function complexDotProduct(a, b){
	let dotProduct = new Complex();
	let len        = Math.min(a.length, b.length);
	
	for(let i = 0; i < len; i++){
		dotProduct.add((new Complex()).set(a[i]).conj().mul(b[i]));
	}
	
	return dotProduct;
}