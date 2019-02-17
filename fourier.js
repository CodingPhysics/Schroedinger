/**************************************** FAST FOURIER TRANSFORM ********************************************/

/*************************************************************************************************************
From: 	Coding Phyisics
email:	CodingPhysicSimulation@gmail.com
github:	https://github.com/CodingPhysics
*************************************************************************************************************/

/************************************************* USAGE *****************************************************

  - discrete Fourier transform of an array 'amplitudes' (length must be a power of 2)
  - 'amplitudes' is overwritten by the complex fourier coefficients
  - inverse fourier transform for inverse = true
  - requires complex number library (complex.js)

  *************************************************************************************************************/

/***************************************** DISCRETE FOURIER TRANFORM *****************************************

  - given time series a(t) = a(j Δ) = a_j with discrete time t = [0, Δ, ... , (N-1)Δ] 
  - discrete Fourier transform is defined by: b_k = ∑_j a_j exp(-2 π i j k/N)
  - fourier coefficients b(f) = b_k correspond to frequencies 
    f = [0, 1/(N Δ), ... , (N/2 - 1)/Δ, ± 1/Δ, -(N/2 - 1)/Δ, ... , -1/(N Δ)]
  - inverse Fourier transform: a_j = (1/N) * ∑_k b_k exp(2 π i j k/N)

  *************************************************************************************************************/

/******************************************** FFT PSEUDO CODE ************************************************

function fft(N, ampl):
	
	if (N = 1):
	
		return ampl
		
	else:
	
		even = fft(N/2, [ampl[0], ampl[2], ... , ampl[N-2]])
		odd  = fft(N/2, [ampl[1], ampl[3], ... , ampl[N-1]])
		
		for k = 0 to N/2 - 1:
			res[k]       = even[k] + odd[k] * exp(-2*PI*i*k/N)
			res[k + N/2] = even[k] - odd[k] * exp(-2*PI*i*k/N)
		
		return res
		
*************************************************************************************************************/

let memo = []; //data structure for complex exponential functions

function fft(amplitudes, inverse){
	let N = amplitudes.length;
	let Nhalf = N/2;
	if(N == 1){
		return amplitudes;
	} else {
		if(inverse){
			conjugate(amplitudes);
		}
		let even = [];
		let odd  = [];
		for(let i = 0; i < Nhalf; i++){
			even.push(amplitudes[i*2]);
			odd.push(amplitudes[i*2+1]);
		}
		even = fft(even);
		odd  = fft(odd);
		
		let delta = -2*Math.PI/N;
		for(let i = 0; i < Nhalf; i++){
			if(memo[N] === undefined){
				memo[N] = [];
			} 
			if(memo[N][i] === undefined){
				let phase = delta*i;
				memo[N][i] = new Complex(Math.cos(phase), Math.sin(phase));
			}
			let oddPhaseFactor    = (new Complex()).set(memo[N][i]).mul(odd[i]);
			amplitudes[i]         = (new Complex()).set(even[i]).add(oddPhaseFactor);
			amplitudes[i + Nhalf] = (new Complex()).set(even[i]).sub(oddPhaseFactor);
		}
		if(inverse){
			conjugate(amplitudes);
			let Ninverse = 1/N;
			for(let i = 0; i < N; i++){
				amplitudes[i].mul(Ninverse);
			}
		}
		return amplitudes;
	}
}
