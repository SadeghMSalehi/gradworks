## Prun Config File
* particles
	* number-of-subjects = 2
	* number-of-labels = 1
	* number-of-particles = []
	* ignore-sampler-input = true/false
	* ignore-particle-input = true/false
	* use-local-repulsion = true/false
	* use-affine-transform = true/false
	* parameters = ( [ sampling-weight, ensemble-weight, sampling-sigma, sampling-cutoff, intensity-weight ] )
	* bspline-transform
		* control-point-spacing = 8
	* sampler-time-steps = [ initial-time, time-step, final-time ]
	* sampler-initial-cache = "initial-output"
	* sampler-cache = "distribution-output"
	* sampler
		* labels = ( [ intersection-mask, intersection-distance-map, intersection-gradient-map ] )
	* time-steps = [ initial-time, time-step, final -time]