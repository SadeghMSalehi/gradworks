./prun2 with dimension = 2
## ParticleRun Command Line Options
* -o
	* Specify a filename for an output image
* --fusion
	* label fusion from a config
	* *ex)* `--fusion config-file output-file target-image`
* --overlap
	* Compute the overlap ratio (dice|jaccard). This option can take two or arguments. The first argument is a gold standard, and other arguments are multiple number of label images to be compared.
	* *ex)* --overlap dice output-text ref1 ref2-1 ref2-2 ... ref2-n
* --p2mat
	* point list to matrix
* --slice
	* extract a slice from 3d volume
	* *ex)* --slice dim index imagefile outputfile
* --imageMerge
	* merge 2D images into a 3d volume (--imageMerge output input1 input2 ...)
* --qa
	* extract a slice with a label map
* --config
	* [file] use this config file
* --demons
	* run Demons registration
* --separate
	* [input] [x] [y] [z] ... separate vector images into individual image files
* --rx
	* registration experiments 
* --dots
	* --rx --dots generate a series of gaussian dot images
* --sigma
	* sigma value [double]
	* *ex)* --sigma 0.8
* --entropyImage
	* Compute an entropy image from a set of given images
	* *ex)* `--entropyImage -o output.nrrd input1.nrrd input2.nrrd ...`
* --test
	* Run in a test mode. The test mode is context sensitive depending on the given argument. For example, if `--entropyImage --test` is given, it will automatically provide a set of input images and produce an output into a specific directory.
* --help
	* print this message
