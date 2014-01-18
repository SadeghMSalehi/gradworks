## kcalc usage 
   kcalc [-e equation] [-o output-file] input1:A input2:B ...

The kcalc performs a pixel-wise arithmetic. The pixel value of each input image is given as variables, A,B,C,D, and E. Several functions implemented in [MuParser](http://muparser.beltoforion.de/) includes +,-,*,/, and ?, as well as, trigonometric functions.

Also, there are the min, max values of each input image for the use of scaling and other purposes, which are given as AMIN, AMAX, BMIN, BMAX, and etc.

Note that the output data type is the same with the last input file. The order of images may produce different results, if images with different types are used.

Some examples are:
* **Addition**: -e (A+B)
* **Averaging**: -e (A+B/2)
* **Thresholding**: -e (A>10?1:0)
* **Scaling**: -e (A-AMIN)/AMAX*255
* **Masking**: -e (A==8?B:0)
* ...

### Options
* -e
	* The equation to compute each output pixel.
	* *ex)* -e (A+B)
* -o
	* output filename (the same data type with the last input)
	* *ex)* -o output.nrrd
* -h
	* print this message

