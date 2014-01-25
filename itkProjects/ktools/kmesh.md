## *kmesh* Usage
* -exportScalars
	* Export scalar values to a text file
	* *ex)* -exportScalars [in-mesh] [scalar.txt]
* -importScalars
	* Add scalar values to a mesh [in-mesh] [scalar.txt] [out-mesh]
* -smoothScalars
	* Gaussian smoothing of scalar values of a mesh. [in-mesh] [out-mesh]
* -appendData
	* Append input meshes into a single data [output-mesh]
* -sigma
	* sigma value [double]
* -scalarName
	* scalar name [string]
* -outputScalarName
	* scalar name for output [string]
* -iter
	* number of iterations [int]
* -attrDim
	* The number of components of attribute
	* *ex)* -attrDim 3 (vector)
* -vti
	* Convert an ITK image to VTI format (VTKImageData)
	* *ex)* -vti imageFile outputFile [-attrDim 3] [-maskImage mask]
* -vtu
	* Convert an ITK image to VTU format (vtkUnstructuredGrid). This is useful when masking is needed.
	* *ex)* -vtu imageFile outputFile -maskImage maskImage
* -maskImage
	* A mask image for the use of -vtu
	* *ex)* -maskImage mask.nrrd
* -traceStream
	* Trace a stream line from a given point set
	* *ex)* -traceStream input-vtu-field input-vtk output-lines output-points
* -traceDirection
	* Choose the direction of stream tracing (both, forward, backward)
	* *ex)* -traceStream ... -traceDirection (both|forward|backward)
* -zrotate
	* Rotate all the points along the z-axis. Change the sign of x and y coordinate.
	* *ex)* -traceStream ... -zrotate
* -filterStream
	* Filter out stream lines which are lower than a given threshold
	* *ex)* -filterStream stream-line-input stream-seed-input stream-line-output -scalarName scalar -threshold xx
* -thresholdMin
	* Give a minimum threshold value for -filterStream
	* *ex)* -threshold 10 (select a cell whose attriubte is greater than 10)
* -thresholdMax
	* Give a maximum threshold value for -filterStream
	* *ex)* -threshold 10 (select a cell whose attriubte is lower than 10)
* -fitting
	* Fit a model into a binary image
	* *ex)* -fitting input-model binary-image output-model
* -ellipse
	* Create an ellipse with parameters []
	* *ex)* -ellipse 101 101 101 51 51 51 20 20 20 -o ellipse.nrrd
* -o
	* Specify a filename for output; used with other options
	* *ex)* -o filename.nrrd
* -h
	* print help message
