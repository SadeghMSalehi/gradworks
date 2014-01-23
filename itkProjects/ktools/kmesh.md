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
	* *ex)* -vti imageFile outputFile [-attrDim 3]
* -vtu
	* Convert an ITK image to VTU format (vtkUnstructuredGrid). This is useful when masking is needed.
	* *ex)* -vtu imageFile outputFile -maskImage maskImage
* -maskImage
	* A mask image for the use of -vtu
	* *ex)* -maskImage mask.nrrd
* -traceStream
	* Trace a stream line from a given point set
	* *ex)* -traceStream input-vtu-field input-vtk output-vtu
* -h
	* print help message
