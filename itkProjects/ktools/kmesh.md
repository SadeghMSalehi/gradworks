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
	* *ex)* -vti imageFile
* -h
	* print help message
