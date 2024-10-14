import vtkreader_merope
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN
import numpy as np


fileName = "composite_normals.vtk"
fileName_ref = "ref_composite_normals.vtk"
normalName = "composite_normals"

def read_normals(fileName):
	all_data = vtkreader_merope.read_file(fileName)
	data = VN.vtk_to_numpy(all_data.GetCellData().GetArray(normalName))
	return data


def difference(data1, data2):
	res = 0
	for i in range(len(data1)):
		vect = data1[i] -  data2[i]
		res += np.dot(vect, vect)
	return res


if difference(read_normals(fileName), read_normals(fileName_ref)) > 1.:
	raise Exception("Impossible")
