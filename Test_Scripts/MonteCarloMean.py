import vtk
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN
import numpy as np

N_R = 0
endtimestep = 0
for i in range(endtimestep):
    for j in range(N_R):
        if(i < 10):
            filename = "Realization_Nr_"+str(j)+".0000"+str(i)+".vtk"
        elif(i >= 10 and i < 100):
            filename = "Realization_Nr_"+str(j)+".000"+str(i)+".vtk"
        elif(i >= 100):
            filename = "Realization_Nr_"+str(j)+".00"+str(i)+".vtk"

        reader = vtk.vtkGenericDataObjectReader()
        reader.SetFileName(filename)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        p = np.array(data.GetPointData().GetArray('p'))
        u1 = np.array(data.GetPointData().GetArray('u1'))
        u2 = np.array(data.GetPointData().GetArray('u2'))
        
        