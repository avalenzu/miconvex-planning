function d = nparray2double2(nparray)
d = reshape(double(py.array.array('d', nparray.flat)), double(py.array.array('d', nparray.shape)));