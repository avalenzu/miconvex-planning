function d = nparray2double3(nparray)
d = reshape(cell2mat(cell(nparray.ravel.tolist())), cell2mat(cell(nparray.shape)));