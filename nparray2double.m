function d = nparray2double(nparray)
new_shape = cell2mat(cell(nparray.shape));
if numel(new_shape) < 2
  new_shape(2) = 1;
end
d = permute(reshape(double(py.array.array('d', nparray.flat)), fliplr(new_shape)), numel(new_shape):-1:1);