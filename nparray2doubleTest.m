function nparray2doubleTest()
N = 1e4;
size = [3, 50];
% generate the numpy arrays
numpyArrays = cell(N, 1);
for i = 1:N
  numpyArrays{i} = py.numpy.random.random(py.tuple(size));
end

% convert
timer1 = tic;
for i = 1:N
  d = nparray2double(numpyArrays{i});
end
time1 = toc(timer1);

timer2 = tic;
for i = 1:N
  d = nparray2double2(numpyArrays{i});
end
time2 = toc(timer2);

timer3 = tic;
for i = 1:N
  d = nparray2double3(numpyArrays{i});
end
time3 = toc(timer3);

fprintf('Method 1: %8.5f s\n', time1)
fprintf('Method 2: %8.5f s\n', time2)
fprintf('Method 3: %8.5f s\n', time3)
