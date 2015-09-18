if ~exist('vis_only', 'var'), vis_only = false; end
if ~vis_only
import py.numpy.array
import py.numpy.matrix
import py.list

np.array = @(mat) py.numpy.array(list(mat));
np.matrix = @(str) py.numpy.matrix(str);
dict = @(varargin) py.dict(pyargs(varargin{:}));
timeout = 40;
threads = 8;
N = 20;
tf = 2*1.6;
leg_length = 0.3;
dt = tf/N/sqrt(leg_length/9.81);
%dt = 0.08/sqrt(leg_length/9.81);
r0 = [0; leg_length/2];
th0 = 0.0;
rf = [2; leg_length];
v0 =  [0; 0];
w0 = 0;
hip_offset_x = 0.5;
hip_in_body = [[-hip_offset_x; -0.25], [hip_offset_x; -0.25]];
% hip_in_body = [-hip_offset_x ; -0.25];
%%
options = struct();
options.floating = true;
options.use_new_kinsol = true;
urdf = fullfile(getDrakePath(), 'solvers', 'test', 'littleBrick.urdf');
particle_urdf = fullfile(getDrakePath(), 'systems', 'plants', 'test', 'PointMass.urdf');
rbm = RigidBodyManipulator(urdf, options);
rbm_vis = rbm;

step_height = 0.1;
platform1_start = -0.5;
platform1_end = 0.3;
platform1_height = 0*step_height;
platform2_start = 0.7;
platform2_end = 1.3;
platform2_pitch = 0*pi/180;
platform2_height = step_height;
platform3_start = 1.7;
% platform3_start = -0.5;
platform3_end = 9.5;
platform3_height = 2*step_height;
platform1 = RigidBodyBox([platform1_end-platform1_start; 1; 0.1], [(platform1_end+platform1_start)/2; 0; platform1_height-0.05], [0; 0; 0]);
platform1.c = [0.2; 0.2; 0.2];
platform2 = RigidBodyBox([platform2_end-platform2_start; 1; 0.1], [(platform2_end+platform2_start)/2; 0; platform2_height-0.05], [0; platform2_pitch; 0]);
platform2.c = [0.2; 0.2; 0.2];
platform3 = RigidBodyBox([platform3_end-platform3_start; 1; 0.1], [(platform3_end+platform3_start)/2; 0; platform3_height-0.05], [0; 0; 0]);
platform3.c = [0.2; 0.2; 0.2];
rbm_vis = rbm_vis.addVisualGeometryToBody(1, platform1);
rbm_vis = rbm_vis.addVisualGeometryToBody(1, platform2);
rbm_vis = rbm_vis.addVisualGeometryToBody(1, platform3);

colormap('lines')
colors = colormap';
colormap('default')
options.collision = false;
leg = RigidBodyCapsule(0.01, leg_length, [0; 0; leg_length/2], [0; 0; 0]);

for j = 1:size(hip_in_body,2)
  % Add feet
  rbm_vis = rbm_vis.addRobotFromURDF(particle_urdf, [], [], options);
  body = rbm_vis.body(end);
  body.visual_geometry{1} = body.visual_geometry{1}.setColor(colors(:,j));
  body.visual_geometry{1}.radius = 0.02;
  rbm_vis = rbm_vis.setBody(rbm_vis.getNumBodies(), body);
  rbm_vis = rbm_vis.addVisualGeometryToBody(rbm_vis.getNumBodies(), leg);

  %Add hips
  rbm_vis = rbm_vis.addRobotFromURDF(particle_urdf, [], [], options);
  body = rbm_vis.body(end);
  body.visual_geometry{1} = body.visual_geometry{1}.setColor(colors(:,j));
  body.visual_geometry{1}.radius = 0.03;
  rbm_vis = rbm_vis.setBody(rbm_vis.getNumBodies(), body);
end
rbm_vis = rbm_vis.compile();
v = HopperVisualizer(rbm_vis.constructVisualizer());
m = rbm.getMass();
I = rbm.body(2).inertia(2,2);
Istar = I/(m*leg_length^2);
%%
hopper = py.hopper.Hopper(N);
hopper.mdt_precision = int32(2);
hopper.dt = dt;
hopper.momentOfInertia = Istar;
hopper.nOrientationSectors = int32(1);
hopper.rotationMax = pi/8;
hopper.velocityMax = 3;
hopper.positionMax = 10;
hopper.forceMax = 3;
hopper.footnames = list({'front', 'hind'});
hopper.hipOffset = dict('front', dict('x', hip_offset_x, 'z', -0.25), ...
                        'hind', dict('x', -hip_offset_x, 'z', -0.25));
                                
hopper.addRegion(pyargs('A', np.matrix('0., -1.,; 1., 0.'), ...
                        'b', np.matrix(sprintf('%f; %f', ...
                                               [-platform1_height, platform2_start]./leg_length))));
hopper.addRegion(pyargs('A', np.matrix('0., -1.,; -1., 0.; 1., 0.'), ...
  'b', np.matrix(sprintf('%f; %f; %f', ...
  [-platform2_height, -platform1_end, max(platform2_end, platform3_start)]./leg_length))));
hopper.addRegion(pyargs('A', np.array([0, -1]), 'b', -platform3_height/leg_length));

hopper.addRegion(pyargs('A', np.matrix('-1., 0.,; 1., 0.'), ...
                        'b', np.matrix(sprintf('%f; %f', ...
                                               [-platform1_start, platform1_end]./leg_length)), ...
                        'Aeq', np.array([0, 1]), ...
                        'beq', platform1_height/leg_length, ...
                        'normal', np.matrix('0.; 1.'), 'mu', 1));
hopper.addRegion(pyargs('A', np.matrix('-1., 0.,; 1., 0.'), ...
                        'b', np.matrix(sprintf('%f; %f', ...
                                               [-platform2_start, platform2_end]./leg_length)), ...
                        'Aeq', np.array([0, 1]), ...
                        'beq', platform2_height/leg_length, ...
                        'normal', np.matrix('0.; 1.'), 'mu', 1));
hopper.addRegion(pyargs('A', np.matrix('-1., 0.,; 1., 0.'), ...
                        'b', np.matrix(sprintf('%f; %f', ...
                                               [-platform3_start, platform3_end]./leg_length)), ...
                        'Aeq', np.array([0, 1]), ...
                        'beq', platform3_height/leg_length, ...
                        'normal', np.matrix('0.; 1.'), 'mu', 1));

%%
m_nlp = py.hopper.testHopper(hopper, list(r0'), list(rf'), leg_length, pyargs('timeout', timeout, 'threads',threads)); 
m = py.hopperUtil.constructMDTModel(m_nlp, hopper.mdt_precision);

gurobi_options = pyargs('TimeLimit', timeout, 'Threads', threads);
opt = py.hopperUtil.constructGurobiSolver(gurobi_options);
opt_nlp = py.hopperUtil.constructCouenneSolver();
%%
results = opt.solve(m, pyargs('tee', true));

%%
results = opt_nlp.solve(m_nlp, pyargs('tee', true));

end
%%
r_data = leg_length*nparray2double(py.hopperUtil.extractPostition(m));
r_hip_data = leg_length*nparray2double(py.hopperUtil.extractHipPosition(m));
p_data = leg_length*nparray2double(py.hopperUtil.extractRelativeFootPosition(m));
f_data = nparray2double(py.hopperUtil.extractFootForce(m));
th_data = nparray2double(py.hopperUtil.extractOrientation(m));
T_data = nparray2double(py.hopperUtil.extractTotalTorque(m));
leg_pitch_data = zeros(1, hopper.N, double(py.len(hopper.footnames)));
% force_angle_data = zeros(1, hopper.N, double(py.len(hopper.footnames)));
for j = 1:double(py.len(hopper.footnames))
  leg_pitch_data(:,:,j) = -atan2(p_data(1,:,j), -p_data(2,:,j));
%   force_angle_data(:,:,j) = atan2(prog_prev.vars.F.value(1,:,j), prog_prev.vars.F.value(2,:,j));
  %force_angle_data(sqrt(sum(prog_prev.vars.F.value(:,:,j).^2)) < 1e-6) = 0;
end

q_data = zeros(rbm_vis.getNumPositions(), hopper.N);
q_data([1,3], :) = r_data;
q_data(5, :) = th_data;
for j = 1:double(py.len(hopper.footnames))
  q_data(6 + 12*(j-1) + [1,3], :) = r_data + r_hip_data(:,:,j) + p_data(:,:,j);
  %q_data([7,9], :) = r_data + p_data;
  q_data(6 + 12*(j-1) + 5, :) = leg_pitch_data(:,:,j);
  q_data(6 + 12*(j-1) + [7,9], :) = r_data + r_hip_data(:,:,j);
end

t = sqrt(leg_length/9.81)*(0:dt:(N-1)*dt);

qtraj = PPTrajectory(foh(t, q_data));
qtraj = qtraj.setOutputFrame(rbm_vis.getPositionFrame());
Ftraj = PPTrajectory(foh(t, reshape(permute(f_data, [1, 3, 2]),[],hopper.N)));

T_actual = sum((p_data(1,:,:)+r_hip_data(1,:,:)).*f_data(2,:,:) - (p_data(2,:,:)+r_hip_data(2,:,:)).*f_data(1,:,:),3);
