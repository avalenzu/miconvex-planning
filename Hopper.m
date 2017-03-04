classdef Hopper < handle
  properties
    rbm
    rbm_vis
    littleDog
    hopper
    platforms = struct('start', {},  'end', {}, 'height', {}, 'geom', {}) 
    regions
    hip_in_body = struct('front', struct('x', 0.5, 'z', -0.25), ...
                         'hind', struct('x', -0.5, 'z', -0.25))
    leg_length
    v
    t_data     
    r_data     
    v_data     
    F_data     
    r_hip_data 
    p_data     
    f_data     
    f_data_vis
    th_data    
    T_data     
    k_data     
    q_data
    T_actual
    region_indicators
    body_region_indicators
    qtraj
    Ftraj
  end

  methods
    function obj = Hopper(leg_length, hip_in_body)
      urdf = fullfile(getDrakePath(), 'solvers', 'test', 'littleBrick.urdf');
      particle_urdf = fullfile(getDrakePath(), 'systems', 'plants', 'test', 'PointMass.urdf');
      options = struct();
      options.floating = true;
      options.use_new_kinsol = true;
      obj.littleDog = LittleDog();
      obj.rbm = RigidBodyManipulator(urdf, options);
      obj.rbm = obj.rbm.setBody(2, obj.littleDog.getBody(2));
      obj.rbm_vis = obj.rbm;
      obj.leg_length = leg_length;
      obj.hip_in_body = obj.constructHipInBody();

      colormap('lines')
      colors = colormap';
      colormap('default')
      options.collision = false;
      leg = RigidBodyCapsule(0.01, leg_length, [0; 0; leg_length/2], [0; 0; 0]);

      for j = 1:2
        % Add feet
        obj.rbm_vis = obj.rbm_vis.addRobotFromURDF(particle_urdf, [], [], options);
        body = obj.rbm_vis.body(end);
        body.visual_geometry{1} = body.visual_geometry{1}.setColor(colors(:,j));
        body.visual_geometry{1}.radius = 0.02;
        obj.rbm_vis = obj.rbm_vis.setBody(obj.rbm_vis.getNumBodies(), body);
        obj.rbm_vis = obj.rbm_vis.addVisualGeometryToBody(obj.rbm_vis.getNumBodies(), leg);

        %Add hips
        obj.rbm_vis = obj.rbm_vis.addRobotFromURDF(particle_urdf, [], [], options);
        body = obj.rbm_vis.body(end);
        body.visual_geometry{1} = body.visual_geometry{1}.setColor(colors(:,j));
        body.visual_geometry{1}.radius = 0.03;
        obj.rbm_vis = obj.rbm_vis.setBody(obj.rbm_vis.getNumBodies(), body);
      end
      obj.rbm_vis = obj.rbm_vis.compile();
      obj.littleDog = obj.littleDog.removeCollisionGroupsExcept('body');
      obj.littleDog = obj.littleDog.compile();
    end

    function regions = addLittleDogTerrain(obj, terrain_name, num_ccw_rotations, resolution, threshold, min_region_size, visualize)
      checkDependency('iris');
      import iris.terrain_grid.*;
      import iris.inflate_region.*;
      import iris.thirdParty.polytopes.*;
      import iris.cspace.*;
      import iris.drawing.*;

      if nargin < 3 || isempty(num_ccw_rotations)
        num_ccw_rotations = zeros(size(terrain_name));
      elseif iscell(num_ccw_rotations)
        num_ccw_rotations = cell2mat(num_ccw_rotations);
      end
      if nargin < 4 || isempty(resolution), resolution = 0.01; end
      if nargin < 5 || isempty(threshold), threshold = 4e-3; end
      if nargin < 6 || isempty(min_region_size), min_region_size = 50; end
      if nargin < 7 || isempty(visualize), visualize = true; end
      terrain = LittleDogTerrain(terrain_name, num_ccw_rotations);
      obj.littleDog = obj.littleDog.setTerrain(terrain);
      obj.littleDog = obj.littleDog.compile();
      obj.rbm_vis = obj.rbm_vis.setTerrain(terrain);
      obj.rbm_vis = obj.rbm_vis.compile();
      x_min = terrain.x_positions(1);
      x_max = terrain.x_positions(end);
      y_min = 0;
      y_max = 0.5;
      n_x_samples = ceil((x_max - x_min)/resolution);
      n_y_samples = ceil((y_max - y_min)/resolution);
      x = linspace(x_min, x_max, n_x_samples);
      y = linspace(y_min, y_max, n_y_samples);
      [X, Y] = meshgrid(x, y);
      heights = reshape(terrain.getHeight([X(:), Y(:)]'), ...
                        n_y_samples, n_x_samples);
      scale = 1;
      %heights = imresize(heights, scale);

      threshold = threshold/scale;
      Q =     imfilter(heights, [1, -1]) - (threshold) > 0;
      Q = Q | imfilter(heights, [-1, 1]) - (threshold) > 0;
      Q = Q | imfilter(heights, [1; -1]) - (threshold) > 0;
      Q = Q | imfilter(heights, [-1; 1]) - (threshold) > 0;
      Q(isnan(heights)) = 1;

      % grid = Q(85:125,25:85);
      grid = ~Q;
      if visualize
        lcmgl = LCMGLClient('Hopper');
        figure(16)
        %surf(heights);
        imshow(grid, 'InitialMagnification', 'fit')
        hold on;
      end
      t0 = tic;
      [obstacles, A, b, mask] = iris.terrain_grid.segment_grid(grid);
      b = cellfun(@(bi) resolution*bi, b, 'UniformOutput', false);
      toc(t0);
      cull_index = cellfun(@(m) sum(m(:)) < min_region_size, mask);
      obstacles(cull_index) = [];
      A(cull_index) = [];
      b(cull_index) = [];
      mask(cull_index) = [];
      % profile viewer

      x_coords = cell(size(obstacles));
      y_coords = cell(size(obstacles));
      regions = cell(2*numel(obstacles), 1);
      %regions = cell(numel(obstacles), 1);
      for j = 1:length(obstacles)
        obs = obstacles{j};
        if visualize
          figure(16)
          plot(obs(2,:), obs(1,:),'r-');
          lcmgl.glColor3f(0,0,0);
          lcmgl.glLineWidth(5);
          lcmgl.glBegin(lcmgl.LCMGL_LINES);
          for i = 1:size(obs,2)-1
            lcmgl.glVertex3d(resolution*obs(2,i),resolution*obs(1,i),terrain.getHeight([resolution*obs(2,i); resolution*obs(1,i)]));
            lcmgl.glVertex3d(resolution*obs(2,i+1),resolution*obs(1,i+1),terrain.getHeight([resolution*obs(2,i+1); resolution*obs(1,i+1)]));
          end
          lcmgl.glVertex3d(resolution*obs(2,end),resolution*obs(1,end),terrain.getHeight([resolution*obs(2,end); resolution*obs(1,end)]));
          lcmgl.glVertex3d(resolution*obs(2,1),resolution*obs(1,1),terrain.getHeight([resolution*obs(2,1); resolution*obs(1,1)]));
          lcmgl.glEnd();
        end
        [y_coords{j}, x_coords{j}] = find(mask{j});
        zj = heights(sub2ind(size(heights), y_coords{j}, x_coords{j}));
        xj = x(x_coords{j});
        yj = y(y_coords{j});
        f = fit([xj(:), yj(:)], zj, 'poly11');
        % z = p00 +  p10*x + p01*y
        % -p00 == -z + p10*x + p01*y
        % p00/norm([-1; p10; p01]) == normalized([-1; p10; p01])'*[-1; p10; p01]
        n = [f.p10, f.p01, -1];
        c = -f.p00;
        c = c/norm(n);
        n = n./norm(n);
        if n(3) < 0
          n = -n;
          c = -c;
        end
        A{j} = [A{j}(:,2), A{j}(:,1), zeros(size(A{j},1), 1)];
        Aeq{j} = n;
        beq{j} = c;
        regions{j} = struct('A', A{j}, 'b', b{j}./obj.leg_length, 'Aeq', Aeq{j}, 'beq', beq{j}./obj.leg_length, 'mu', 1, 'normal', n');
        regions{numel(obstacles) + j} = struct('A', [A{j}; -n], 'b', [b{j}; -c]./obj.leg_length);
      end
      regions{end+1} = struct('A', [0, 0, -1; -1, 0, 0], 'b', [-max(heights(:)); -x_min]./obj.leg_length);
      if visualize
        lcmgl.switchBuffers();
        hold off
      end
    end

    function q_data = getQData(obj)
      q_data = obj.q_data;
    end

    function constructVisualizer(obj)
      obj.rbm_vis = obj.rbm_vis.compile();
      obj.v = HopperVisualizer(obj.rbm_vis.constructVisualizer());
    end

    function num_positions = getNumPositions(obj)
      num_positions = obj.rbm_vis.getNumPositions();
    end

    function loadResults(obj, data)
      m = obj.littleDog.getMass();
      obj.t_data          = sqrt(obj.leg_length/9.81)*data.t;
      obj.r_data          = obj.leg_length*data.r;
      obj.v_data          = obj.leg_length/sqrt(obj.leg_length/9.81)*data.v;
      obj.F_data          = m*9.81*data.F;
      obj.r_hip_data      = obj.leg_length*data.r_hip;
      obj.p_data          = obj.leg_length*data.p;
      obj.f_data          = data.f;%*obj.littleDog.getMass()*9.81;
      obj.f_data_vis      = data.f*2*obj.leg_length;%*obj.littleDog.getMass()*9.81;
      obj.th_data         = data.th;
      obj.T_data          = data.T*obj.littleDog.getMass()*9.81*obj.leg_length;
      obj.k_data          = data.k*(obj.littleDog.getMass()*obj.leg_length^2)*sqrt(9.81/obj.leg_length);
      obj.region_indicators = round(data.region_indicators);
      obj.body_region_indicators = round(data.body_region_indicators);
      N = size(obj.r_data, 2);
      n_feet = size(obj.p_data, 3);
      leg_pitch_data  = zeros(1, N, double(n_feet));
      for j = 1:double(n_feet)
        leg_pitch_data(:,:,j) = -atan2(obj.p_data(1,:,j), -obj.p_data(3,:,j));
      end
      q_data = zeros(obj.rbm_vis.getNumPositions(), N);
      q_data([1,2,3], :) = obj.r_data;
      q_data(5, :) = obj.th_data;
      for j = 1:n_feet
        q_data(6 + 12*(j-1) + [1,2,3], :) = obj.r_data + obj.r_hip_data(:,:,j) + obj.p_data(:,:,j);
        q_data(6 + 12*(j-1) + 5, :) = leg_pitch_data(:,:,j);
        q_data(6 + 12*(j-1) + [7,8,9], :) = obj.r_data + obj.r_hip_data(:,:,j);
      end
      obj.q_data = q_data;

      t = obj.t_data;

      qtraj = PPTrajectory(foh(t, q_data));
      obj.qtraj = qtraj.setOutputFrame(obj.rbm_vis.getPositionFrame());
      obj.Ftraj = PPTrajectory(zoh(t, reshape(permute(obj.f_data_vis, [1, 3, 2]),[],N)));

      obj.T_actual = -sum((obj.p_data(1,:,:)+obj.r_hip_data(1,:,:)).*obj.f_data(2,:,:) - (obj.p_data(2,:,:)+obj.r_hip_data(2,:,:)).*obj.f_data(1,:,:),3);
    end

    function playback(obj, speed, varargin)
      obj.v.playback_speed = speed;
      obj.v.playback(obj.qtraj, obj.Ftraj, varargin{:})
    end

    function Istar = getDimensionlessMomentOfInertia(obj)
      A = obj.littleDog.centroidalMomentumMatrix(obj.littleDog.home());
      I = A(2,5);
      m = obj.littleDog.getMass();
      Istar = I/(m*obj.leg_length^2);
    end

    function hip_in_body = getHipInBody(obj)
      hip_in_body = obj.hip_in_body;
    end

    function addRegion(obj, region)
      if isempty(obj.regions)
        obj.regions = region;
      else
        obj.regions(end+1) = region;
      end
    end

    function addPlatform(obj, platform_start, platform_end, platform_height, platform_left, platform_right)
      platform.start = platform_start*obj.leg_length;
      platform.end = platform_end*obj.leg_length;
      platform.height = platform_height*obj.leg_length;
      platform.geom = RigidBodyBox([platform.end-platform.start; platform_left - platform_right; 0.1], [(platform.end+platform.start)/2; 0.5*(platform_left+platform_right); platform.height-0.05], [0; 0; 0]);
      obj.platforms(end + 1) = platform;
      obj.rbm_vis = obj.rbm_vis.addVisualGeometryToBody(1, platform.geom);
      obj.rbm_vis = obj.rbm_vis.compile();
      obj.littleDog = obj.littleDog.addGeometryToBody(1, platform.geom);
      obj.littleDog = obj.littleDog.compile();
    end

    function hip_in_body = constructHipInBody(obj)
      for name = {'front', 'back'}
        str = name{1};
        body = obj.littleDog.parseBodyOrFrameID(sprintf('%s_left_hip', str));
        Ttree = obj.littleDog.getBody(body).Ttree;
        hip_in_body.(str).x = Ttree(1,4)/obj.leg_length;
        hip_in_body.(str).z = Ttree(3,4)/obj.leg_length;
      end
    end

    function [sol, prog] = solveCDFKP(obj, seed, options)
      if nargin < 3 || isempty(options), options = struct(); end
      if nargin < 2, seed = []; end
      options = obj.parseOptionsStruct(options); 
      assert(~isempty(obj.r_data), 'You must load data first!');
      robot = obj.littleDog;
      nq = robot.getNumPositions();
      min_distance = 0.03;
      foot = struct('id',[],'in_stance',[]);
      foot(1,1).id = robot.findFrameId('front_left_foot_center');
      foot(1,2).id = robot.findFrameId('front_right_foot_center'); 
      foot(2,1).id = robot.findFrameId('back_left_foot_center');
      foot(2,2).id = robot.findFrameId('back_right_foot_center');
      foot_positions = bsxfun(@plus, obj.r_data, obj.r_hip_data + obj.p_data);
      for i = 1:2
        knee = RigidBodySphere(0.01, [(-1)^i*0.0265; 0; 0], [0; 0; 0]);
        for j = 1:2
          robot = robot.addCollisionGeometryToBody(robot.getFrame(foot(i,j).id).body_ind, knee);
        end
      end
      robot = robot.compile();
      N = size(obj.r_data, 2);

      % Load nominal data
      xstar = home(robot);
      qstar = xstar(1:nq);
      q0 = qstar;
      qf = qstar;

      % Set up time parameters
      dt = diff(obj.t_data);
      dt_range = [min(dt), max(dt)];
      tf_range = [0.5*N*dt_range(1), 3*N*dt_range(2)];
      %tf_range = [obj.t_data(end), obj.t_data(end)];
      
      % Compute q_nom
      q_nom = zeros(nq, N);
      foot_constraints = cell(2,2);
      for n = 1:N
        for i = 1:2 %front-back
          for j = 1:2 %left-right
            lb = NaN(3,1);
            lb([1,3]) = foot_positions([1,3], n, i);
            ub = lb;
            foot_constraints{i, j} = WorldPositionConstraint(robot, foot(i,j).id, zeros(3,1), lb, ub);
          end
        end
        lb = obj.r_data(:, n);
        ub = lb;
        com_constraint = WorldCoMConstraint(robot, lb, ub);
        hip_inds = robot.findPositionIndices('hip_roll');
        posture_constraint = PostureConstraint(robot);
        posture_constraint = posture_constraint.setJointLimits([4;5;6;hip_inds], [0; -pi/8; 0; zeros(size(hip_inds))], [0; pi/8; 0; zeros(size(hip_inds))]);
        min_distance_constraint = MinDistanceConstraint(robot, min_distance);
        ikoptions = IKoptions(robot);
        if n == 1
          qseed = qstar;
        else
          qseed = q_nom(:,n-1);
        end
        [q_nom(:, n), info, infeasible_constraint] = robot.inverseKin(qseed, qstar, foot_constraints{:}, com_constraint, posture_constraint,  ikoptions);
        assert(info < 10)
      end
      %keyboard


      % Set up cost variables
      %q_nom = bsxfun(@times,qstar,ones(1,N));
      %q_nom(5,:) = obj.th_data;
      %q_nom([1,3],:) = obj.r_data;
      state_cost = Point(getStateFrame(robot),ones(getNumStates(robot),1));
      state_cost.base_x = 0;
      state_cost.base_y = 0;
      state_cost.base_roll = 10;
      state_cost.base_yaw = 10;
      state_cost.front_left_hip_roll = 5;
      state_cost.front_right_hip_roll = 5;
      state_cost.back_left_hip_roll = 5;
      state_cost.back_right_hip_roll = 5;
      state_cost = double(state_cost);
      Q = 1e3*diag(state_cost(1:nq)); 
      state_cost(nq+1:nq+6) = 0;
      Qv = 1e3*diag(state_cost(nq+1:end));
      Q_comddot = diag([1,1,1]);
      Q_contact_force = 5*eye(3);

      % Set up common friction-cone variables
      num_edges = 3;
      FC_angles = linspace(0,2*pi,num_edges+1);FC_angles(end) = [];

      % Set up contact wrench structure
      contact_wrench_struct = struct('active_knot', {}, 'cw', {});
      for i = 1:2 % front-back
        for j = 1:numel(obj.regions)
          mu = double(obj.regions(j).mu);
          if mu > 0
            idx = find(obj.region_indicators(j, :, i));
            end_idx = find(diff(idx) ~= 1);
            start_idx = end_idx + 1;
            start_idx = [1, start_idx];
            end_idx = [end_idx,numel(idx)];
            time_index = cell(numel(start_idx),1);
            for l = 1:numel(start_idx)
              time_index{l} = idx(start_idx(l):(end_idx(l)-1));
            end
            if ~isempty(idx)
              for k = 1:2 % left-right
                FC_axis = obj.regions(j).normal;
                FC_perp1 = rotx(pi/2)*FC_axis;
                FC_perp2 = cross(FC_axis, FC_perp1);
                FC_edge = bsxfun(@plus, FC_axis, mu*(bsxfun(@times,cos(FC_angles),FC_perp1) + ...
                  bsxfun(@times,sin(FC_angles),FC_perp2)));
                %contact_wrench_struct(end+1).active_knot = setdiff(idx, _idx);
                contact_wrench_struct(end+1).active_knot = idx;
                contact_wrench_struct(end).active_knot(end_idx) = [];
                %contact_wrench_struct(end+1).active_knot = idx;
                %contact_wrench_struct(end).cw = LinearFrictionConeWrench(robot,foot(i,k).id,zeros(3,1),FC_edge);
                contact_wrench_struct(end).cw = FrictionConeWrench(robot,foot(i,k).id,zeros(3,1),mu,FC_axis);
              end
            end
          end
        end
      end

      options.time_option = 2;
      prog = ComDynamicsFullKinematicsPlanner(robot,N,tf_range,Q_comddot,Qv,Q,q_nom,Q_contact_force,contact_wrench_struct,options);
      %prog = prog.setCheckGrad(true);

      % Add velocity constraints
      max_joint_velocity = 4*pi;
      lb = -max_joint_velocity*ones(nq-6, N);
      ub = max_joint_velocity*ones(nq-6, N);
      prog = prog.addConstraint(BoundingBoxConstraint(lb, ub), prog.v_inds(7:end, :));

      % Add collision avoidance
      prog = prog.addRigidBodyConstraint(MinDistanceConstraint(robot, min_distance),1:N);

      % Add Timestep bounds
      h_min = 0.5*dt_range(1); h_max = 3*dt_range(2);
      prog = prog.addBoundingBoxConstraint(BoundingBoxConstraint(h_min*ones(N-1,1),h_max*ones(N-1,1)),prog.h_inds(:));
      %prog = prog.addConstraint(ConstantConstraint(dt), prog.h_inds(:));
      %prog = prog.addCost(QuadraticConstraint(-Inf, Inf, eye(numel(prog.h_inds)), zeros(numel(prog.h_inds),1)), prog.h_inds(:));

      % Add symmetry constraints
       prog = prog.addConstraint(obj.symmetryConstraint(obj.littleDog, 2:N), prog.q_inds(:,2:end));
      %
      % Constrain hip roll
      posture_constraint = PostureConstraint(robot);
      hip_inds = robot.findPositionIndices('hip_roll');
      knee_inds = robot.findPositionIndices('knee');
      joint_inds = [hip_inds; knee_inds];
      knee_lb = 60*pi/180;
      posture_constraint = posture_constraint.setJointLimits(joint_inds, [zeros(size(hip_inds)); -Inf(2,1); knee_lb*ones(2,1)], [zeros(size(hip_inds)); -knee_lb*ones(2,1); Inf(2,1)]);
      prog = prog.addRigidBodyConstraint(posture_constraint, 1:N);


      % Add initial conditions
      %prog = prog.addConstraint(ConstantConstraint(q_nom(:,1)), prog.q_inds(:,1));
      prog = prog.addConstraint(ConstantConstraint(zeros(3,1)), prog.H_inds(:,1));
      prog = prog.addConstraint(ConstantConstraint(zeros(3,1)), prog.Hdot_inds(:,1));
      prog = prog.addConstraint(ConstantConstraint(obj.r_data(1,1)), prog.com_inds(1,1));
      prog = prog.addConstraint(ConstantConstraint(zeros(3,1)), prog.comdot_inds(:,1));
      %prog = prog.addConstraint(ConstantConstraint(zeros(3,1)), prog.comddot_inds(:,1));
      prog = prog.addConstraint(ConstantConstraint(zeros(nq,1)), prog.v_inds(:,1));

      % Add final conditions
      %prog = prog.addConstraint(ConstantConstraint(qstar(7:end)), prog.q_inds(7:end,N));
      prog = prog.addConstraint(ConstantConstraint(zeros(3,1)), prog.H_inds(:,N));
      prog = prog.addConstraint(ConstantConstraint(zeros(3,1)), prog.Hdot_inds(:,N));
      prog = prog.addConstraint(ConstantConstraint(obj.r_data(1,N)), prog.com_inds(1,N));
      prog = prog.addConstraint(ConstantConstraint(zeros(3,1)), prog.comdot_inds(:,N));
      %prog = prog.addConstraint(ConstantConstraint(zeros(3,1)), prog.comddot_inds(:,N));
      prog = prog.addConstraint(ConstantConstraint(zeros(nq,1)), prog.v_inds(:,N));

      % Add constraints on base
      tol = 0.5;
      %prog = prog.addConstraint(BoundingBoxConstraint(obj.r_data(1,:)-tol, obj.r_data(1,:)+tol), prog.com_inds(1,:));
      %prog = prog.addConstraint(BoundingBoxConstraint(obj.r_data(2,:)-tol, obj.r_data(2,:)+tol), prog.com_inds(3,:));
      %prog = prog.addConstraint(BoundingBoxConstraint(obj.k_data-tol, obj.k_data+tol), prog.H_inds(2,:));
      prog = prog.addConstraint(BoundingBoxConstraint(obj.r_data(2,:), obj.r_data(2,:)), prog.q_inds(2,:));
      prog = prog.addConstraint(BoundingBoxConstraint(zeros(N,1), zeros(N,1)), prog.q_inds(4,:));
      prog = prog.addConstraint(BoundingBoxConstraint(zeros(N,1), zeros(N,1)), prog.q_inds(6,:));
      prog = prog.addConstraint(BoundingBoxConstraint(-pi/8*ones(N,1), pi/8*ones(N,1)), prog.q_inds(5,:));

      % Foot region constraints
      foot_position_fcn = cell(2);
      for i = 1:2
        for k = 1:2
          foot_position_fcn{i, k} = drakeFunction.kinematic.WorldPosition(robot, foot(i,k).id);
          for n = 1:N
            tol = 1e-1;
            lb = -tol*ones(2,1);
            ub = tol*ones(2,1);
            if k == 1
              lb = [lb(1); 0; lb(2)];
              ub = [ub(1); Inf; ub(2)];
            else
              lb = [lb(1); -Inf; lb(2)];
              ub = [ub(1); 0; ub(2)];
            end
            foot_position = foot_positions(:, n, i);
            xz_error_fcn = drakeFunction.Affine([1, 0, 0; 0, 0, 1], -foot_positions([1,3], n, i));
            %xz_error_fcn = drakeFunction.Affine(eye(3), -foot_position);
            norm_squared_fcn = drakeFunction.euclidean.NormSquared(2);
            cost = DrakeFunctionConstraint(-Inf, Inf, norm_squared_fcn(xz_error_fcn(foot_position_fcn{i,k})));
            %constraint = DrakeFunctionConstraint(lb, ub, xz_error_fcn(foot_position_fcn{i,k}));
            cnstr_inds = prog.q_inds(:,n);
            prog = prog.addCost(cost, cnstr_inds);
            %prog = prog.addConstraint(constraint, cnstr_inds);
          end
        end
      end
      for j = 1:numel(obj.regions)
        A = [obj.regions(j).A; obj.regions(j).Aeq; -obj.regions(j).Aeq]./obj.leg_length;
        b = [obj.regions(j).b; obj.regions(j).beq; -obj.regions(j).beq];
        region_fcn = drakeFunction.Affine(A, -b);
        lb = -Inf(size(b));
        ub = zeros(size(b));

        % body region constraint
        idx = find(obj.body_region_indicators(j, :));
        constraint = DrakeFunctionConstraint(lb, ub-0.25, region_fcn);
        for time_index = reshape(idx, 1, [])
          cnstr_inds = prog.q_inds(1:3,time_index);
          prog = prog.addConstraint(constraint,cnstr_inds);
        end

        % leg region constraints
        for i = 1:2 % front-back
          idx = find(obj.region_indicators(j, :, i));
          if ~isempty(idx)
            for k = 1:2 % left-right
              fcn = region_fcn(foot_position_fcn{i,k});
%               debug = drakeFunction.Debug(fcn.dim_output, @(f,df) display(f(f>0)));
%               fcn = compose(debug, fcn);
              constraint = DrakeFunctionConstraint(lb, ub, fcn);
              for time_index = reshape(idx, 1, [])
                cnstr_inds = prog.q_inds(:,time_index);
                prog = prog.addConstraint(constraint,cnstr_inds);
              end
              if obj.regions(j).mu > 0
                % stance feet stationary constraints
                end_idx = find(diff(idx) ~= 1);
                start_idx = end_idx + 1;
                start_idx = [1, start_idx];
                end_idx = [end_idx,numel(idx)];
                time_index = cell(numel(start_idx),1);
                for l = 1:numel(start_idx)
                  time_index{l} = idx(start_idx(l):end_idx(l));
                  if time_index{l}(1) ~= 1
                    time_index{l} = [time_index{l}(1)-1, time_index{l}];
                  end
                end
                prog = prog.addRigidBodyConstraint(WorldFixedPositionConstraint(robot,foot(i,k).id, zeros(3,1)),time_index);
              else
                % swing feet collision avoidance
                swing_idx = setdiff([idx-1, idx+1], idx);
                swing_idx(swing_idx < 1 | swing_idx > N) = [];
                for time_index = reshape(swing_idx, 1, [])
                  cnstr_inds = prog.q_inds(:,time_index);
                  prog = prog.addConstraint(constraint,cnstr_inds);
                end
              end
            end
          end
        end
      end

      % Set up seed
      if isempty(seed)
        x_seed = zeros(prog.num_vars,1);
        q_seed = q_nom;
        %q_seed([1, 3], :) = obj.r_data;
        %q_seed(5, :) = obj.th_data;
        v_seed = gradient(q_seed);
        com_seed = obj.r_data;
        comdot_seed = zeros(3, N);
        comdot_seed([1, 3], :) = obj.v_data;
        comddot_seed = zeros(3, N);
        comddot_seed([1, 3], :) = (1/obj.littleDog.getMass())*obj.F_data;
        I = obj.getDimensionlessMomentOfInertia()*obj.littleDog.getMass()*obj.leg_length^2;
        H_seed = zeros(3, N);
        H_seed(2,:) = obj.k_data;
        Hdot_seed = gradient(H_seed);
        x_seed(prog.H_inds(:)) = reshape(H_seed,[],1);
        x_seed(prog.Hdot_inds(:)) = reshape(Hdot_seed,[],1);
        x_seed(prog.h_inds) = dt;
        x_seed(prog.q_inds(:)) = reshape(q_seed,[],1);
        x_seed(prog.v_inds(:)) = reshape(v_seed,[],1);
        x_seed(prog.com_inds(:)) = reshape(com_seed,[],1);
        x_seed(prog.comdot_inds(:)) = reshape(comdot_seed,[],1);
        x_seed(prog.comddot_inds(:)) = reshape(comddot_seed,[],1);
        for i = 1:2
          for j = 1:2
            f_seed = zeros(3, N);
            f_seed([1, 3], :) = obj.f_data(:, :, i);
            x_seed(prog.lambda_inds{2*(i-1) + j}(:)) = 0.5*f_seed;
            %prog = prog.addCost(QuadraticConstraint(-Inf, Inf, eye(numel(f_seed)), -2*f_seed(:)), prog.lambda_inds{2*(i-1) + j}(:));
          end
        end
      else
        x_seed = seed.x_sol;
      end

      % Set up solver options
      %prog = prog.setSolver('ipopt');
      prog = prog.setSolverOptions('snopt','iterationslimit',1e6);
      prog = prog.setSolverOptions('snopt','majoriterationslimit', options.major_iteration_limit);
      prog = prog.setSolverOptions('snopt','majorfeasibilitytolerance',1e-5);
      prog = prog.setSolverOptions('snopt','majoroptimalitytolerance',6e-4);
      prog = prog.setSolverOptions('snopt','superbasicslimit',2000);
      prog = prog.setSolverOptions('snopt','linesearchtolerance',0.99);
      prog = prog.setSolverOptions('snopt','scaleoption',1);
      %prog = prog.setSolverOptions('snopt','sense','Feasible point');
      prog = prog.setSolverOptions('snopt','print','snopt.out');

      % Solve trajectory optimization
      tic
      %profile on;
      [x_sol,~,~] = prog.solve(x_seed);
      %profile off;
      toc

      % Parse trajectory optimization output
      sol.x_sol = x_sol;
      sol.q = reshape(x_sol(prog.q_inds(:)),nq,N);
      sol.v = reshape(x_sol(prog.v_inds(:)),nq,N);
      sol.h = reshape(x_sol(prog.h_inds),1,[]);
      sol.t = cumsum([0 sol.h]);
      sol.com = reshape(x_sol(prog.com_inds),3,[]);
      sol.comdot = reshape(x_sol(prog.comdot_inds),3,[]);
      sol.comddot = reshape(x_sol(prog.comddot_inds),3,[]);
      sol.H = reshape(x_sol(prog.H_inds),3,[]);
      sol.Hdot = reshape(x_sol(prog.Hdot_inds),3,[])*prog.torque_multiplier;
      sol.lambda = cell(2,1);
      for i = 1:numel(prog.lambda_inds)
        sol.lambda{i} = reshape(x_sol(prog.lambda_inds{i}),size(prog.lambda_inds{i},1),[],N);
      end
      sol.xtraj= PPTrajectory(foh(sol.t,[sol.q;sol.v]));
      sol.xtraj= sol.xtraj.setOutputFrame(robot.getStateFrame);
    end

    function delete(obj)
      obj.rbm = [];
      obj.rbm_vis = [];
      obj.littleDog = [];
      obj.hopper = [];
      obj.platforms = [];
      obj.hip_in_body = [];
      obj.leg_length = [];
      obj.v = [];
      obj.t_data = [];
      obj.r_data = [];
      obj.r_hip_data = [];
      obj.p_data = [];
      obj.f_data = [];
      obj.th_data = [];
      obj.T_data = [];
      obj.q_data = [];
      obj.T_actual = [];
      obj.region_indicators = [];
      obj.regions = [];
      obj.qtraj = [];
      obj.Ftraj = [];
    end
  end

  methods (Static)
    function options = parseOptionsStruct(options_in)
      options = Hopper.defaultOptionsStruct();
      for fieldname_cell = fields(options_in)'
        fieldname = fieldname_cell{1};
        if isfield(options,fieldname)
          options.(fieldname) = options_in.(fieldname);
        end
      end
    end

    function options = defaultOptionsStruct()
      options.visualize = true;
      options.major_iteration_limit = 200;
      options.suffix = '';
    end

    function symmetry_cnstr = symmetryConstraint(robot,t_idx)
      N = numel(t_idx);
      num_symmetry = 6;
      symmetry_matrix = zeros(num_symmetry,robot.getNumPositions);
      function mat = addSymmetricPair(mat, rows, joint1_name, joint2_name)
        joint1_idx = robot.getBody(robot.findJointId(joint1_name)).position_num;
        joint2_idx = robot.getBody(robot.findJointId(joint2_name)).position_num;
        mat(rows, [joint1_idx, joint2_idx]) = [-1, 1];
      end
      function mat = addAntiSymmetricPair(mat, rows, joint1_name, joint2_name)
        joint1_idx = robot.getBody(robot.findJointId(joint1_name)).position_num;
        joint2_idx = robot.getBody(robot.findJointId(joint2_name)).position_num;
        mat(rows, [joint1_idx, joint2_idx]) = [1, 1];
      end
      symmetry_matrix = addAntiSymmetricPair(symmetry_matrix, 1, ...
                                              'front_left_hip_roll', ...
                                              'front_right_hip_roll');
      symmetry_matrix = addSymmetricPair(symmetry_matrix, 2, ...
                                              'front_left_hip_pitch', ...
                                              'front_right_hip_pitch');
      symmetry_matrix = addSymmetricPair(symmetry_matrix, 3, ...
                                              'front_left_knee', ...
                                              'front_right_knee');
      symmetry_matrix = addAntiSymmetricPair(symmetry_matrix, 4, ...
                                              'back_left_hip_roll', ...
                                              'back_right_hip_roll');
      symmetry_matrix = addSymmetricPair(symmetry_matrix, 5, ...
                                              'back_left_hip_pitch', ...
                                              'back_right_hip_pitch');
      symmetry_matrix = addSymmetricPair(symmetry_matrix, 6, ...
                                              'back_left_knee', ...
                                              'back_right_knee');
      symmetry_cnstr = LinearConstraint(zeros(num_symmetry*N,1),zeros(num_symmetry*N,1),kron(speye(N),symmetry_matrix));
      cnstr_name = cell(num_symmetry*N,1);
      for i = 1:N
        cnstr_name{(i-1)*num_symmetry+1} = sprintf('front_hip_roll_symmetry[%d]',t_idx(i));
        cnstr_name{(i-1)*num_symmetry+2} = sprintf('front_hip_pitch_symmetry[%d]',t_idx(i));
        cnstr_name{(i-1)*num_symmetry+3} = sprintf('front_knee_symmetry[%d]',t_idx(i));
        cnstr_name{(i-1)*num_symmetry+4} = sprintf('back_hip_roll_symmetry[%d]',t_idx(i));
        cnstr_name{(i-1)*num_symmetry+5} = sprintf('back_hip_pitch_symmetry[%d]',t_idx(i));
        cnstr_name{(i-1)*num_symmetry+6} = sprintf('back_knee_symmetry[%d]',t_idx(i));
      end
      symmetry_cnstr = symmetry_cnstr.setName(cnstr_name);
    end
  end
end
