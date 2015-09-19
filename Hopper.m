classdef Hopper < handle
  properties
    rbm
    rbm_vis
    platforms = struct('start', {},  'end', {}, 'height', {}, 'geom', {}) 
    hip_in_body = struct('front', struct('x', 0.5, 'z', -0.25), ...
                         'hind', struct('x', -0.5, 'z', -0.25))
    leg_length
    v
    t_data     
    r_data     
    r_hip_data 
    p_data     
    f_data     
    th_data    
    T_data     
    q_data
    T_actual
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
      obj.rbm = RigidBodyManipulator(urdf, options);
      obj.rbm_vis = obj.rbm;
      obj.hip_in_body = hip_in_body;
      obj.leg_length = leg_length;

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
      obj.t_data          = sqrt(obj.leg_length/9.81)*data.t;
      obj.r_data          = obj.leg_length*data.r;
      obj.r_hip_data      = obj.leg_length*data.r_hip;
      obj.p_data          = obj.leg_length*data.p;
      obj.f_data          = data.f;%*obj.rbm_vis.getMass()*9.81;
      obj.th_data         = data.th;
      obj.T_data          = data.T*obj.rbm_vis.getMass()*9.81*obj.leg_length;
      N = size(obj.r_data, 2);
      n_feet = size(obj.p_data, 3);
      leg_pitch_data  = zeros(1, N, double(n_feet));
      for j = 1:double(n_feet)
        leg_pitch_data(:,:,j) = -atan2(obj.p_data(1,:,j), -obj.p_data(2,:,j));
      end
      q_data = zeros(obj.rbm_vis.getNumPositions(), N);
      q_data([1,3], :) = obj.r_data;
      q_data(5, :) = obj.th_data;
      for j = 1:n_feet
        q_data(6 + 12*(j-1) + [1,3], :) = obj.r_data + obj.r_hip_data(:,:,j) + obj.p_data(:,:,j);
        q_data(6 + 12*(j-1) + 5, :) = leg_pitch_data(:,:,j);
        q_data(6 + 12*(j-1) + [7,9], :) = obj.r_data + obj.r_hip_data(:,:,j);
      end
      obj.q_data = q_data;

      t = obj.t_data;

      qtraj = PPTrajectory(foh(t, q_data));
      obj.qtraj = qtraj.setOutputFrame(obj.rbm_vis.getPositionFrame());
      obj.Ftraj = PPTrajectory(foh(t, reshape(permute(obj.f_data, [1, 3, 2]),[],N)));

      obj.T_actual = sum((obj.p_data(1,:,:)+obj.r_hip_data(1,:,:)).*obj.f_data(2,:,:) - (obj.p_data(2,:,:)+obj.r_hip_data(2,:,:)).*obj.f_data(1,:,:),3);
    end

    function playback(obj, speed)
      obj.v.playback_speed = speed;
      obj.v.playback(obj.qtraj, obj.Ftraj)
    end

    function Istar = getDimensionlessMomentOfInertia(obj)
      m = obj.rbm.getMass();
      I = obj.rbm.body(2).inertia(2,2);
      Istar = I/(m*obj.leg_length^2);
    end

    function hip_in_body = getHipInBody(obj)
      hip_in_body = obj.hip_in_body;
    end

    function addPlatform(obj, platform_start, platform_end, platform_height)
      platform.start = platform_start*obj.leg_length;
      platform.end = platform_end*obj.leg_length;
      platform.height = platform_height*obj.leg_length;
      platform.geom = RigidBodyBox([platform.end-platform.start; 1; 0.1], [(platform.end+platform.start)/2; 0; platform.height-0.05], [0; 0; 0]);
      obj.platforms(end + 1) = platform;
      obj.rbm_vis = obj.rbm_vis.addVisualGeometryToBody(1, platform.geom);
      obj.rbm_vis = obj.rbm_vis.compile();
    end
  end
end
