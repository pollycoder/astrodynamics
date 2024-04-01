function [best_params, best_value] = pso_parallel(obj_function, param_ranges, num_particles, max_iter)
    % obj_function: 要优化的目标函数句柄，接受参数向量作为输入，返回一个标量数值
    % param_ranges: 参数范围矩阵，每行代表一个参数的范围 [min_value, max_value]
    % num_particles: 粒子数量
    % max_iter: 最大迭代次数
    if nargin < 4
        max_iter = 1000;
    end
    if nargin < 3
        num_particles = 1000;
    end

    num_params = size(param_ranges, 1);
    
    % 初始化粒子的位置和速度
    particles_pos = rand(num_particles, num_params);
    for i = 1:num_params
        particles_pos(:, i) = particles_pos(:, i) * (param_ranges(i, 2) - param_ranges(i, 1)) + param_ranges(i, 1);
    end

    particles_vel = zeros(num_particles, num_params);

    % 初始化粒子的最佳位置和最佳适应度值
    particles_best_pos = particles_pos;
    particles_best_value = zeros(num_particles, 1);

    parfor i = 1:num_particles
        particles_best_value(i) = obj_function(particles_best_pos(i, :));
    end

    % 全局最佳位置和最佳适应度值
    [global_best_value, global_best_index] = min(particles_best_value);
    global_best_pos = particles_best_pos(global_best_index, :);

    % PSO 参数
    w = 0.5; % 惯性权重
    c1 = 2.0; % 个体学习因子
    c2 = 2.0; % 社会学习因子

    % PSO 主循环
    for iter = 1:max_iter
        parfor i = 1:num_particles
            % 更新速度
            r1 = rand(1, num_params);
            r2 = rand(1, num_params);
            particles_vel(i, :) = w * particles_vel(i, :) + c1 * r1 .* (particles_best_pos(i, :) - particles_pos(i, :)) + c2 * r2 .* (global_best_pos - particles_pos(i, :));

            % 限制速度在参数范围内
            particles_vel(i, :) = min(max(particles_vel(i, :), param_ranges(:, 1)'), param_ranges(:, 2)');

            % 更新位置
            particles_pos(i, :) = particles_pos(i, :) + particles_vel(i, :);

            % 限制位置在参数范围内
            particles_pos(i, :) = min(max(particles_pos(i, :), param_ranges(:, 1)'), param_ranges(:, 2)');

            % 计算适应度值
            value = obj_function(particles_pos(i, :));

            % 更新个体最佳位置和最佳适应度值
            if value < particles_best_value(i)
                particles_best_value(i) = value;
                particles_best_pos(i, :) = particles_pos(i, :);
            end
        end

        % 更新全局最佳位置和最佳适应度值
        [min_value, min_index] = min(particles_best_value);
        if min_value < global_best_value
            global_best_value = min_value;
            global_best_pos = particles_best_pos(min_index, :);
        end
        if mod(iter,100)==0
            disp(['Iteration ', num2str(iter), ', Best Value: ', num2str(global_best_value)]);
        end
    end

    % 返回优化结果
    best_params = global_best_pos;
    best_value = global_best_value;
end