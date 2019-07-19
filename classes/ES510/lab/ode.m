clc, clear all, close all

iex = 5;
%% Example 1: dy = cos(x)
if iex == 1
    c = [-3:1:3]; % a subsetted range of C (in Reality this could be any real value)
    x = [-5:0.1:5]; % pick x range (function actually extends from -inf to + inf)
    
    dy = @(x) cos(x); 
    y = @(x) sin(x);  % KEY
    
    figure; hold on; grid on;
    for ii=1:length(c)
        f = y(x) + c(ii);
        plot(x,f)
    end
    xlabel('x');
    ylabel('y');
    title('dy = cos x');
end

%% Example 2: dy = 2y
if iex == 2
    c = [-5:1:5]; % a subsetted range of C (in Reality this could be any real value)
    x = [-0.5:0.1:1.5]; % pick x range (function actually extends from -inf to + inf)
    
    y = @(x) exp(2*x); % KEY
    
    figure; hold on; grid on;
    for ii=1:length(c)
        f = y(x) + c(ii);
        %if c(ii)==1 
        %    plot(x,f,'LineWidth',2)
        %end
    end
    xlabel('x');
    ylabel('y');
    title('dy = 2y');
end

%% Example 3: Direction field (falling object)
if iex == 3
    t = [0:0.5:5];   % sec
    v0 = [0:5:50];   % m/s
    gravity = 9.8;   % m/s^2
    air_drag = 0.47; % gamma
    arrow_len = 0.5;
    
    figure; hold on
    for ii=1:length(t)
        for jj = 1:length(v0)
            dv = gravity - (air_drag * v0(jj)); % KEY: tangent (derivative)
            % plot
            u = cos(atan(dv)) * arrow_len; % x-component of arrow
            v = sin(atan(dv)) * arrow_len; % y-component of arrow
            quiver(t(ii),v0(jj),u,v,'k'); % plot arrow (XXX: something not right)
        end
    end
    xlabel('t');
    ylabel('v');
    title('Direction field for a falling object; $$\frac{dv}{dt} = g - \frac{\gamma v}{m}$$','interpreter','latex')
end

%% Example 4: Direction field [ y' = cos(t) cos(y) ]
if iex == 4
    x = [0:0.4:2*pi];
    y = [-pi:0.4:pi];
    arrow_len = 0.4;
    
    figure; hold on
    for ii=1:length(x)
        for jj = 1:length(y)
            dy = 2 * cos(x(ii)) * cos(y(jj)); % KEY
            % plot
            u = cos(atan(dy)) * arrow_len; % x-component of arrow
            v = sin(atan(dy)) * arrow_len; % y-component of arrow
            quiver(x(ii),y(jj),u,v,'k'); % plot arrow
        end 
    end
    xlabel('x');
    ylabel('y');
    title('Direction field for $$y'' = 2 cos(x) cos(y)$$','interpreter','latex')
end

%% Example 5: Euler's method [ dy/dx = -y ]
if iex == 5
    % initial condition
    x0 = 0;
    y0 = 1;
    
    % initialize
    x = [0:0.1:2];
    
    % Plot actual solution
    y = @(x) exp(-x); % analytical solution    
    plot(x,y(x),'k','LineWidth',2); % plot analytical solution
    hold on
    
    % Compute numerical solution
    deltax = 0.4; % discretization in x
    xvec = [x0:deltax:x(end)];
    yvec(1) = y0;
    for ii=1:length(xvec)
        yvec(ii+1) = yvec(ii) + (-yvec(ii))*deltax; % KEY
        plot(xvec(ii),yvec(ii),'o','MarkerFaceColor','r','MarkerEdgeColor','k')
    end
    
    deltax = 0.1; % smalled discretization (more accurate results)
    xvec = [x0:deltax:x(end)];
    yvec(1) = y0;
    for ii=1:length(xvec)
        yvec(ii+1) = yvec(ii) + (-yvec(ii))*deltax; % KEY
        plot(xvec(ii),yvec(ii),'o','MarkerFaceColor','b','MarkerEdgeColor','k')
    end
    
    xlabel('x'); ylabel('y');
    title('Euler solution to $\frac{dy}{dx} = -y$','interpreter','latex')
end
