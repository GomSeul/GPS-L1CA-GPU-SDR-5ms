function plot_3d_reftrj(llh,ar,hori_unit,vert_unit,varargin)
%   Plot 3D reference trajectory
%
%	plot_3d_reftrj(llh,ar,hori_unit,vert_unit,...)
%
%   INPUTS
%       llh = latitude (deg), longitude (deg), height (m)
%       ar = axis ratio (1: equal, 2: tight)
%       hori_unit = unit of horizontal axis
%              ('m': meter, 'km': kilo-meter)
%       vert_unit = unit of vertical axis
%              ('m': meter, 'km': kilo-meter, 'ft': feet)
%       other arguments = same as plot function
%
%   OUTPUTS
%       none
%
%

pos_enu = convert_position_llh2enu(llh,llh(1,:));
pos_enu = [pos_enu(:,1:2),llh(:,3)];

switch hori_unit
    case 'm'
        plot_trj(pos_enu,vert_unit,varargin{:});
        xlabel('East (m)');
        ylabel('North (m)');
        grid;
        if ar == 1
            axis equal;
        else
            axis tight;
        end
    case 'km'
        pos_enu(:,1:2) = pos_enu(:,1:2) * 1e-3;
        plot_trj(pos_enu,vert_unit,varargin{:});
        xlabel('East (km)');
        ylabel('North (km)');
        grid;
        if ar == 1
            axis equal;
        else
            axis tight;
        end
    otherwise
        error('Not supported horizontal axis unit');
end
end

function plot_trj(pos_enu,vert_unit,varargin)
figure;
switch vert_unit
    case 'm'
        plot3(pos_enu(:,1),pos_enu(:,2),pos_enu(:,3),varargin{:});
%         plot(pos_enu(:,1),pos_enu(:,2),pos_enu(1,1),pos_enu(1,2),'o',pos_enu(end,1),pos_enu(end,2),'o');       %
        zlabel('Altitude (m)');
    case 'km'
        plot3(pos_enu(:,1),pos_enu(:,2),pos_enu(:,3) * 1e-3,varargin{:});
        zlabel('Altitude (km)');
    case 'ft'
        plot3(pos_enu(:,1),pos_enu(:,2),pos_enu(:,3) * 3.2808,varargin{:});
        zlabel('Altitude (ft)');
    otherwise
        error('Not supported vertical axis unit');
end
end
