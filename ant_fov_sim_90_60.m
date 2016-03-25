function [] = ant_fov_sim_90_60(inc,alt,dl_loc,launch_date,n)
% title = ANT_FOV_SIM_(A)_(B).m
% A = antenna boresight w.r.t. FIDI boresight
% B = +/- antenna FOV
% simulates a circular orbit and finds the points at which the satellite is
% over a specified ground station (phoenix) when launched at 12:00 (noon)
% from cape canaveral, florida
% info = [phi,eclipse,t]
% phi = right ascension and declination of the sun with respect to the Y
% axis
% eclipse = 1 for in eclipse, 0 for not in eclipse
% t = the time from equatorial pass (X-Y plane) of the phi angle readings
% alt = 600km
% inc = inclination of the simulated orbit from the Z axis (0 to 180 deg)
% launch_date = date of launch in day # of the year (day #1-365)
% n = number of points to calculate for 1 orbit
% NOTE: units are in metric, km, kg, s
% NOTE: Launch Complex 41 = 28d,35s N / 80d,34s,58as W
% NOTE: Earth's axis is oblique 23.4 deg
close all;
cl = clock;

% dates and time info
launch_loc_celef = [28.388333,-80.603611];
j2000_equinox = datenum(2000,03,20);
mission_start = datenum(launch_date);
year = 365;
launch_toy_eci = (mission_start-j2000_equinox)/year-...
    floor((mission_start-j2000_equinox)/year)-launch_loc_celef(2)/(360*365);
dl_toy_eci = launch_toy_eci-(launch_loc_celef(2)-dl_loc(2))/360;

launch_loc_ecef = [cosd(launch_loc_celef(1)),...
    0,...
    sind(launch_loc_celef(1))]/...
    norm([cosd(launch_loc_celef(1)),...
    0,...
    sind(launch_loc_celef(1))]);
launch_loc_ecef = launch_loc_ecef';
Q_tol = [cosd(launch_toy_eci*360), -sind(launch_toy_eci*360), 0;
         sind(launch_toy_eci*360), cosd(launch_toy_eci*360), 0;
         0, 0, 1];
launch_loc_eci = Q_tol*launch_loc_ecef;
launch_loc_eci = launch_loc_eci*(6378);

dl_loc_ecef = [cosd(dl_loc(1)),0,sind(dl_loc(1))]/...
    norm([cosd(dl_loc(1)),0,sind(dl_loc(1))]);
dl_loc_ecef = dl_loc_ecef';

%% calculation ground station FOV projection onto the celestial sphere
% dist_to_orbit: distance from gnd sta to SDS-1 in orbit
% see 9_COMM/Ground Station/Satellite Position Determination/ for reference
re = 6378.1;
gnd_station_fov = 160;
dist_to_orbit = (-2.*re.*cosd(gnd_station_fov/2)+(4.*re.^2.*...
    cosd(gnd_station_fov/2).^2-4.*(re^2-(re+alt)^2)).^(.5))./2;
z = sind(gnd_station_fov/2)*dist_to_orbit;
gnd_station_fov_proj = asind(z/(re+alt));

%% SDS-1 Antenna FOV Setup
% setup #1: antenna boresight at 90 deg to FIDI boresight, +/- 30 deg
% setup #2: antenna boresight opposite to FIDI boresight, +/- 30 deg

sc_ant_fov = [90,60];

%% spherical right triangle to find RA of Node; sin(b) = tan(DEC)*cot(inc)
b = asind(tand(launch_loc_celef(1))*cotd(inc));
ra_N = launch_toy_eci*360 - b;
omega = acosd(cosd(launch_loc_celef(1))*cosd(b));

X = [1,0,0];
Y = [0,1,0];
Z = [0,0,1];

%% determining the transformation matrix from perifocal to geocentric
R_ra_N = [cosd(ra_N),sind(ra_N),0; -sind(ra_N),cosd(ra_N),0; 0,0,1];
R_inc = [1,0,0; 0,cosd(inc),sind(inc); 0,-sind(inc),cosd(inc)];
R_omega = [cosd(omega),sind(omega),0; -sind(omega),cosd(omega),0; 0,0,1];
Q = R_omega*R_inc*R_ra_N;
Q_t = Q';

inc_sun = 23.4;
R_ra_N_sun = [1,0,0; 0,1,0; 0,0,1];
R_inc_sun = [1,0,0; 0,cosd(inc_sun),sind(inc_sun); 0,-sind(inc_sun),cosd(inc_sun)];
R_omega_sun = [cosd(launch_toy_eci*360), sind(launch_toy_eci*360), 0;
                -sind(launch_toy_eci*360), cosd(launch_toy_eci*360), 0;
                0, 0, 1];
Q_sun = R_omega_sun*R_inc_sun*R_ra_N_sun;
Q_t_sun = Q_sun';

%% determining the perifocal orbit elements (r and v)
mu = 398600;
r_e = 6378.1;
r = r_e+alt;
h = sqrt(r*mu);
h = h(1);
T = 2*pi/sqrt(mu)*r^(3/2);
theta_dot = 360/T;
T_mission = round(120*24*60*60/T)*T;
n_orbits = T_mission/T;
t_orbit = linspace(0,T,n);
theta = t_orbit*theta_dot;

t_sat = linspace(0,T*n_orbits,n*n_orbits);
theta_sat = t_sat*theta_dot;
theta_earth_dot = (7.2921150*10^(-5))/(pi)*180;

% change from CELESTIAL ECI to CELESTIAL EARTH FIXED (CELEF)
% CELECI + correction = CELEF

celeci_to_celef = -launch_toy_eci*360-theta_earth_dot*t_sat+launch_loc_celef(2);

% determining the sun's period, and then it's angular rotation about the
% ECI reference frame, aka the angular rotation of earth's orbit around the
% sun; unfortunately the earth is NOT the center of the universe...
t_sun = linspace(0,T_mission,n*n_orbits);
T_sun = 365.25*86400;
theta_dot_sun = 360/T_sun;
theta_sun = t_sun*theta_dot_sun;
r_sun = 10000;

%% determining the geocentric orbit elements (r and v)
r_pf = zeros(3,n);
v_pf = zeros(3,n);
r_gc = zeros(3,n);
v_gc = zeros(3,n);
r_pf_sun = zeros(3,n*n_orbits);
r_gc_sun = zeros(3,n*n_orbits);
r_sat = zeros(3,n*n_orbits);
r_sat_celef = zeros(2,n*n_orbits);
r_sat_celef_c = zeros(2,n*n_orbits);
r_sat_info = struct('long',zeros(1,n*n_orbits),'lat',zeros(1,n*n_orbits),...
    'dl',zeros(1,n*n_orbits),'scpv',zeros(3,n*n_orbits),'gnd2scv',zeros(3,n*n_orbits));
r_gnd_sta = zeros(3,n*n_orbits);

for i = 1:n
    r_pf(1:3,i) = h^2/mu*[cosd(theta(i));sind(theta(i));0];
    v_pf(1:3,i) = mu/h*[-sind(theta(i));cosd(theta(i));0];
    r_gc(1:3,i) = Q_t*r_pf(1:3,i);
    v_gc(1:3,i) = Q_t*v_pf(1:3,i);
end


for i = 1:n*n_orbits
    waitbar(i/(n*n_orbits));
    
    r_pf_sun(1:3,i) = r_sun*[cosd(theta_sun(i));sind(theta_sun(i));0];
    r_gc_sun(1:3,i) = Q_t_sun*r_pf_sun(1:3,i);
    
    r_sat(1:3,i) = Q_t*[cosd(theta_sat(i));sind(theta_sat(i));0];
    if r_sat(2,i) <= 0 && r_sat(1,i) < 0
        r_sat_celef(1:2,i) = [asind(r_sat(3,i)/norm(r_sat(1:3,i)));atand(r_sat(2,i)/r_sat(1,i))-180];
    elseif r_sat(2,i) > 0 && r_sat(1,i) < 0
        r_sat_celef(1:2,i) = [asind(r_sat(3,i)/norm(r_sat(1:3,i)));atand(r_sat(2,i)/r_sat(1,i))+180];
    else
        r_sat_celef(1:2,i) = [asind(r_sat(3,i)/norm(r_sat(1:3,i)));atand(r_sat(2,i)/r_sat(1,i))];
    end
    
    % adjusting the satellite celestial earth fixed (r_sat_celef)
    % coordinates to account for earth's rotation and GMT's location
    r_sat_celef_c(2,i) = r_sat_celef(2,i)+ 360 + celeci_to_celef(i);

    while r_sat_celef_c(2,i) < -180
            r_sat_celef_c(2,i) = r_sat_celef_c(2,i) + 360;
    end
    while r_sat_celef_c(2,i) > 180
            r_sat_celef_c(2,i) = r_sat_celef_c(2,i) - 360;
    end
    
    % adding downlink capability information to the r_sat_celef_c matrix
    % r_sat_info = {[lat, long, D/L Avail.?; 
    %               SC Pointing Vector; 
    %               GND STA to SC Pointing Vector;]}
    
    Q_earth = [cosd(dl_toy_eci*360+theta_earth_dot*t_sat(i)-360), -sind(dl_toy_eci*360+theta_earth_dot*t_sat(i)-360), 0;
             sind(dl_toy_eci*360+theta_earth_dot*t_sat(i)-360), cosd(dl_toy_eci*360+theta_earth_dot*t_sat(i)-360), 0;
             0, 0, 1];
    r_gnd_sta(1:3,i) = Q_earth*dl_loc_ecef;
    
    % testing to see if SC is within the ground stations FOV
    % gnd_station_fov_proj = furthest true anomaly of SDS-1 away from GND
    %   STAT, in the GND STAT's F.O.V.
    % r_sat: vector from E.C. to SDS-1
    % r_gnd_sta: vector from E.C. to GND STAT
    
    r_sat_info.lat(i) =       r_sat_celef(1,i); 
    r_sat_info.long(i) =        r_sat_celef_c(2,i);
    r_sat_info.scpv(1:3,i) =       r_gc_sun(1:3,i);
    r_sat_info.gnd2scv(1:3,i) =    r_sat(1:3,i)-r_gnd_sta(1:3,i);
    
    if acosd(dot(r_sat(1:3,i),r_gnd_sta(1:3,i))/...
            (norm(r_sat(1:3,i))*norm(r_gnd_sta(1:3,i)))) <= gnd_station_fov_proj
        if acosd(dot(r_sat_info.gnd2scv(1:3,i),r_sat_info.scpv(1:3,i))/...
            (norm(r_sat_info.gnd2scv(1:3,i))*norm(r_sat_info.scpv(1:3,i)))) <= sc_ant_fov(1)+sc_ant_fov(2)...
             && acosd(dot(r_sat_info.gnd2scv(1:3,i),r_sat_info.scpv(1:3,i))/...
             (norm(r_sat_info.gnd2scv(1:3,i))*norm(r_sat_info.scpv(1:3,i)))) >= sc_ant_fov(1)-sc_ant_fov(2)
            r_sat_info.dl(i) =         1;
        else
            r_sat_info.dl(i) =         0;
        end
    else
        r_sat_info.dl(i) =         0;
    end
end

% determining eclipse angle between sun/gnd-sc vectors
eclipse_angle = 180-asind(re/(re+alt));

%% post processing on downlinking capabilities -----------------------------
dl_count = 0;
timestep = t_sat(2);
downlink_times = 0;
downlink_info = 0;
downlink_angle = cell(1);

for i = 1:n*n_orbits
    if r_sat_info.dl(i) == 1 && r_sat_info.dl(i+1) == 1
        dl_count = dl_count+1;
        downlink_info(end+1,1) = acosd(dot(r_sat_info.scpv(1:3,i),r_sat_info.gnd2scv(1:3,i))/...
            (norm(r_sat_info.scpv(1:3,i))*norm(r_sat_info.gnd2scv(1:3,i))));
        
        % testing to see if SC is in eclipse
        if acosd(dot(r_sat_info.scpv(1:3,i),r_sat(1:3,i))/...
            (norm(r_sat_info.scpv(1:3,i))*norm(r_sat(1:3,i)))) >= eclipse_angle
            downlink_info(end,2) = 1;
        end
        
    elseif r_sat_info.dl(i) == 1 && r_sat_info.dl(i+1) == 0
        dl_count = dl_count+1;
        downlink_times(end+1) = dl_count*timestep;
        downlink_info(end+1,1) = acosd(dot(r_sat_info.scpv(1:3,i),r_sat_info.gnd2scv(1:3,i))/...
            (norm(r_sat_info.scpv(1:3,i))*norm(r_sat_info.gnd2scv(1:3,i))));
        
        % testing to see if SC is in eclipse
        if acosd(dot(r_sat_info.scpv(1:3,i),r_sat(1:3,i))/...
            (norm(r_sat_info.scpv(1:3,i))*norm(r_sat(1:3,i)))) >= eclipse_angle
            downlink_info(end,2) = 1;
        end
        
        downlink_info = downlink_info(2:end,:);
        downlink_angle(end+1) = {downlink_info};
        downlink_info = 0;
        dl_count = 0;
    elseif i == n*n_orbits && r_sat_info.dl(i) == 1
        dl_count = dl_count+1;
        downlink_times(end+1) = dl_count*timestep;
        downlink_info(end+1,1) = acosd(dot(r_sat_info.scpv(1:3,i),r_sat_info.gnd2scv(1:3,i))/...
            (norm(r_sat_info.scpv(1:3,i))*norm(r_sat_info.gnd2scv(1:3,i))));
        
        % testing to see if SC is in eclipse
        if acosd(dot(r_sat_info.scpv(1:3,i),r_sat(1:3,i))/...
            (norm(r_sat_info.scpv(1:3,i))*norm(r_sat(1:3,i)))) >= eclipse_angle
            downlink_info(end,2) = 1;
        end
        
        downlink_info = downlink_info(2:end,:);
        downlink_angle(end+1) = {downlink_info};
        downlink_info = 0;
        dl_count = 0;
    end
end
downlink_times = downlink_times(2:end);
downlink_angle = downlink_angle(2:end);
all_downlink_angles = 0;

for i = 1:length(downlink_angle);
    angles = cell2mat(downlink_angle(i));
    all_downlink_angles(end+1:end+length(angles(:,1))) = angles(:,1);
    downlink_angle_ranges(i) = mean(angles(:,1));
end
all_downlink_angles = all_downlink_angles(2:end);

mean_downlink_time = mean(downlink_times);
mean_downlink_count = length(downlink_times)/120;


%% Saving...
savefilename = ['ant_fov_sim_90_60_',num2str(cl(2)),'-',num2str(cl(3)),'-',num2str(cl(4)),'-',num2str(cl(5))];
save(savefilename);

%% plotting
figure;
hist(downlink_times,20)
title('Downlink Time Histogram')
xlabel('Downlink Time (s)')
ylabel('Count')

figure;
hist(all_downlink_angles,50)
title('Downlink Angle Histogram')
xlabel('Downlink Angle (deg)')
ylabel('Count')

figure;
hist(downlink_angle_ranges,50)
title('Downlink Angle Range Histogram')
xlabel('Downlink Angle Range (deg)')
ylabel('Count')

figure
plot3(r_gc(1,:),r_gc(2,:),r_gc(3,:),'r','linewidth',2);
hold
plot3(r_gc_sun(1,:),r_gc_sun(2,:),r_gc_sun(3,:),'y','linewidth',2);
plot3([0,X(1)]*8000,[0,X(2)]*8000,[0,X(3)]*8000,'r','linewidth',2)
plot3([0,Y(1)]*8000,[0,Y(2)]*8000,[0,Y(3)]*8000,'g','linewidth',2)
plot3([0,Z(1)]*8000,[0,Z(2)]*8000,[0,Z(3)]*8000,'b','linewidth',2)
plot3([0,launch_loc_eci(1)],[0,launch_loc_eci(2)],[0,launch_loc_eci(3)],'r','linewidth',2)
earth_sphere(launch_toy_eci*360-launch_loc_celef(2));
axis equal square

r_sat_celef_plot(1,:) = r_sat_celef(1,:)+90;
r_sat_celef_plot(2,:) = r_sat_celef_c(2,:)+180;

figure
ratio = 1/4;
convert = (1250/180)*ratio;
B = imread('earth_map.jpg');
C = imresize(B, ratio);
C = flipdim(C,1);
imshow(C);
hold
plot(r_sat_celef_plot(2,1:360)*convert,r_sat_celef_plot(1,1:360)*convert,'g.');
plot(r_sat_celef_plot(2,361:720)*convert,r_sat_celef_plot(1,361:720)*convert,'b.');
plot(r_sat_celef_plot(2,721:5760)*convert,r_sat_celef_plot(1,721:5760)*convert,'r.');
for i = 1:5760
    if r_sat_info.dl(i) == 1
        plot(r_sat_celef_plot(2,i)*convert,r_sat_celef_plot(1,i)*convert,'y.');
    end
end
title('SDS-1 Groundtrack Orbits #1-16');
xlabel('Right Ascension (deg)');
ylabel('Declination (deg)');
set(gca, 'ydir', 'normal');
