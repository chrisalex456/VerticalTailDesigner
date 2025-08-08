%% Aircraft Stability Analysis
% This script analyzes the effects
% of vertical tail parameters on aircraft stability.
% Assumptions: symmetric aircraft
% Constants and sweep parameters are defined for each section


%% Givens
clc; clear all;
W = 4000; % aircraft weight in lbs
b = 35; % wingspan in ft
mean_aChord = 6.5; % mean aerodynamic chord in ft
S = 227; % wing area in ft^2
Ix = 2500; % moment of inertia about x-axis in slug*ft^2
Iy = 20500; % moment of inertia about y-axis in slug*ft^2
Iz = 25000; % moment of inertia about z-axis in slug*ft^2
Ixy = 0; % product of inertia of x- and y- axis in slug*ft^2
h = 10000; % altitude in ft
rho = 0.001755; % air density in slug/ft^3
aSpeed_knots = 638.35; % speed of sound in knots
pressure = 1455.33; % psf
temp = 483; % Rankine
speed_knots = 200; % aircraft speed in knots
v = speed_knots*1.69; % aircraft velocity in ft/s
Cla_vertTail = 5.9; % 2D coefficient of lift
q = (rho*v^2)/2; % dynamic pressure in lb/ft^2

%% Lift Curve Slope/Lift

% sweep parameters
AR_vertTail = linspace(1.3, 2, 15); % vertical tail aspect ratio
S_vertTail = linspace(10, 25, 15); % vertical tail area (ft^2)

% initialized variables
CLa_corrected = zeros(length(AR_vertTail), length(S_vertTail));
L = zeros(length(AR_vertTail), length(S_vertTail));

for i = 1:length(AR_vertTail)
    for j = 1:length(S_vertTail)
        CLa_num = Cla_vertTail/(1+(Cla_vertTail/(pi*AR_vertTail(i))));
        L_num = CLa_num*0.5*rho*v^2*S_vertTail(i);
        CLa_corrected(i, j) = CLa_num;
        L(i, j) = L_num;
    end
end

% Graph of Lift as a Function of Vertical Tail Parameters
figure(1)
mesh(AR_vertTail, S_vertTail, L)
xlabel('Area Ratio')
ylabel('Vertical Tail Area (ft^2)')
zlabel('Vertical Tail Lift')
title('Effect of Vertical Tail Parameters on Lift')

%% Roll Static Stability
% Clb < 0 for roll static stability
% z_v should be positive

% constants
Clbw = 0; % wing position and wing sweep contribution
lambda = 0.4; % wing taper ratio
gamma = 3; % wing geometric dihedral (degrees)
CLa_wing = 1.5; % lift-curve slope for the wing
S_vertTail_f = 15.4; % vertical tail area (ft^2) (determined from lift section)
CLa_vertTail_c = 2.67; % corrected vertical tail coefficient of lift

% sweep parameters
h_vertTail = linspace(5, 10, 15); % vertical tail height (ft)
lambda_vertTail = linspace(0.2, 0.5, 15); % vertical tail taper ratio


% initialized variables
z_v = zeros(length(h_vertTail), length(S_vertTail));
Clb = zeros(length(h_vertTail), length(S_vertTail));

for j = 1:length(h_vertTail)
    for k = 1:length(lambda_vertTail)
        CLb_dihedral = -(1/6)*CLa_wing*gamma*((1+2*lambda)/(1+lambda));
        z_v_num = (1/3)*(h_vertTail(j))*((1+2*lambda_vertTail(k))/(1+lambda_vertTail(k))); % vertical distance of ac of vt above cg
        Cy_vertTail = -CLa_vertTail_c*(S_vertTail_f/S); % assumed the dynamic pressure ratio is one
        Clb_vertTail = (z_v_num/b)*Cy_vertTail;
        Clb_num = Clbw + Clb_vertTail + CLb_dihedral;
        z_v(j, k) = z_v_num;
        Clb(j, k) = Clb_num;
    end
end

% Graph of Roll Static Stability
figure(2)
mesh(h_vertTail, lambda_vertTail, Clb)
xlabel("Vertical Tail Height (ft)")
ylabel("Taper Ratio")
zlabel("Roll Stability")
title("Effect of Vertical Tail Parameters on Roll Stability")
grid on;

%% Directional Stability
% CnB > 0 for directional stability

% constants
CnB_wf = -0.10; % wing and fuselage contribution to CnB

% sweep parameters
x_v = linspace(10, 25, 15); % the x distance from the ac of the vt to the aircraft cg

%initialized variables
CnB = zeros(length(x_v), 1);

for i = 1:length(x_v)
    CnB_vertTail_num = CLa_vertTail_c*(S_vertTail_f*x_v(i)/(S*b));
    CnB(i, 1) = CnB_vertTail_num + CnB_wf;
end

%% Yawing Moment
% typically negative

% constants
Cn0 = 0; % for symmetric aircraft
B = 15; % max sideslip angle
Cnda = 0; % aileron contribution to yawing moment
da = 0; % max aileron deflection
dr = -35; % rudder deflection angle

% sweep parameters
% CnB from previous section

% initialized variables
Cn = zeros(length(CnB), length(x_v));
vc = 0.02326; % volume coefficient
deflection = 0.75; % rudder deflection efficiency

tau = 0.25;
Cndr =-vc*deflection*tau*CLa_vertTail_c;

for i = 1:length(CnB)
    for j = 1:length(x_v)
        Cn(i, j) = Cn0 + CnB(i)*B + Cnda*da + Cndr*dr;
    end
end

% Graph of Yawing Moment
figure(3)
mesh(x_v, CnB, Cn)
xlabel('x_v')
ylabel('CnB')
zlabel('Cn')
title('Cn as a Function of x_v and CnB')

%% Yaw Damping Derivative
% should be negative

% constants
x_v_num = 12; % feet
h_vertTail_num = 5; % feet
nv = 1; % number of vertical tails

Cnr_wing = -0.25*CLa_wing*(b/v);
inter = (2*x_v_num^2)/b^2;
S_ratio = S_vertTail_f/S;
Cnr_VT = -CLa_vertTail_c*inter*S_ratio;
Cnr = Cnr_wing + Cnr_VT;

%% Moments

% constants
CnB_num = 0.0182933;

NB = ((q*S*b)/Iz)*CnB_num; % yawing moment due to sideslip
Nr = ((q*S*b^2)/(2*Iz*v)*Cnr); % yawing moment due to yaw rate
Ndr = ((q*S*b)/Iz)*Cndr; % yawing moment due to rudder deflection

w = sqrt(NB); % natural frequency
zeta = -Nr/(2*w); % damping ratio
