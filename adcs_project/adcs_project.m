function varargout = adcs_project(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @adcs_project_OpeningFcn, ...
                   'gui_OutputFcn',  @adcs_project_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% --- Executes just before adcs_project is made visible.
function adcs_project_OpeningFcn(hObject, eventdata, handles, varargin)
set(handles.pitch0,'string',25)
set(handles.roll0,'string',90)
set(handles.yaw0,'string',0)

set(handles.pitch,'string',0 )
set(handles.roll,'string',0)
set(handles.yaw,'string',0)
% w0=[0.53;0.53;0.053]*pi/180;

set(handles.w1,'string',0.53)
set(handles.w2,'string',0.53)
set(handles.w3,'string',0.053)

set(handles.T,'string',300)
set(handles.kp,'string',50)
set(handles.kd,'string',500)

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% --- Outputs from this function are returned to the command line.
function varargout = adcs_project_OutputFcn(hObject, eventdata, handles) 
% pitch0

% 
% axes(1)
% axis([-10 10 0 10])


varargout{1} = handles.output;


% --- Executes on button press in START.
function START_Callback(hObject, eventdata, handles)
rad2deg=180/pi;

roll0_deg=str2double(get(handles.roll0,'string'));
yaw0_deg=str2double(get(handles.yaw0,'string'));
pitch0_deg=str2double(get(handles.pitch0,'string'));

roll_deg=str2double(get(handles.roll,'string'));
pitch_deg=str2double(get(handles.pitch,'string'));
yaw_deg=str2double(get(handles.yaw,'string'));

w1_deg=str2double(get(handles.w1,'string'));
w2_deg=str2double(get(handles.w2,'string'));
w2_deg=str2double(get(handles.w3,'string'));


tf=str2double(get(handles.T,'string'));
kp=str2double(get(handles.kp,'string'));
kd=str2double(get(handles.kd,'string'));


% % code
dt=1;
t=[0:dt:tf]';
m=length(t);
j=diag([10000 9000 12000]);
invj=inv(j);

% Initial and Desired Quaternion, and Initial Rate
eul0_rad=[roll0_deg  pitch0_deg yaw0_deg ]*pi/180;
eul_rad=[roll_deg  pitch_deg yaw_deg ]*pi/180;


quatd = eul2quat(eul_rad,'XYZ')
quat0 = eul2quat(eul0_rad,'XYZ')

qd=quatd';
q0=quat0';
q0=q0/norm(q0);
w0=[0.53;0.53;0.053]*pi/180;
x=zeros(m,7);
x(1,:)=[q0' w0'];

% Gains
kp=str2double(get(handles.kp,'string'));
kd=str2double(get(handles.kd,'string'));

%% Main Loop
for i = 1:m-1
 f1=dt*euler_quat_reg_fun(x(i,:),j,invj,qd,kp,kd);
 f2=dt*euler_quat_reg_fun(x(i,:)+0.5*f1',j,invj,qd,kp,kd);
 f3=dt*euler_quat_reg_fun(x(i,:)+0.5*f2',j,invj,qd,kp,kd);
 f4=dt*euler_quat_reg_fun(x(i,:)+f3',j,invj,qd,kp,kd);
 x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
end

% Error Quaternion
qerr=quat_err(x(:,1:4),[qd(1)*ones(m,1) qd(2)*ones(m,1) qd(3)*ones(m,1) qd(4)*ones(m,1)]);

% Torque
%torque=-kp*qerr(:,1:3)-kd*x(:,5:7);
torque=-kp*kron(sign(qerr(:,4)),[1 1 1]).*qerr(:,1:3)-kd*x(:,5:7);


%% Plot Results
% clf
axes(handles.axes2);
% plot(t,qerr(:,1))
set(gca,'FontSize',4)
axis([0 300 -0.2 0.8])
set(gca,'Ytick',[-0.2 0 0.2 0.4 0.6 0.8])
set(gca,'fontsize',12)
ylabel('dq1')

axes(handles.axes3);
% plot(t,qerr(:,2))
set(gca,'FontSize',.5)
axis([0 300 -0.2 0.8])
set(gca,'Ytick',[-0.2 0 0.2 0.4 0.6 0.8])
set(gca,'fontsize',12)
ylabel('dq2')

axes(handles.axes4);
% plot(t,qerr(:,3))
set(gca,'FontSize',6)
axis([0 300 -0.2 0.8])
set(gca,'Ytick',[-0.2 0 0.2 0.4 0.6 0.8])
set(gca,'fontsize',12)
ylabel('dq3')

axes(handles.axes5);
% plot(t,qerr(:,4))
axis([0 300 0 1.2])
set(gca,'Ytick',[-0.2 0 0.2 0.4 0.6 0.8])
set(gca,'fontsize',12)
ylabel('dq4')


% % % % % % % % % % % % % % 
axes(handles.axes6);
set(gca,'FontSize',3)
axis([0 300 -1 1])
% plot(t,x(:,5))
set(gca,'fontsize',12)
ylabel('w1 (rad/sec)')

axes(handles.axes7);
set(gca,'FontSize',3)
% plot(t,x(:,6))
set(gca,'fontsize',12)
ylabel('w2 (rad/sec)')

axes(handles.axes8);
set(gca,'FontSize',3)
% plot(t,x(:,6))
set(gca,'fontsize',12)
ylabel('w3 (rad/sec)')


% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
axes(handles.axes9);
% plot(t,torque(:,1))
set(gca,'fontsize',12)
ylabel('L1 (N-m)')

axes(handles.axes12);
% plot(t,torque(:,2))
set(gca,'fontsize',12)
ylabel('L2 (N-m)')

axes(handles.axes13);
% plot(t,torque(:,3))
set(gca,'fontsize',12)
ylabel('L3 (N-m)')

%% live plot
axes(handles.axes1)
targetHandle = DrawAxes(0.1);
targetHandle.Matrix = makehgtform('zrotate', eul_rad(3),'yrotate', eul_rad(2),'xrotate',eul_rad(1));
currentHandle = DrawAxes(1);
%axis auto
%axis equal
%axis vis3d 
axis([-1 1 -1 1 -1 2])

    for ii = 1 : length(t)
    currentHandle.Matrix = quat2tform(x(ii, 1:4));
    drawnow
    
    axes(handles.axes2);
        set(gca,'FontSize',7)
    hold on
    plot([1:ii*dt],qerr([1:ii],1))

     axes(handles.axes3);
    hold on
    set(gca,'FontSize',7)
    plot([1:ii*dt],qerr([1:ii],2))
    
    axes(handles.axes4);
        set(gca,'FontSize',7)
    hold on
    plot([1:ii*dt],qerr([1:ii],3))
    
    axes(handles.axes5);
        set(gca,'FontSize',7)
    hold on
    plot([1:ii*dt],qerr([1:ii],4))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    axes(handles.axes6);
    hold on
    set(gca,'FontSize',5)
    axis([0 tf   -1.1*abs(min(x(:,5)))   1.1*max(x(:,5))]) ;
    plot([1:ii*dt],x([1:ii],5))
    
    axes(handles.axes7);
    hold on
    set(gca,'FontSize',5)
    axis([0 tf   -1.1*abs(min(x(:,6)))   1.1*max(x(:,6))]) ;
    plot([1:ii*dt],x([1:ii],6))
    
    axes(handles.axes8);
    hold on
    set(gca,'FontSize',7)
    axis([0 tf   -1.1*abs(min(x(:,7)))   1.1*max(x(:,7))]) ;
    plot([1:ii*dt],x([1:ii],7))
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    axes(handles.axes9);
    axis([0 tf   -1.1*abs(min(torque(:,1)))   1.1*max(torque(:,1))]) ;
    hold on
    plot([1:ii*dt],torque([1:ii],1))
    
    axes(handles.axes12);
    axis([0 tf   -1.1*abs(min(torque(:,2)))   1.1*max(torque(:,2))]) 
    hold on
    plot([1:ii*dt],torque([1:ii],2))
   
    axes(handles.axes13);
    axis([0 tf   -1.1*abs(min(torque(:,3)))   1.1*max(torque(:,3))]) 
    hold on
    plot([1:ii*dt],torque([1:ii],3))
    

    
    end
