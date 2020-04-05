function varargout = NSFGUI(varargin)
% NSFGUI MATLAB code for NSFGUI.fig
%      NSFGUI, by itself, creates a new NSFGUI or raises the existing
%      singleton*.
%
%      H = NSFGUI returns the handle to a new NSFGUI or the handle to
%      the existing singleton*.
%
%      NSFGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NSFGUI.M with the given input arguments.
%
%      NSFGUI('Property','Value',...) creates a new NSFGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NSFGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NSFGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NSFGUI

% Last Modified by GUIDE v2.5 13-Aug-2019 16:03:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NSFGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NSFGUI_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before NSFGUI is made visible.
function NSFGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NSFGUI (see VARARGIN)

% Choose default command line output for NSFGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.termination_state,'Visible','off');
set(handles.table1,'Visible','off')
set(handles.table2,'Visible','off')
set(handles.save_dir,'Visible','off');
set(handles.browse_push,'Visible','off');
set(handles.edit11,'Visible','off');
set(handles.text9,'Visible','off');
% UIWAIT makes NSFGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NSFGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function N_val_Callback(hObject, eventdata, handles)
% hObject    handle to N_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N_val as text
%        str2double(get(hObject,'String')) returns contents of N_val as a double


% --- Executes during object creation, after setting all properties.
function N_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function M_val_Callback(hObject, eventdata, handles)
% hObject    handle to M_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of M_val as text
%        str2double(get(hObject,'String')) returns contents of M_val as a double


% --- Executes during object creation, after setting all properties.
function M_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function K_val_Callback(hObject, eventdata, handles)
% hObject    handle to K_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of K_val as text
%        str2double(get(hObject,'String')) returns contents of K_val as a double


% --- Executes during object creation, after setting all properties.
function K_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w_val_Callback(hObject, eventdata, handles)
% hObject    handle to w_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w_val as text
%        str2double(get(hObject,'String')) returns contents of w_val as a double


% --- Executes during object creation, after setting all properties.
function w_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in alpha_check.
function alpha_check_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of alpha_check


% --- Executes on button press in beta_check.
function beta_check_Callback(hObject, eventdata, handles)
% hObject    handle to beta_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of beta_check



function alpha_val_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_val as text
%        str2double(get(hObject,'String')) returns contents of alpha_val as a double


% --- Executes during object creation, after setting all properties.
function alpha_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beta_val_Callback(hObject, eventdata, handles)
% hObject    handle to beta_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta_val as text
%        str2double(get(hObject,'String')) returns contents of beta_val as a double


% --- Executes during object creation, after setting all properties.
function beta_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in start_push.
function start_push_Callback(hObject, eventdata, handles)
% hObject    handle to start_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(handles.alpha_val,'enable','off');
% set(handles.beta_val,'enable','off');
global stopit
stopit=0;
set(handles.N_val,'enable','off');
set(handles.M_val,'enable','off');
set(handles.K_val,'enable','off');
set(handles.w_val,'enable','off');
set(handles.alpha_val,'enable','off');
set(handles.beta_val,'enable','off');
set(handles.max_run_time,'enable','off');
set(handles.save_check,'enable','off');
set(handles.reset_push,'Enable','off');
set(handles.stop_push,'Enable','on');

% drawnow;
% set(handles.termination_state,'Visible','on');
% drawnow;
% set(handles.termination_state,'String','Please enter the input');
% drawnow;
% set(handles.table1,'Visible','on');
% drawnow;
% set(handles.table2,'Visible','on');



if isempty((handles.N_val.String)) ||  isempty((handles.M_val.String)) || isempty((handles.K_val.String)) || isempty((handles.w_val.String))
    
    set(handles.termination_state,'Visible','on');
    drawnow;
    set(handles.termination_state,'String','Please enter all the inputs');
    drawnow;
    set(handles.stop_push,'Enable','on');
    set(handles.reset_push,'Enable','on');
else
    set(handles.termination_state,'Visible','on');
    drawnow;
    set(handles.termination_state,'String','intlinprog is running');
    drawnow;
end


N = str2num(handles.N_val.String);
M = str2num(handles.M_val.String);
K = str2num(handles.K_val.String);
w = str2num(handles.w_val.String);
alpha = str2num(handles.alpha_val.String);
e = str2num(handles.beta_val.String);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = w';

ycoeff = c(:);
f = ycoeff;
for k=1:K
    f = [f;zeros(N*M,1)];
end
f = [f;zeros(N*M*M,1)];

X = [];
Q = [];
for k = 1:K+1
    P = zeros(M,N*M);                  %2 3 4 5 6 RHS
    Y = [];
    for j=1:M
        y = zeros(N,M);
        for i=1:N
            if w(i,j)~=0
                y(i,j) = 1;
            end
        end
        y = y';
        y = y(:);
        Y = [Y;y'];
    end
    P(:,:,k) = Y;
    Q((k-1)*M+1:k*M,:,k) = P(:,:,k);
end
X = reshape(Q,(K+1)*M,(K+1)*N*M);
b = 1*ceil(N*K/M)*ones(M,1);
b = [b;1*ceil(N/M)*ones(K*M,1)];


Q = [];
for k = 1:K+1
    P = zeros(M,N*M);                  %2 3 4 5 6 LHS
    Y = [];
    for j=1:M
        y = zeros(N,M);
        for i=1:N
            if w(i,j)~=0
                y(i,j) = -1;
            end
        end
        y = y';
        y = y(:);
        Y = [Y;y'];
    end
    
    P(:,:,k) = Y;
    Q((k-1)*M+1:k*M,:,k) = P(:,:,k);
end
dd = reshape(Q,(K+1)*M,(K+1)*N*M);
X = [X;dd];
b = [b;-1*floor(N*K/M)*ones(M,1)];
b = [b;-1*floor(N/M)*ones(K*M,1)];

Y=[];                                      %13,14,15,16
Q = [];
P = [];
for k = 1:K
    S = [];
    for i=1:N
        for j=1:M
            if w(i,j)~=0
                s = zeros(N,M);
                s(i,j) = 1;
                s = s';
                s = s(:);
                S = [S;s'];
                
                y = zeros(N,M);
                y(i,j) = -1;
                y = y';
                y = y(:);
                Y = [Y;y'];
            else
                s = zeros(1,N*M);
                S = [S;s];
                
                y = zeros(1,N*M);
                Y = [Y;y];
            end
        end
    end
    P(:,:,k) = S;
    Q((k-1)*M*N+1:k*N*M,(k-1)*M*N+1:k*N*M,k) = P(:,:,k);
end

PP = 0;
for k=1:K
    PP = PP+Q(:,:,k);
end
Q = [Y,PP];
X = [X;Q];


b = [b;zeros((K)*N*M,1)];


X = [X,zeros((N*M*K)+(M*(K+1)*2) , N*M*M)];
Q = [];                                    %17,18,19
L = [];
S = [];
a=1;
for i = 1:N
    for j = 1:M
        for jdash = 1:M
            if j~= jdash && w(i,j)~=0 && w(i,jdash)~=0
                l = zeros(N,M);
                s = zeros(N,M);
                l(i,j)= 1;
                s(i,jdash)=1;
                l = l';
                l = l(:);
                s = s';
                s = s(:);
                L = [L;l'];
                S = [S;s'];
                xx(a,a) = 1;
                pp(a) = 1;
                a=a+1;
            else
                l = zeros(1,N*M);
                s = zeros(1,N*M);
                L = [L;l];
                S = [S;s];
                xx(a,a) = 0;
                pp(a)=0;
                a=a+1;
            end
        end
    end
end
ZZ = [L,S];
for k=1:K-1
    Q((k-1)*M*M*N+1:(k)*M*M*N,(k-1)*M*N+1:(k+1)*M*N) = ZZ;
end

y = zeros(N*M*M*(K-1),N*M);


Q = [y,Q];
xxx = [];
ppp = [];
for k = 1:K-1
    xxx = [xxx;xx];
    ppp = [ppp;pp'];
end
Q = [Q,-1*xxx];
X = [X;Q];

b = [b;ppp];


Xeq = [];
Qeq = [];
for k = 1:K+1
    P = zeros(N,N*M);                  %7 8 9 10 11
    Y = [];
    for i=1:N
        y = zeros(N,M);
        for j=1:M
            if w(i,j)~=0
                y(i,j) = 1;
            end
        end
        y = y';
        y = y(:);
        Y = [Y;y'];
    end
    P(:,:,k) = Y;
    Qeq((k-1)*N+1:k*N,:,k) = P(:,:,k);
end
Xeq = reshape(Qeq,(K+1)*N,(K+1)*N*M);
beq = K*ones(N,1);
beq = [beq;1*ones(K*N,1)];


QQ = [];                             %12
for i=1:N
    for j=1:M
        if w(i,j)~=0
            qq = zeros(N,M);
            qq(i,j) = 1;
            qq = qq';
            qq = qq(:);
            QQ = [QQ;qq'];
        else
            qq = zeros(1,N*M);
            QQ = [QQ;qq];
        end
    end
end
RR = [];
for k = 1:K
    RR = [RR,QQ];
end
RR = [-1*QQ,RR];
Xeq = [Xeq;RR];
beq = [beq;zeros(N*M,1)];


Xeq = [Xeq,zeros(N*(K+1)+N*M,N*M*M)];


e = 0.001;
alpha = 1000;
[rows1,~] = size(X);
X = [X,zeros(rows1,N*M*M)];

[rows2,~] = size(Xeq);
Xeq = [Xeq,zeros(rows2,N*M*M)];


a = 1;
ww1 = zeros(1,N*M*M);
ww2 = zeros(1,N*M*M);
sl = zeros(N*M*M);
for i = 1:N
    for j = 1:M
        for jdash = 1:M
            if j~=jdash && w(i,j)~=0 && w(i,jdash)~=0
                ww1(a) = w(i,j);
                ww2(a) = w(i,jdash);
                sl(a,a) = 0;
                epsilon(a) = e;
                a = a+1;                                                                                                            kks=0;
            else
                ww1(a) = 0;
                ww2(a) = 0;
                sl(a,a) = 0;
                epsilon(a) = 0;
                a = a+1;
            end
        end
    end
end
xx2 = N*xx;
lhs1 = [zeros(N*M*M,N*M*(K+1)),xx2,-1*sl];
rhs1 = -1*ww1'+ww2'+N;
X = [X;lhs1];
b = [b;rhs1];

lhs2 = [zeros(N*M*M,N*M*(K+1)),-1*xx2,sl];
rhs2 = ww1'-ww2'-epsilon';

X = [X;lhs2];
b = [b;rhs2];



intcon = (1:((K+1)*N*M+N*M*M));
lb = zeros(1,((K+1)*N*M)+2*N*M*M);
ub = ones(1,((K+1)*N*M+N*M*M));
f = [f;alpha*pp'];
if ~stopit
    options = optimoptions('intlinprog','IntegerTolerance',1e-6,'MaxTime',str2num(handles.max_run_time.String),'OutputFcn',@stopmyiter);
    [c,fval,exitflag] = intlinprog(f,intcon,X,b,Xeq,beq,lb,ub,options);
else
    set(handles.termination_state,'String','Intlinprog has stopped');
    drawnow;
end
switch exitflag
    case 2
        exitval = 'intlinprog stopped prematurely. Integer feasible point found.';
    case 1
        exitval = 'intlinprog converged to the solution x.';
    case 0
        exitval = 'intlinprog stopped prematurely. No integer feasible point found.';
    case -1
        exitval = 'intlinprog stopped by user.';
    case -2
        exitval = 'No feasible point found in model 1. Running model 2 ';
    case -3
        exitval = 'Root LP problem is unbounded.';
end

set(handles.termination_state,'String',exitval);
drawnow;

if exitflag==-2
    X = [];
    Q = [];
    for k = 1:K+1
        P = zeros(M,N*M);                  %2 3 4 5 6 RHS
        Y = [];
        for j=1:M
            y = zeros(N,M);
            for i=1:N
                if w(i,j)~=0
                    y(i,j) = 1;
                end
            end
            y = y';
            y = y(:);
            Y = [Y;y'];
        end
        P(:,:,k) = Y;
        Q((k-1)*M+1:k*M,:,k) = P(:,:,k);
    end
    X = reshape(Q,(K+1)*M,(K+1)*N*M);
    b = 1*ceil(N*K/M)*ones(M,1);
    b = [b;1*ceil(N/M)*ones(K*M,1)];
    
    
    Q = [];
    for k = 1:K+1
        P = zeros(M,N*M);                  %2 3 4 5 6 LHS
        Y = [];
        for j=1:M
            y = zeros(N,M);
            for i=1:N
                if w(i,j)~=0
                    y(i,j) = -1;
                end
            end
            y = y';
            y = y(:);
            Y = [Y;y'];
        end
        
        P(:,:,k) = Y;
        Q((k-1)*M+1:k*M,:,k) = P(:,:,k);
    end
    dd = reshape(Q,(K+1)*M,(K+1)*N*M);
    X = [X;dd];
    b = [b;-1*floor(N*K/M)*ones(M,1)];
    b = [b;-1*floor(N/M)*ones(K*M,1)];
    %X*y-b;
    
    
    Y=[];                                      %13,14,15,16
    Q = [];
    P = [];
    for k = 1:K
        S = [];
        for i=1:N
            for j=1:M
                if w(i,j)~=0
                    s = zeros(N,M);
                    s(i,j) = 1;
                    s = s';
                    s = s(:);
                    S = [S;s'];
                    
                    y = zeros(N,M);
                    y(i,j) = -1;
                    y = y';
                    y = y(:);
                    Y = [Y;y'];
                else
                    s = zeros(1,N*M);
                    S = [S;s];
                    
                    y = zeros(1,N*M);
                    Y = [Y;y];
                end
            end
        end
        P(:,:,k) = S;
        Q((k-1)*M*N+1:k*N*M,(k-1)*M*N+1:k*N*M,k) = P(:,:,k);
    end
    
    PP = 0;
    for k=1:K
        PP = PP+Q(:,:,k);
    end
    Q = [Y,PP];
    X = [X;Q];
    
    
    b = [b;zeros((K)*N*M,1)];
    
    
    X = [X,zeros((N*M*K)+(M*(K+1)*2) , N*M*M)];
    Q = [];                                    %17,18,19
    L = [];
    S = [];
    a=1;
    xx = [];
    for i = 1:N
        for j = 1:M
            for jdash = 1:M
                if j~= jdash && w(i,j)~=0 && w(i,jdash)~=0
                    l = zeros(N,M);
                    s = zeros(N,M);
                    l(i,j)= 1;
                    s(i,jdash)=1;
                    l = l';
                    l = l(:);
                    s = s';
                    s = s(:);
                    L = [L;l'];
                    S = [S;s'];
                    xx(a,a) = 1;
                    pp(a) = 1;
                    a=a+1;
                else
                    l = zeros(1,N*M);
                    s = zeros(1,N*M);
                    L = [L;l];
                    S = [S;s];
                    xx(a,a) = 0;
                    pp(a)=0;
                    a=a+1;
                end
            end
        end
    end
    ZZ = [L,S];
    for k=1:K-1
        Q((k-1)*M*M*N+1:(k)*M*M*N,(k-1)*M*N+1:(k+1)*M*N) = ZZ;
    end
    
    y = zeros(N*M*M*(K-1),N*M);
    
    
    Q = [y,Q];
    xxx = [];
    ppp = [];
    for k = 1:K-1
        xxx = [xxx;xx];
        ppp = [ppp;pp'];
    end
    Q = [Q,-1*xxx];
    X = [X;Q];
    
    b = [b;ppp];
    
    
    Xeq = [];
    Qeq = [];
    for k = 1:K+1
        P = zeros(N,N*M);                  %7 8 9 10 11
        Y = [];
        for i=1:N
            y = zeros(N,M);
            for j=1:M
                if w(i,j)~=0
                    y(i,j) = 1;
                end
            end
            y = y';
            y = y(:);
            Y = [Y;y'];
        end
        P(:,:,k) = Y;
        Qeq((k-1)*N+1:k*N,:,k) = P(:,:,k);
    end
    Xeq = reshape(Qeq,(K+1)*N,(K+1)*N*M);
    beq = K*ones(N,1);
    beq = [beq;1*ones(K*N,1)];
    
    
    QQ = [];                             %12
    for i=1:N
        for j=1:M
            if w(i,j)~=0
                qq = zeros(N,M);
                qq(i,j) = 1;
                qq = qq';
                qq = qq(:);
                QQ = [QQ;qq'];
            else
                qq = zeros(1,N*M);
                QQ = [QQ;qq];
            end
        end
    end
    RR = [];
    for k = 1:K
        RR = [RR,QQ];
    end
    RR = [-1*QQ,RR];
    Xeq = [Xeq;RR];
    beq = [beq;zeros(N*M,1)];
    
    
    Xeq = [Xeq,zeros(N*(K+1)+N*M,N*M*M)];
    
    
    e = 0.001;
    alpha = 1000;
    [rows1,~] = size(X);
    X = [X,zeros(rows1,N*M*M)];
    
    [rows2,~] = size(Xeq);
    Xeq = [Xeq,zeros(rows2,N*M*M)];
    
    
    a = 1;
    ww1 = zeros(1,N*M*M);
    ww2 = zeros(1,N*M*M);
    sl = zeros(N*M*M);
    for i = 1:N
        for j = 1:M
            for jdash = 1:M
                if j~=jdash && w(i,j)~=0 && w(i,jdash)~=0
                    ww1(a) = w(i,j);
                    ww2(a) = w(i,jdash);
                    sl(a,a) = 1;
                    epsilon(a) = e;
                    a = a+1;                                                                                                            kks=0;
                else
                    ww1(a) = 0;
                    ww2(a) = 0;
                    sl(a,a) = 0;
                    epsilon(a) = 0;
                    a = a+1;
                end
            end
        end
    end
    xx2 = N*xx;
    lhs1 = [zeros(N*M*M,N*M*(K+1)),xx2,-1*sl];
    rhs1 = -1*ww1'+ww2'+N;
    X = [X;lhs1];
    b = [b;rhs1];
    
    lhs2 = [zeros(N*M*M,N*M*(K+1)),-1*xx2,sl];
    rhs2 = ww1'-ww2'-epsilon';
    
    X = [X;lhs2];
    b = [b;rhs2];
    
    intcon = (1:((K+1)*N*M+N*M*M));
    lb = zeros(1,((K+1)*N*M)+2*N*M*M);
    ub = ones(1,((K+1)*N*M+N*M*M));
    if ~stopit
        options = optimoptions('intlinprog','IntegerTolerance',1e-6,'MaxTime',str2num(handles.max_run_time.String),'OutputFcn',@stopmyiter);
        [c,fval,exitflag] = intlinprog(f,intcon,X,b,Xeq,beq,lb,ub,options);
    else
        set(handles.termination_state,'String','Intlinprog has stopped');
        drawnow;
    end
    switch exitflag
        case 2
            exitval = 'intlinprog stopped prematurely. Integer feasible point found.';
        case 1
            exitval = 'intlinprog converged to the solution x.';
        case 0
            exitval = 'intlinprog stopped prematurely. No integer feasible point found.';
        case -1
            exitval = 'intlinprog stopped by an output function or plot function.';
        case -2
            exitval = 'No feasible point found. Please enter alpha and beta';
        case -3
            exitval = 'Root LP problem is unbounded.';
    end
    
    set(handles.termination_state,'String',exitval);
    if ~stopit
        set(handles.table1,'Visible','on');
        drawnow;
        set(handles.table2,'Visible','on');
        drawnow;
        set(handles.edit11,'Visible','on');
        set(handles.text9,'Visible','on');
        set(handles.edit11,'enable','on','String',num2str(fval));
        
        kkk = reshape(c(1:(K+1)*N*M),[N*M,(K+1)]);
        val = kkk';
        
        final1 = zeros(N,M);
        for k=1:K
            vec = reshape(val(k+1,:),[M,N]);
            vec = k*vec';
            final1 = final1+vec;
        end
        final1 = round(final1);
        
        
        
        proposal_num = {};
        for i=1:N
            proposal_num{end+1} =strcat('Prop',num2str(i));
        end
        
        reviewer_num = {};
        for i=1:M
            reviewer_num{end+1} =strcat('rev',num2str(i));
        end
        
        handles.table1.RowName = proposal_num;
        handles.table1.ColumnName = reviewer_num;
        handles.table1.Data = final1;
        
        final2 = [];
        for k=1:K
            final2 = [final2;sum(final1==k)];
        end
        final2 = [final2;sum(final2)];
        rown = {'Lead';'Scribe'};
        for k=1:K-2
            rown{end+1} = strcat('R',num2str(k));
        end
        rown{end+1} = 'Total';
        handles.table2.RowName = rown;
        handles.table2.ColumnName = reviewer_num;
        handles.table2.Data = final2;
    end
else
    if ~stopit
        set(handles.table1,'Visible','on');
        drawnow;
        set(handles.table2,'Visible','on');
        drawnow;
        set(handles.edit11,'Visible','on');
        set(handles.text9,'Visible','on');
        set(handles.edit11,'enable','on','String',num2str(fval));
        
        kkk = reshape(c(1:(K+1)*N*M),[N*M,(K+1)]);
        val = kkk';
        
        final1 = zeros(N,M);
        for k=1:K
            vec = reshape(val(k+1,:),[M,N]);
            vec = k*vec';
            final1 = final1+vec;
        end
        final1 = round(final1);
        
        
        
        proposal_num = {};
        for i=1:N
            proposal_num{end+1} =strcat('Prop',num2str(i));
        end
        
        reviewer_num = {};
        for i=1:M
            reviewer_num{end+1} =strcat('rev',num2str(i));
        end
        
        handles.table1.RowName = proposal_num;
        handles.table1.ColumnName = reviewer_num;
        handles.table1.Data = final1;
        
        final2 = [];
        for k=1:K
            final2 = [final2;sum(final1==k)];
        end
        final2 = [final2;sum(final2)];
        rown = {'Lead';'Scribe'};
        for k=1:K-2
            rown{end+1} = strcat('R',num2str(k));
        end
        rown{end+1} = 'Total';
        handles.table2.RowName = rown;
        handles.table2.ColumnName = reviewer_num;
        handles.table2.Data = final2;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DirPath
check_save = handles.save_check.Value;
if check_save ==1
    clk = clock;
    name = [DirPath '/' num2str(clk(4)) num2str(clk(5)) num2str(ceil(clk(6))) '.mat'];
    save(name,'c','fval');  
end  
set(handles.reset_push,'Enable','on');
set(handles.stop_push,'Enable','off');


% --- Executes on button press in reset_push.
function reset_push_Callback(hObject, eventdata, handles)
% hObject    handle to reset_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stopit
set(handles.N_val,'enable','on','String','');
set(handles.M_val,'enable','on','String','');
set(handles.K_val,'enable','on','String','');
set(handles.w_val,'enable','on','String','');
set(handles.max_run_time,'enable','on','String','7200');
set(handles.alpha_val,'enable','on','String','1000');
set(handles.beta_val,'enable','on','String','0.001');
set(handles.save_check,'enable','on');
set(handles.save_check,'Value',0);
set(handles.save_dir,'Visible','off');
set(handles.browse_push,'Visible','off');
set(handles.termination_state,'Visible','off');
set(handles.table1,'Visible','off');
set(handles.table2,'Visible','off');
set(handles.edit11,'Visible','off');
set(handles.text9,'Visible','off');


stopit = 0;

% --- Executes on button press in stop_push.
function stop_push_Callback(hObject, eventdata, handles)
% hObject    handle to stop_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stopit
stopit=1;
set(handles.reset_push,'Enable','on');


function termination_state_Callback(hObject, eventdata, handles)
% hObject    handle to termination_state (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of termination_state as text
%        str2double(get(hObject,'String')) returns contents of termination_state as a double


% --- Executes during object creation, after setting all properties.
function termination_state_CreateFcn(hObject, eventdata, handles)
% hObject    handle to termination_state (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_check.
function save_check_Callback(hObject, eventdata, handles)
% hObject    handle to save_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_check
if handles.save_check.Value
    set(handles.save_dir,'Visible','on');
    set(handles.browse_push,'Visible','on');
else
    set(handles.save_dir,'Visible','off')
    set(handles.browse_push,'Visible','off');
end


function save_dir_Callback(hObject, eventdata, handles)
% hObject    handle to save_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of save_dir as text
%        str2double(get(hObject,'String')) returns contents of save_dir as a double


% --- Executes during object creation, after setting all properties.
function save_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_push.
function browse_push_Callback(hObject, eventdata, handles)
% hObject    handle to browse_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DirPath
path=uigetdir();
if any(path)~=0
    DirPath = path;
    set(handles.save_dir,'String',DirPath);
else
    DirPath = pwd;
end



function max_run_time_Callback(hObject, eventdata, handles)
% hObject    handle to max_run_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_run_time as text
%        str2double(get(hObject,'String')) returns contents of max_run_time as a double


% --- Executes during object creation, after setting all properties.
function max_run_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_run_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
