function varargout = call_duration_gui(varargin)
% CALL_DURATION_GUI MATLAB code for call_duration_gui.fig
%      CALL_DURATION_GUI, by itself, creates a new CALL_DURATION_GUI or raises the existing
%      singleton*.
%
%      H = CALL_DURATION_GUI returns the handle to a new CALL_DURATION_GUI or the handle to
%      the existing singleton*.
%
%      CALL_DURATION_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALL_DURATION_GUI.M with the given input arguments.
%
%      CALL_DURATION_GUI('Property','Value',...) creates a new CALL_DURATION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before call_duration_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to call_duration_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help call_duration_gui

% Wu-Jung Lee | leewujung@gmail.com
% 2015/12/21  Code overhaul using guide, use new file format
%             (separate data and detection files)


% Last Modified by GUIDE v2.5 22-Dec-2015 10:32:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @call_duration_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @call_duration_gui_OutputFcn, ...
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


% --- Executes just before call_duration_gui is made visible.
function call_duration_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to call_duration_gui (see VARARGIN)

% Choose default command line output for call_duration_gui
handles.output = hObject;

% Link plotting axes
linkaxes([handles.axes_spectrogram,handles.axes_time_series],'x');

% Set zoom and pan motion
gui_op.hzoom = zoom;
gui_op.hpan = pan;

% Update handles structure
guidata(hObject, handles);

% Save global var
data = [];
setappdata(0,'data_duration',data);
setappdata(0,'gui_op_duration',gui_op);

% UIWAIT makes call_duration_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = call_duration_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_load_sig.
function button_load_sig_Callback(hObject, eventdata, handles)
% hObject    handle to button_load_sig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Set path
gui_op = getappdata(0,'gui_op_duration');
if ~isfield(gui_op,'sig_path')
    if ispref('call_duration_gui') && ispref('call_duration_gui','sig_path')
        sig_path = getpref('call_duration_gui','sig_path');
    else
        sig_path = '';
    end
    gui_op.sig_path = uigetdir(sig_path,'Select the folder containing *_mic_data.mat files');
    if isequal(gui_op.sig_path,0)
        return
    end
    setpref('call_duration_gui','sig_path',gui_op.sig_path);
end
if ~isfield(gui_op,'duration_path')
    if ispref('call_duration_gui') && ispref('call_duration_gui','duration_path')
        duration_path = getpref('call_duration_gui','duration_path');
    else
        duration_path = gui_op.sig_path;
    end
    gui_op.duration_path = uigetdir(duration_path,'Select the folder containing *_detect.mat files');
    if isequal(gui_op.duration_path,0)
        return
    end
    setpref('call_duration_gui','duration_path',gui_op.duration_path);
end

% Select file to load
[fname,pname] = uigetfile(fullfile(gui_op.sig_path,'*.mat'),'Select signal file');
if isequal(pname,0)
    return
else
    k = strfind(fname,'_detect');
    if isempty(k)  % if not '_detect' file
        ss = strsplit(fname,'.');
        fname_det_pre = strjoin(ss(1:length(ss)-1),'.');
        fname_det = [fname_det_pre,'_detect.mat'];
        fname_sig = fname;
    else
        fname_det = fname;
        fname_sig = fname(1:k-1);
    end
    gui_op.fname_sig = fname_sig;
    gui_op.fname_det = fname_det;
end

setappdata(0,'gui_op_duration',gui_op);


% Load detection results
if exist(fullfile(gui_op.duration_path,fname_det),'file')  % if detect file exists in current path
    D_det = load(fullfile(gui_op.duration_path,fname_det));
    % copy all fields into current data structure
    field_in_file = fieldnames(D_det);
    for iF=1:length(field_in_file)
        data.(field_in_file{iF}) = D_det.(field_in_file{iF});
    end
    disp('Previous call marking results loaded');
else
    disp('No previous call marking results');
end

% Load mic signal
if exist(fullfile(gui_op.sig_path,fname_sig),'file')  % if sig file exists in current path
    D_sig = load(fullfile(gui_op.sig_path,fname_sig));
    data.sig = D_sig.sig;  % only load sig from the raw file
    disp('Signal loaded');
else
    disp('Mic file not found');
    return
end

% patch for num_ch_in_file
data.num_ch_in_file = size(data.sig,2);
disp(['Number of channels in file: ',num2str(data.num_ch_in_file)]);

% patch for sig_t
if ~isfield(data,'sig_t')
    data.sig_t = (0:length(data.sig)-1)/data.fs;
else
    if isempty(data.sig_t)
        data.sig_t = (0:length(data.sig)-1)/data.fs;
    end
end

data.fname = fname_sig;
data.pname = pname;

% change GUI window name to loaded file
set(handles.figure1,'Name',data.fname);

% gui_op stuff
if isfield(data,'call') % if previously saved results
    gui_op.curr_ch = data.call(1).channel_marked;  % keep track of the current channel being displayed
end
gui_op.curr_call_num = 1;
gui_op.curr_call_crange = [NaN NaN]; % keep track of current coloraxis range

% update displays
set(handles.edit_curr_ch,'String',num2str(gui_op.curr_ch));
set(handles.edit_curr_call_num,'String',num2str(gui_op.curr_call_num));
set(handles.text_total_call_num,'String',['/',num2str(length(data.call))]);

setappdata(0,'data_duration',data);
setappdata(0,'gui_op_duration',gui_op);

% Clear current display
cla(handles.axes_spectrogram);
cla(handles.axes_time_series);

% Plot call
showcall(handles)


% --- Executes on button press in button_next.
function button_next_Callback(hObject, eventdata, handles)
% hObject    handle to button_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data_duration');
gui_op = getappdata(0,'gui_op_duration');
tmp = gui_op.curr_call_num+1;  % update current call index
if tmp<=length(data.call)
    data.call(gui_op.curr_call_num).channel_marked = str2double(get(handles.edit_curr_ch,'String'));  % save the channel used for duration marking
    gui_op.curr_call_num = tmp;  % increment call number
    gui_op.curr_ch = data.call(gui_op.curr_call_num).channel_marked;  % set current display channel
    
    % save data
    setappdata(0,'gui_op_duration',gui_op);
    setappdata(0,'data_duration',data);
    
    % update displays
    set(handles.checkbox_low_quality,'Val',data.call(gui_op.curr_call_num).low_quality);
    set(handles.edit_curr_call_num,'String',num2str(gui_op.curr_call_num));
    set(handles.edit_curr_ch,'String',num2str(gui_op.curr_ch));
    showcall(handles);
else
    msgbox('All calls processed!');
end

% --- Executes on button press in button_previous.
function button_previous_Callback(hObject, eventdata, handles)
% hObject    handle to button_previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data_duration');
gui_op = getappdata(0,'gui_op_duration');
tmp = gui_op.curr_call_num-1;  % update current call index
if tmp>=1
    data.call(gui_op.curr_call_num).channel_marked = str2double(get(handles.edit_curr_ch,'String'));  % save the channel used for duration marking
    gui_op.curr_call_num = tmp;
    gui_op.curr_ch = data.call(gui_op.curr_call_num).channel_marked;  % set current display channel
    
    % save data
    setappdata(0,'gui_op_duration',gui_op);
    setappdata(0,'data_duration',data);
    
    % update displays
    set(handles.checkbox_low_quality,'Val',data.call(gui_op.curr_call_num).low_quality);
    set(handles.edit_curr_call_num,'String',num2str(gui_op.curr_call_num));
    set(handles.edit_curr_ch,'String',num2str(gui_op.curr_ch));
    showcall(handles);
else
    msgbox('At the beginning of call sequence!');
end


% --- Executes on button press in button_delete_call.
function button_delete_call_Callback(hObject, eventdata, handles)
% hObject    handle to button_delete_call (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data_duration');
gui_op = getappdata(0,'gui_op_duration');

button = questdlg('Confirm delete this call?','Delete this call?','Yes','No','No');
if strcmp(button,'Yes')
    data.call(gui_op.curr_call_num) = [];
    data.aux_data(gui_op.curr_call_num) = [];
    data.total_call_num = length(data.call);
    if gui_op.curr_call_num>data.total_call_num
        gui_op.curr_call_num = data.total_call_num;
    end
    set(handles.text_total_call_num,'String',['/',num2str(data.total_call_num)]);
    set(handles.edit_curr_call_num,'String',num2str(gui_op.curr_call_num));
    setappdata(0,'data_duration',data);
    setappdata(0,'gui_op_duration',gui_op);
    showcall(handles);
end


% --- Executes on button press in button_ch_plus.
function button_ch_plus_Callback(hObject, eventdata, handles)
% hObject    handle to button_ch_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data_duration');
gui_op = getappdata(0,'gui_op_duration');
tmp = mod(gui_op.curr_ch+1,data.num_ch_in_file);
if tmp==0
    gui_op.curr_ch = data.num_ch_in_file;
else
    gui_op.curr_ch = tmp;
end
setappdata(0,'gui_op_duration',gui_op);
set(handles.edit_curr_ch,'String',num2str(gui_op.curr_ch));
showcall(handles);


% --- Executes on button press in button_ch_minus.
function button_ch_minus_Callback(hObject, eventdata, handles)
% hObject    handle to button_ch_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data_duration');
gui_op = getappdata(0,'gui_op_duration');
tmp = mod(gui_op.curr_ch-1,data.num_ch_in_file);
if tmp==0
    gui_op.curr_ch = data.num_ch_in_file;
else
    gui_op.curr_ch = tmp;
end
setappdata(0,'gui_op_duration',gui_op);
set(handles.edit_curr_ch,'String',num2str(gui_op.curr_ch));
showcall(handles);


% --- Executes on button press in button_save.
function button_save_Callback(hObject, eventdata, handles)
% hObject    handle to button_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = getappdata(0,'data_duration');
gui_op = getappdata(0,'gui_op_duration');

% prepare for saving results
field_in_file = fieldnames(data);
for iF=1:length(field_in_file)
    A.(field_in_file{iF}) = data.(field_in_file{iF});
end
A.sig = [];

tt = strsplit(A.fname,'.mat');
save_fname = sprintf('%s_detect.mat',tt{1});
[save_fname,save_pname] = uiputfile('*.mat','Save detection results',fullfile(gui_op.duration_path,save_fname));
if isequal(save_pname,0)
    return
else
    save([save_pname,'/',save_fname],'-struct','A');
end



function edit_curr_ch_Callback(hObject, eventdata, handles)
% hObject    handle to edit_curr_ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_curr_ch as text
%        str2double(get(hObject,'String')) returns contents of edit_curr_ch as a double



function edit_curr_call_num_Callback(hObject, eventdata, handles)
% hObject    handle to edit_curr_call_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data_duration');
gui_op = getappdata(0,'gui_op_duration');
tmp = str2double(get(handles.edit_curr_call_num,'String'));  % updated current call index
if tmp>length(data.call)
    tmp = length(data.call);
elseif tmp<1
    tmp = 1;
end
data.call(gui_op.curr_call_num).channel_marked = str2double(get(handles.edit_curr_ch,'String'));  % save the channel used for duration marking
gui_op.curr_call_num = tmp;  % increment call number
gui_op.curr_ch = data.call(gui_op.curr_call_num).channel_marked;  % set current display channel

% save data
setappdata(0,'gui_op_duration',gui_op);
setappdata(0,'data_duration',data);

% update displays
set(handles.checkbox_low_quality,'Val',data.call(gui_op.curr_call_num).low_quality);
set(handles.edit_curr_call_num,'String',num2str(gui_op.curr_call_num));
set(handles.edit_curr_ch,'String',num2str(gui_op.curr_ch));
showcall(handles);


% --- Executes during object creation, after setting all properties.
function edit_curr_call_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_curr_call_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text_total_call_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_total_call_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Plot spectrogram and time series of a single call

% --- Executes on button press in checkbox_low_quality.
function checkbox_low_quality_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_low_quality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data_duration');
gui_op = getappdata(0,'gui_op_duration');
data.call(gui_op.curr_call_num).low_quality = get(handles.checkbox_low_quality,'Val');
setappdata(0,'data_duration',data);


% --- Executes on slider movement.
function slider_caxis_Callback(hObject, eventdata, handles)
% hObject    handle to slider_caxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_op = getappdata(0,'gui_op_duration');
caxis(handles.axes_spectrogram,[gui_op.curr_call_crange(1)+...
             range(gui_op.curr_call_crange)*get(handles.slider_caxis,'value'),...
             gui_op.curr_call_crange(2)]);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_caxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_caxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_call_disp_len_Callback(hObject, eventdata, handles)
% hObject    handle to edit_call_disp_len (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
showcall(handles);
% Hints: get(hObject,'String') returns contents of edit_call_disp_len as text
%        str2double(get(hObject,'String')) returns contents of edit_call_disp_len as a double


% --- Executes during object creation, after setting all properties.
function edit_call_disp_len_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_call_disp_len (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% delete global var

if isappdata(0,'data_duration')
    rmappdata(0,'data_duration');
end
if isappdata(0,'gui_op_duration')
    rmappdata(0,'gui_op_duration');
end

% Hint: delete(hObject) closes the figure
delete(hObject);



function showcall(handles)
data = getappdata(0,'data_duration');
gui_op = getappdata(0,'gui_op_duration');

% extract data for the current call
ch_sel = gui_op.curr_ch;
call_loc = data.call(gui_op.curr_call_num).locs;
ext_time = str2double(get(handles.edit_call_disp_len,'String'))/1e3;  % length of signal to display [sec]
want_idx = call_loc + (-round(ext_time*data.fs)/2:round(ext_time*data.fs)/2);
call = data.sig(want_idx,ch_sel);
call_start_idx = data.call(gui_op.curr_call_num).call_start_idx-want_idx(1)+1;
call_end_idx = data.call(gui_op.curr_call_num).call_end_idx-want_idx(1)+1;
if isnan(call_start_idx)
    call_start_idx = round(1/3*length(want_idx));
end
if isnan(call_end_idx)
    call_end_idx = round(2/3*length(want_idx));
end
call_time = (1:length(call))/data.fs; % time stamp of call [s]

gui_op.want_idx_start = want_idx(1);

% plot spectrogram
[~,F,T,P] = spectrogram(call,128,120,128,data.fs);
axes(handles.axes_spectrogram);
P = 10*log10(abs(P));  % [dB]
imagesc(T*1e3,F/1e3,P);
if isnan(gui_op.curr_call_crange(1))
    gui_op.curr_call_crange = [min(min(P)) max(max(P))];
end
caxis(handles.axes_spectrogram,[gui_op.curr_call_crange(1)+...
             range(gui_op.curr_call_crange)*get(handles.slider_caxis,'value'),...
             gui_op.curr_call_crange(2)]);
hold on
axis xy
ylabel('Frequency (kHz)');
setAxesPanMotion(gui_op.hpan,handles.axes_spectrogram,'horizontal');
setAxesZoomMotion(gui_op.hzoom,handles.axes_spectrogram,'horizontal');

% plot time-series
axes(handles.axes_time_series);
plot(call_time*1e3,call);
call_volt_range = range(call);
yy = [min(call)-call_volt_range*0.1, max(call)+call_volt_range*0.1];
xlim(T([1 end])*1e3);
ylim(yy);
hold on
xlabel('Time (ms)');
ylabel('Volt');
setAxesPanMotion(gui_op.hpan,handles.axes_time_series,'horizontal');
setAxesZoomMotion(gui_op.hzoom,handles.axes_time_series,'horizontal');

% plot vertical lines on time-series
gui_op.h2start = plot(handles.axes_time_series,[1 1]*call_start_idx/data.fs*1e3,yy,...
                  'm-','linewidth',2,'tag','call_start_line');
gui_op.h2end = plot(handles.axes_time_series,[1 1]*call_end_idx/data.fs*1e3,yy,...
                  'm-','linewidth',2,'tag','call_end_line');

% plot box for current start/end and bandwidth estimation on spectrogram
call_idx_len = call_end_idx-call_start_idx+1;
if isnan(data.call(gui_op.curr_call_num).bandwidth(1))
    data.call(gui_op.curr_call_num).bandwidth = [10 80]*1e3;  % initialize bandwidth estimation
end
bw_len = diff(data.call(gui_op.curr_call_num).bandwidth);
gui_op.hbox = imrect(handles.axes_spectrogram,[call_start_idx/data.fs*1e3,...  % convert to [ms] for plotting
                        data.call(gui_op.curr_call_num).bandwidth(1)/1e3,...  % convert to [kHz] for plotting
                        call_idx_len/data.fs*1e3,...  % convert to [ms] for plotting
                        bw_len/1e3]);  % convert to [kHz] for plotting

% gui_op.hbox.addNewPositionCallback(@(pos)saveNewPos(pos,gcbo));
gui_op.hbox.addNewPositionCallback(@(pos)saveNewPos(pos));

hold(handles.axes_spectrogram,'off');
hold(handles.axes_time_series,'off');

% save data
setappdata(0,'data_duration',data);
setappdata(0,'gui_op_duration',gui_op);


function saveNewPos(pos)
data = getappdata(0,'data_duration');
gui_op = getappdata(0,'gui_op_duration');

% move vertical lines on the time-series plot
xnew = sort([pos(1) pos(1)+pos(3)]);
set(gui_op.h2start,'XData',xnew(1)*[1 1]);
set(gui_op.h2end,'XData',xnew(2)*[1 1]);

% save new pos to object
xx = sort([pos(1) (pos(1)+pos(3))]);
data.call(gui_op.curr_call_num).call_start_idx = round(xx(1)/1e3*data.fs)+gui_op.want_idx_start-1;
data.call(gui_op.curr_call_num).call_end_idx = round(xx(2)/1e3*data.fs)+gui_op.want_idx_start-1;
data.call(gui_op.curr_call_num).bandwidth = sort([pos(2) (pos(2)+pos(4))])*1e3;  % call start and end freq [Hz]

% save data
setappdata(0,'data_duration',data);
setappdata(0,'gui_op_duration',gui_op);
