function varargout = chsel_gui(varargin)
% CHSEL_GUI MATLAB code for chsel_gui.fig
%      CHSEL_GUI, by itself, creates a new CHSEL_GUI or raises the existing
%      singleton*.
%
%      H = CHSEL_GUI returns the handle to a new CHSEL_GUI or the handle to
%      the existing singleton*.
%
%      CHSEL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHSEL_GUI.M with the given input arguments.
%
%      CHSEL_GUI('Property','Value',...) creates a new CHSEL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before chsel_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to chsel_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help chsel_gui

% Last Modified by GUIDE v2.5 27-Oct-2015 15:10:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @chsel_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @chsel_gui_OutputFcn, ...
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


% --- Executes just before chsel_gui is made visible.
function chsel_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to chsel_gui (see VARARGIN)

% Choose default command line output for chsel_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Load data from main GUI
data = getappdata(0,'data');

if ~isfield(data,'call') || size(data.sig,2)~=size(data.sig_rough,2) % if not previously saved results or if # channels different
    data.deci_len = 1000;
    data.sig_rough = getroughsig(data.sig,data.deci_len);
    data.sig_rough_t = (0:size(data.sig_rough,1)-1)/data.fs*data.deci_len;
    data.shift_gap = max(max(squeeze(data.sig_rough)));
    data.shift_gap = ceil(data.shift_gap/0.01)*0.01;
end

% Plot all channels
axes(handles.axes_ch_sig);
plot(data.sig_rough_t,data.sig_rough+repmat((1:size(data.sig_rough,2))*data.shift_gap,size(data.sig_rough,1),1));
xlim(data.sig_rough_t([1 end]));
set(gca,'ytick',(1:data.num_ch_in_file)*data.shift_gap,'yticklabel',{num2str((1:data.num_ch_in_file)')});
ylim([0 (data.num_ch_in_file+1)*data.shift_gap]);

% Set default channel number
set(handles.edit_chsel,'String',num2str(round(data.num_ch_in_file/2)));

% Save processed stuff back to global data
setappdata(0,'data2',data);

% If only one channel or calls have been detected --> send channel selection directly
if data.num_ch_in_file==1 || isfield(data,'call')
     button_select_Callback(hObject, eventdata, handles);
     % go back to main GUI
else  % enable the channel selection button
    set(handles.button_select,'enable','on');
    set(handles.edit_chsel,'enable','on');
    set(handles.text1,'enable','on');
end

% change GUI window name to loaded file
set(handles.figure1,'Name',data.fname);

% UIWAIT makes chsel_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function sig_rough = getroughsig(sig,deci_len)
sig_rough = zeros((size(sig,1)-mod(size(sig,1),deci_len))/deci_len,size(sig,2));
for iCH = 1:size(sig,2)
    sig_ch = sig(:,iCH);
    sig_ch = reshape(sig_ch(1:end-mod(length(sig_ch),deci_len)),deci_len,[]);
    sig_rough(:,iCH) = max(sig_ch,[],1);
end


% --- Outputs from this function are returned to the command line.
function varargout = chsel_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_chsel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_chsel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% save global data
data = getappdata(0,'data2');
data.chsel = str2num(get(handles.edit_chsel,'String'));
gui_op = getappdata(0,'gui_op');
gui_op.chsel_current = data.chsel;  % keep track of the current channel being displayed

% disable button and edit
set(handles.button_select,'enable','off');
set(handles.edit_chsel,'enable','off');
set(handles.text1,'enable','off');

% save global data
setappdata(0,'data',data);
setappdata(0,'gui_op',gui_op);
setappdata(0,'handles_ch_sig',handles);

% go back to main GUI
hh = getappdata(0,'handles_op');
axes(hh.axes_spectrogram);

% close(chsel_gui);


% Hints: get(hObject,'String') returns contents of edit_chsel as text
%        str2double(get(hObject,'String')) returns contents of edit_chsel as a double


% --- Executes during object creation, after setting all properties.
function edit_chsel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_chsel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_select.
function button_select_Callback(hObject, eventdata, handles)
% hObject    handle to button_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% save global data
data = getappdata(0,'data2');
data.chsel = str2num(get(handles.edit_chsel,'String'));
gui_op = getappdata(0,'gui_op');
gui_op.chsel_current = data.chsel;  % keep track of the current channel being displayed

% disable button and edit
set(handles.button_select,'enable','off');
set(handles.edit_chsel,'enable','off');
set(handles.text1,'enable','off');

% save global data
setappdata(0,'data',data);
setappdata(0,'gui_op',gui_op);
setappdata(0,'handles_ch_sig',handles);

% go back to main GUI
hh = getappdata(0,'handles_op');
axes(hh.axes_spectrogram);

% close(chsel_gui);
