function varargout = NeuroSeg(varargin)
% NeuroSeg MATLAB code for NeuroSeg.fig
%      NeuroSeg, by itself, creates a new NeuroSeg or raises the existing
%      singleton*.
%      H = NeuroSeg returns the handle to a new NeuroSeg or the handle to
%      the existing singleton*.
%
%      NeuroSeg('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NeuroSeg.M with the given input arguments.
%
%      NeuroSeg('Property','Value',...) creates a new NeuroSeg or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NeuroSeg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NeuroSeg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NeuroSeg

% Last Modified by GUIDE v2.5 13-Dec-2016 20:17:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NeuroSeg_OpeningFcn, ...
                   'gui_OutputFcn',  @NeuroSeg_OutputFcn, ...
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


% --- Executes just before NeuroSeg is made visible.
function NeuroSeg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NeuroSeg (see VARARGIN)

% Choose default command line output for NeuroSeg
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
set(handles.RoI_index,'Max',2);
axes(handles.axes1);
imlog=imread('logo.jpg');
imshow(imlog);
text(0.5*size(imlog,1)-200,0.5*size(imlog,2)-150,'Welcome to NeuroSeg','FontSize',30);

handles.thetaN=8;
set(handles.Frequency,'string',num2str(20));
set(handles.Interval,'string',num2str(2));
set(handles.ROIradius,'string',num2str(5));
set(handles.MinSigma,'string',num2str(9));%9   11
set(handles.MaxSigma,'string',num2str(16));%16  26
set(handles.MinArea,'string',num2str(450));%450 800

set(gca, 'XTick', []);
set(gca, 'YTick', []);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = NeuroSeg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in LoadPicture_pushbutton.
function LoadPicture_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadPicture_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
axes(handles.axes1);
cla;
[filename pathname]=uigetfile('*.tif;*.bmp;*.jpeg','Please select a picture file');
[I_ori,map,alpha]=imread([pathname filename]);
if size(I_ori,3)>1
   I = mean(I_ori,3);
else
   I = I_ori;
end

if isfield(handles,'SeedsScatter')
   handles=rmfield(handles,'SeedsScatter');
end
if isfield(handles,'htext')
   handles=rmfield(handles,'htext');
end
if isfield(handles,'hline') 
    handles=rmfield(handles,'hline');
end

I_adjust=double(I);
I_adjust=(I_adjust-min(I_adjust(:)))/(max(I_adjust(:))-min(I_adjust(:)));
if size(I_ori,3)>1
    imshow(I_ori(:,:,1:3));
else
    imagesc(I_adjust);
    colormap(gray);
    axis equal;
    axis off;
end
hold on;
handles.I_adjust=I_adjust; %imagedata between 0 and 1;
handles.ImageAdjust=uint8(I_adjust*255);%image data between 0 and 255
handles.filename=filename;
handles.pathname=pathname;
handles.IMGxrange=size(I,2);
handles.IMGyrange=size(I,1);
set(handles.edit1,'string',[pathname filename]);
set(handles.edit1,'ForegroundColor',[1 0 0]);
set(handles.RoI_index,'string',[])
handles.parasaved=[];
guidata(hObject,handles);


% --- Executes on button press in LoadRoiAndPicturepushbutton.
function LoadRoiAndPicturepushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadRoiAndPicturepushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
axes(handles.axes1);
cla;
[filename pathname] = uigetfile('*.mat','Select the MAT file for the Image and ROI');
load([pathname filename]);

set(handles.edit1,'string',[pathname filename]);
set(handles.edit1,'ForegroundColor',[1 0 0]);%statusbar(0, 'Desktop status: processing...');

handles.parasaved = ParametersOutput;
I=ImageDataMean;
if size(I,3)>1
   I = rgb2gray(I);
end
I_adjust=double(I);
I_adjust=(I_adjust-min(I_adjust(:)))/(max(I_adjust(:))-min(I_adjust(:)));

if isfield(handles,'SeedsScatter')
   delete(handles.SeedsScatter)%delete the scatter
   handles=rmfield(handles,'SeedsScatter');
end
if isfield(handles,'htext')
   delete(handles.htext);
   handles=rmfield(handles,'htext');
end
if isfield(handles,'hline')
    delete(handles.hline); 
    handles=rmfield(handles,'hline');
end
imagesc(I_adjust);
colormap(gray);
axis equal;
axis off;
hold on;

hline=[];
htext=[];
if iscell(ParametersOutput.xy_all)==1
    ParametersOutput.xy_all=cell2mat(ParametersOutput.xy_all);
end
for k=1:length(ParametersOutput.xypos)
    hline(1,k)=plot(ParametersOutput.xypos{k}(1,:),ParametersOutput.xypos{k}(2,:),'r','LineWidth',2);
    htext(1,k)=text(ParametersOutput.xy_all(1,k),ParametersOutput.xy_all(2,k),num2str(k-1),'color','g','FontSize',12);
end
handles.ImageAdjust=ImageDataMean;%save the image data
handles.htext = htext;
handles.hline = hline;

handles.IMGxrange=size(I,2);
handles.IMGyrange=size(I,1);

cell_no=(0:length(ParametersOutput.xypos)-1);
%set(handles.RoI_index,'Value',[1:10]')
set(handles.RoI_index,'String',num2str(cell_no'),'Value',1)

handles.filename=filename;
handles.pathname=pathname;

guidata(hObject,handles) ;
axes(handles.axes1);


% --- Executes on button press in Save_pushbutton4.
function Save_pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to Save_pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
axes(handles.axes1);
htext = handles.htext;
hline = handles.hline;
ImageDataMean = handles.ImageAdjust;
ParametersOutput=handles.parasaved;
filename = handles.filename;
[FileName,PathName] = uiputfile([filename(1:end-4) '_Image_ROIdata.mat'],'Save Image and ROI');
save([PathName '\' FileName],'ImageDataMean','ParametersOutput','hline','htext');

h=figure('visible','off');
imshow(ImageDataMean);
hold on;
for k=1:length(ParametersOutput.xypos)
    plot(ParametersOutput.xypos{k}(1,:),ParametersOutput.xypos{k}(2,:),'r','LineWidth',2);
    text(ParametersOutput.xy_all(1,k),ParametersOutput.xy_all(2,k),num2str(k-1),'color','g','FontSize',12);
end
title('Average image with ROIs','FontSize',30);
set(gcf,'PaperPositionMode','auto');
set(gcf, 'position', get(0,'ScreenSize'));
print(h,'-dtiff','-r0',[PathName '\' FileName(1:end-4) '_Image&ROI.tif']);
close(h);


% --- Executes on button press in NeuroSeg_pushbutton.
function DrawRoI_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to NeuroSeg_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
axes(handles.axes1);
if ~isfield(handles,'ImageAdjust')
     msgbox('Please load a picure file or a picture&RoI file') ;
else
hf = imfreehand; 
pos = getPosition(hf);
delete(hf);
xc = [pos(:,1);pos(1,1)];
yc = [pos(:,2);pos(1,2)];
xc(xc>handles.IMGxrange)=handles.IMGxrange;
yc(yc>handles.IMGyrange)=handles.IMGyrange;
xc(xc<1)=1;
yc(yc<1)=1;   
A = polyarea(xc,yc);
if A>=1
xypos=[xc';yc'];%xypos
[axisX axisY]=meshgrid(round(min(xc)):round(max(xc)),round(min(yc)):round(max(yc)));
axis_data=[axisX(:) axisY(:)];
in_image = inpolygon(axis_data(:,1),axis_data(:,2),xc,yc);
image_pixel_added=(axis_data(in_image,:))';%Pixels
xy=(mean(image_pixel_added'))';
if isfield(handles.parasaved, 'Pixels')
    handles.parasaved.Pixels=[handles.parasaved.Pixels,image_pixel_added];
    handles.parasaved.xy_all=[handles.parasaved.xy_all,xy];
    handles.parasaved.xypos=[handles.parasaved.xypos,xypos];
    delete(handles.htext);
else
    handles.parasaved.Pixels{1}=image_pixel_added;
    handles.parasaved.xy_all=xy;
    handles.parasaved.xypos{1}=xypos;
    handles.cellnum=1;
end
hline = plot(xc,yc,'r','LineWidth',2);
htext = zeros(1,size(handles.parasaved.xy_all,2));
htext=zeros(1,size(handles.parasaved.xy_all,2));
for k=1:size(handles.parasaved.xy_all,2)    
    htext(1,k)=text(handles.parasaved.xy_all(1,k),handles.parasaved.xy_all(2,k),num2str(k-1),'color','g','FontSize',12);
end

if isfield(handles, 'hline')
    handles.hline=[handles.hline,hline];
else
    handles.hline=hline;
end
    handles.htext=htext;

cell_no=(0:size(handles.parasaved.xy_all,2)-1);
set(handles.RoI_index,'String',num2str(cell_no'),'Value',1)    
else
     msgbox('Please choose ROI with area more than zero') ;
end
     guidata(hObject,handles);
  
end


% --- Executes on button press in Delete_pushbutton.
function Delete_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Delete_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
axes(handles.axes1);
[xpos,ypos] = ginput(1);
xypos=handles.parasaved.xypos;    

for i=1:length(xypos)
    in_image = inpolygon(xpos,ypos,xypos{i}(1,:),xypos{i}(2,:));
    if in_image
        ROIindexToDeleted=i;
        break;
    else
        ROIindexToDeleted=-1;
    end
end

if ROIindexToDeleted>0
   
    delete(handles.hline(ROIindexToDeleted));
    ROItoKeep=setdiff(1:length(xypos),ROIindexToDeleted);
    handles.parasaved.Pixels=handles.parasaved.Pixels(ROItoKeep);
    handles.parasaved.xy_all=handles.parasaved.xy_all(:,ROItoKeep);
    handles.parasaved.xypos=handles.parasaved.xypos(ROItoKeep);
    handles.hline=handles.hline(1,ROItoKeep);
    delete(handles.htext);
    htext=zeros(1,size(handles.parasaved.xy_all,2));
    
    for k=1:size(handles.parasaved.xy_all,2)    
        htext(1,k)=text(handles.parasaved.xy_all(1,k),handles.parasaved.xy_all(2,k),num2str(k-1),'color','g','FontSize',12);
    end
    handles.htext=htext;
    
   cell_no=(0:size(handles.parasaved.xy_all,2)-1);
   if ~isempty(cell_no)
       set(handles.RoI_index,'String',num2str(cell_no'),'Value',1)   
   else
       set(handles.RoI_index,'String',[]) ;    
   end
   
else
    msgbox('Please choose ROI within its region') ;
end
guidata(hObject,handles);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in RoI_index.
function RoI_index_Callback(hObject, eventdata, handles)
% hObject    handle to RoI_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
index_selected = get(hObject,'Value');
ShowPixels=handles.parasaved.Pixels(index_selected);
Showxy_all=handles.parasaved.xy_all(:,index_selected);
Showxypos=cell2mat(handles.parasaved.xypos(index_selected));
for si=1:5
showP=plot(Showxypos(1,:),Showxypos(2,:),'g');
pause(0.1)
delete(showP)
end

% Hints: contents = cellstr(get(hObject,'String')) returns RoI_index contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RoI_index


% --- Executes during object creation, after setting all properties.
function RoI_index_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RoI_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in AddCircle_pushbutton.
function AddCircle_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to AddCircle_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
axes(handles.axes1);
if ~isfield(handles,'ImageAdjust')
    msgbox('Please load a picure file or a picture&RoI file') ;
else
    ang=0:0.01:2*pi;
    r=str2double(get(handles.ROIradius,'String'));
    handles.parasaved.r=r;
    [xpos,ypos] = ginput(1);
    xc=xpos+r*cos(ang);
    yc=ypos+r*sin(ang);
    xc(xc>handles.IMGxrange)=handles.IMGxrange;
    yc(yc>handles.IMGyrange)=handles.IMGyrange;
    xc(xc<1)=1;
    yc(yc<1)=1;
    A = polyarea(xc,yc);
    if A>=1
        xypos=[xc;yc];%xypos
        [axisX axisY]=meshgrid(round(min(xc)):round(max(xc)),round(min(yc)):round(max(yc)));
        axis_data=[axisX(:) axisY(:)];
        in_image = inpolygon(axis_data(:,1),axis_data(:,2),xc,yc);
        image_pixel_added=(axis_data(in_image,:))';%Pixels
        xy=(mean(image_pixel_added'))';
        if isfield(handles.parasaved, 'Pixels')
            handles.parasaved.Pixels=[handles.parasaved.Pixels,image_pixel_added];
            handles.parasaved.xy_all=[handles.parasaved.xy_all,xy];
            handles.parasaved.xypos=[handles.parasaved.xypos,xypos];
            delete(handles.htext);
        else
            handles.parasaved.Pixels{1}=image_pixel_added;
            handles.parasaved.xy_all=xy;
            handles.parasaved.xypos{1}=xypos;
            handles.cellnum=1;
        end
        hline = plot(xc,yc,'r','LineWidth',2);
        htext = zeros(1,size(handles.parasaved.xy_all,2));
        
        for k=1:size(handles.parasaved.xy_all,2)
            htext(1,k)=text(handles.parasaved.xy_all(1,k),handles.parasaved.xy_all(2,k),num2str(k-1),'color','g','FontSize',12);
        end
        
        if isfield(handles, 'hline')
            handles.hline=[handles.hline,hline];
        else
            handles.hline=hline;
        end
        handles.htext=htext;
        
        cell_no=(0:size(handles.parasaved.xy_all,2)-1);
        %set(handles.RoI_index,'Value',[1:10]')
        set(handles.RoI_index,'String',num2str(cell_no'),'Value',1)
    else
        msgbox('Please choose ROI with area more than zero') ;
    end
    guidata(hObject,handles);
end



function ROIradius_Callback(hObject, eventdata, handles)
% hObject    handle to ROIradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ROIradius as text
%        str2double(get(hObject,'String')) returns contents of ROIradius as a double


% --- Executes during object creation, after setting all properties.
function ROIradius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROIradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SeedDetection_pushbutton.
function SeedDetection_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to SeedDetection_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);
axes(handles.axes1);

if ~isfield(handles,'ImageAdjust')
    msgbox('Please load a picure file or a picture&RoI file') ;
else
    if isfield(handles,'SeedsScatter')
       delete(handles.SeedsScatter);
       handles=rmfield(handles,'SeedsScatter');
    end
    if isfield(handles,'htext')
        delete(handles.htext);
        handles=rmfield(handles,'htext');
    end
    if isfield(handles,'hline')
       delete(handles.hline); 
       handles=rmfield(handles,'hline');
    end
    
    I_norm_ori=handles.I_adjust; %imagedata between 0 and 1;
    min_sigma=str2num(get(handles.MinSigma,'string'));
    max_sigma=str2num(get(handles.MaxSigma,'string'));
    Ntheta=handles.thetaN;
    C=0;
    tm = 0;
    p2 = 2*max_sigma*3+1;
    bw = adaptivethreshold(I_norm_ori,p2,C,tm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [new_coordinate_x,new_coordinate_y,absigma_leave]=SeedDetection(I_norm_ori,min_sigma,max_sigma,bw,Ntheta);
    handles.SeedsScatter=scatter(new_coordinate_x,new_coordinate_y,'r+');
    handles.new_coordinate_x = new_coordinate_x;
    handles.new_coordinate_y = new_coordinate_y;
    handles.absigma_leave = absigma_leave;
    handles.Ntheta=Ntheta;
    handles.min_sigma=min_sigma;
    handles.max_sigma=max_sigma;
    guidata(hObject,handles);
end



function MinSigma_Callback(hObject, eventdata, handles)
% hObject    handle to MinSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinSigma as text
%        str2double(get(hObject,'String')) returns contents of MinSigma as a double


% --- Executes during object creation, after setting all properties.
function MinSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxSigma_Callback(hObject, eventdata, handles)
% hObject    handle to MaxSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxSigma as text
%        str2double(get(hObject,'String')) returns contents of MaxSigma as a double


% --- Executes during object creation, after setting all properties.
function MaxSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AutoSegmentation.
function AutoSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to AutoSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
axes(handles.axes1);

if ~isfield(handles,'ImageAdjust')
    msgbox('Please load a picure file or a picture&RoI file') ;
elseif ~isfield(handles,'new_coordinate_x')
    msgbox('Please push the button ''Seed Detection''');
else
    I_norm_ori=handles.I_adjust; %imagedata between 0 and 1;
    MinArea=str2num(get(handles.MinArea,'string'));
    new_coordinate_x=handles.new_coordinate_x;
    new_coordinate_y=handles.new_coordinate_y;
    absigma_leave=handles.absigma_leave;
    Ntheta=handles.Ntheta;
    min_sigma=handles.min_sigma;
    ParametersOutput=fastSegmentation(I_norm_ori,new_coordinate_x,new_coordinate_y,absigma_leave,MinArea,Ntheta,min_sigma);
    hline=[];
    htext=[];
if iscell(ParametersOutput.xy_all)==1
    ParametersOutput.xy_all=cell2mat(ParametersOutput.xy_all);
end

if isfield(handles,'SeedsScatter')
   delete(handles.SeedsScatter)%delete the scatter
   handles=rmfield(handles,'SeedsScatter');
end
if isfield(handles,'htext')
   delete(handles.htext);
   handles=rmfield(handles,'htext');
end
if isfield(handles,'hline')
    delete(handles.hline); 
    handles=rmfield(handles,'hline');
end

for k=1:length(ParametersOutput.xypos)
    hline(1,k)=plot(ParametersOutput.xypos{k}(1,:),ParametersOutput.xypos{k}(2,:),'r','LineWidth',2);
 %   htext(1,k)=text(ParametersOutput.xy_all(1,k),ParametersOutput.xy_all(2,k),num2str(k-1),'color','g','FontSize',12);
end

handles.htext = htext;
handles.hline = hline;

cell_no=(0:length(ParametersOutput.xypos)-1);
set(handles.RoI_index,'String',num2str(cell_no'),'Value',1)
handles.parasaved=ParametersOutput;
guidata(hObject,handles); 
end


function MinArea_Callback(hObject, eventdata, handles)
% hObject    handle to MinArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinArea as text
%        str2double(get(hObject,'String')) returns contents of MinArea as a double


% --- Executes during object creation, after setting all properties.
function MinArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in loadVideo.
function loadVideo_Callback(hObject, eventdata, handles)
% hObject    handle to loadVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
axes(handles.axes1);
cla;
[filename pathname]=uigetfile('*.avi','Please select a avi file');
xyloObj = VideoReader([pathname filename]);
   nFrames = xyloObj.NumberOfFrames;
   vidHeight = xyloObj.Height;
   vidWidth = xyloObj.Width;
   mov_data = zeros(vidHeight,vidWidth,nFrames);%record movie data
  %  mov_data = zeros(460,500,nFrames);%record movie data
  
  hwaitbar = waitbar(0,'Loading video,please wait...');
  
  steps = nFrames;
  width=20;
  
  
   mov_data_temp=read(xyloObj,1);
   if length(size(mov_data_temp))==2
       for movi=1:nFrames
            waitbar(movi / steps)
            mov_data_temp=read(xyloObj,movi);
            mov_data(:,:,movi)=mov_data_temp;%(width+1:end-20,:);
       end
   else
       for movi=1:nFrames
            waitbar(movi / steps)
            mov_data_temp=rgb2gray(read(xyloObj,movi));
            mov_data(:,:,movi)=mov_data_temp;
       end
   end
close(hwaitbar); 
I=mean(double(mov_data),3);
I_adjust=double(I);
I_adjust=(I_adjust-min(I_adjust(:)))/(max(I_adjust(:))-min(I_adjust(:)));

if isfield(handles,'SeedsScatter')
   
   handles=rmfield(handles,'SeedsScatter');
end
if isfield(handles,'htext')
  
   handles=rmfield(handles,'htext');
end
if isfield(handles,'hline')
   try 
       delete(handles.hline);
   end
    handles=rmfield(handles,'hline');
end

imagesc(imadjust(I_adjust));
colormap(gray);
axis equal;
axis off;
hold on;
handles.I_adjust=I_adjust; %imagedata between 0 and 1;
handles.ImageAdjust=uint8(I_adjust*255);%image data between 0 and 255
handles.filename=filename;
handles.pathname=pathname;
handles.IMGxrange=size(I,2);
handles.IMGyrange=size(I,1);
set(handles.edit1,'string',[pathname filename]);
set(handles.edit1,'ForegroundColor',[1 0 0]);
set(handles.RoI_index,'string',[])
handles.parasaved=[];
handles.mov_data=mov_data;
guidata(hObject,handles);


% --- Executes on button press in plot_df_f.
function plot_df_f_Callback(hObject, eventdata, handles)
% hObject    handle to plot_df_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
Frequency = str2num(get(handles.Frequency,'string'));
Interval = str2num(get(handles.Interval,'string'));
mov_data = handles.mov_data;
ParametersOutput = handles.parasaved;
Pixels=ParametersOutput.Pixels;
nFrames=size(mov_data,3);
Cellnum=length(ParametersOutput.xypos);
F=zeros([nFrames Cellnum]);
vidHeight = handles.IMGyrange;
vidWidth  = handles.IMGxrange;
for fi=1:nFrames
    mov_cdata=mov_data(:,:,fi);
  for ci=1:Cellnum
      xi = Pixels{ci}(1,:);
      yi = Pixels{ci}(2,:);
      ind=sub2ind([vidHeight vidWidth],yi,xi);
      F(fi,ci)=mean(mov_cdata(ind));
  end
end
F_ave=mean(F);
y_labe=num2cell([0:5:Cellnum]);
 
 for cci=1:length(y_labe)
     y_labe{cci}=num2str(y_labe{cci});
 end
 
figure,
title('Calcium traces')
hold on;
ymax=0;
ymin=0;
dF_data=zeros(size(F));
for ii=1:Cellnum
    F_temp=(F(:,ii)-F_ave(ii))/F_ave(ii);
    plot((1:nFrames)/Frequency,F_temp+(ii)*Interval);
    dF_data(:,ii)=F_temp;
    ymin = min(min(F_temp+(ii)*Interval),ymin);
    ymax = max(max(F_temp+(ii)*Interval),ymax);
end

 xlabel('t(s)');
 ylabel('Cell NO');
 xmin = 0;
 xmax = nFrames/Frequency;
 ymin=0;
 axis([xmin xmax ymin ymax+Interval]);

 set(gca, 'YTick', [Interval Interval*5:Interval*5:Cellnum*Interval]);
 set(gca,'YTickLabel',y_labe);
 set(gca,'FontName','Times New Roman','FontSize',14);
 
 
 handles.F_data = F; %save F,not deltaF/F;
 handles.dF_data=dF_data;%save deltaF/F;
 guidata(hObject,handles);


% --- Executes on button press in save_df_f.
function save_df_f_Callback(hObject, eventdata, handles)
% hObject    handle to save_df_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
filename = handles.filename;
[FileName,PathName] = uiputfile([filename(1:end-4) '_Image_ROIdata.txt'],'Save Image and ROI');
data_output=handles.dF_data;
StrFormat=[];
for k=1:size(data_output,2)-1
    StrFormat=[StrFormat, '%5.8f '];
end
StrFormat=[StrFormat, '%5.8f'];

StrFormat=[StrFormat '\n'];
FileID = fopen([PathName '\' FileName],'w');
fprintf(FileID,StrFormat,data_output);
fclose(FileID);

[FileName2,PathName2] = uiputfile([filename(1:end-4) '_f_Image_ROIdata.txt'],'Save Image and ROI');
data_output2=handles.F_data;
StrFormat=[];
for k=1:size(data_output2,2)-1
    StrFormat = [StrFormat, '%5.8f '];
end
StrFormat=[StrFormat, '%5.8f'];
StrFormat=[StrFormat '\n'];
FileID2 = fopen([PathName2 '\' FileName2],'w');
fprintf(FileID2,StrFormat,data_output2);
fclose(FileID2);



function Interval_Callback(hObject, eventdata, handles)
% hObject    handle to Interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Interval as text
%        str2double(get(hObject,'String')) returns contents of Interval as a double


% --- Executes during object creation, after setting all properties.
function Interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Frequency_Callback(hObject, eventdata, handles)
% hObject    handle to Frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Frequency as text
%        str2double(get(hObject,'String')) returns contents of Frequency as a double


% --- Executes during object creation, after setting all properties.
function Frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
