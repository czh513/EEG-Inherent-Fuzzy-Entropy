clc;
clear all;
close all;

file_structure = dir('/home/zcao/Documents/migraine_xiandao_project/dataset/CMFM/');
m=size(file_structure,1);

channel = {  '1'  '2'  '3'  '4'   '5' '6'  '7'  '8'  '9'  '10'  '11' '12'  '13'  '14' '15' '16'  '17' '18' '19' '20'...
    '21' '22' '23' '24' '25' '26' '27' '28' '29'  '30' } ;
color=colormap(jet);

addpath('/home/zcao/Documents/MSE_ALL/fuzzy_MSE/');
addpath('/home/zcao/Documents/HHT_runcode/');


for ii=3:m
    
% clearvars -except file_structure channel color   ii
close all;

%% Pro-processing
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% ------------------------------------------------
location1=strcat('/home/zcao/Documents/migraine_xiandao_project/dataset/CMFM/',file_structure(ii).name,'/open/');
location2=strcat('/home/zcao/Documents/migraine_xiandao_project/dataset/CMFM/',file_structure(ii).name,'/close/');

% SAVES LOCATION
mkdir('/home/zcao/Documents/migraine_xiandao_project/results/CMFM/entropy/',file_structure(ii).name);
locationsave31=strcat('/home/zcao/Documents/migraine_xiandao_project/results/CMFM/entropy/',file_structure(ii).name,'/open_entropy');
locationsave32=strcat('/home/zcao/Documents/migraine_xiandao_project/results/CMFM/entropy/',file_structure(ii).name,'/close_entropy');

%location:open
file_structure2 = dir(location1);
m2=size(file_structure2,1);
for i1=3:m2
location_open{i1-2,:}=strcat(location1,file_structure2(i1).name);
end

%location:close
file_structure3 = dir(location2);
m3=size(file_structure3,1);
for i2=3:m3
location_close{i2-2,:}=strcat(location2,file_structure3(i2).name);
end

% open
for dd=1:3
EEG = pop_loadcnt( location_open{dd,1},'dataformat', 'auto','memmapfile', '');
EEG = eeg_checkset( EEG ); eeglab redraw;

EEG = pop_resample(EEG, 250);eeglab redraw;

EEG = pop_chanedit(EEG, 'lookup','/home/zcao/Documents/eeglab13_5_4b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
EEG = eeg_checkset( EEG );eeglab redraw;
EEG = pop_select( EEG,'nochannel',{'EKG(L)' 'EKG(R)'  'A1' 'A2' 'VEOU' 'VEOL'});
eeglab redraw;

EEG = pop_eegfilt( EEG, 0, 50,  0, 0, 'fir1', 0); %30-50Hz gramma
EEG = pop_eegfilt( EEG, 1, 0,   0, 0, 'fir1', 0);
eeglab redraw;

for m=1:30 % interplating channels
  e=find(EEG.data(m,:)==0)
  if ~isempty(e);
  EEG=eeg_interp(EEG,m);
end
end     

EEGtest = pop_rejcont(EEG, 'elecrange',[1:30] ,'freqlimit',[20 40] ,'threshold',10,'epochlength',0.5,'contiguous',4,'addlength',0.25,'taper','hamming');
newdata=double(EEGtest.data(:,end-13000:end-1));

eeglab redraw;
EEG_open(:,:,dd)=newdata;
end

% close
for dd=1:3
EEG = pop_loadcnt( location_close{dd,1},'dataformat', 'auto','memmapfile', '');
EEG = eeg_checkset( EEG ); eeglab redraw;

EEG = pop_resample(EEG, 250);eeglab redraw;

EEG=pop_chanedit(EEG, 'lookup','/home/zcao/Documents/eeglab13_5_4b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
EEG = eeg_checkset( EEG );eeglab redraw;
EEG = pop_select( EEG,'nochannel',{'EKG(L)' 'EKG(R)'  'A1' 'A2' 'VEOU' 'VEOL'});
eeglab redraw;

EEG = pop_eegfilt( EEG, 0, 50,  0, 0, 'fir1', 0);
EEG = pop_eegfilt( EEG, 1, 0,   0, 0, 'fir1', 0);
eeglab redraw;

for m=1:30 % interplating channels
  e=find(EEG.data(m,:)==0)
  if ~isempty(e);
  EEG=eeg_interp(EEG,m);
end
end     

EEGtest = pop_rejcont(EEG, 'elecrange',[1:30] ,'freqlimit',[20 40] ,'threshold',10,'epochlength',0.5,'contiguous',4,'addlength',0.25,'taper','hamming');
newdata=double(EEGtest.data(:,end-13000:end-1));

EEG_close(:,:,dd)=newdata;
end


%% Inherent Fuzzy Entropy

% EMD
for ddd=1:6
for n=1:30
  
datafile=EEG.data(n,:,ddd);
    %==============Run EEMD=================================  
    %give noise level
noiselevel=0;
%   give Number of ensemble 
Nensemble=1;
    %do eemd  and show result
        
EEMDIMF=eemd(datafile,noiselevel,Nensemble);

y=double(sum(EEMDIMF(:,4:10),2))';%1 is rawdata;
project_data=y;
EMD_data(n,:)=project_data;
end
EMD_EEG_data(:,:,ddd)=EMD_data;
end

% fuzzy entropy

m=2;
r=0.15;

for ddd=1:6
for nk=1:30

idata=squeeze(EMD_EEG_data(nk,:,ddd));

%normalization
y=idata;
y=y-mean(y);
st=std(y);  
y=y/st;

% fuzzy entropy
sampe = FuzzyEn(y,m,r*st,2);
IFE(nk,ddd)=sampe; % result

end
end

IFE_OE=mean(IFE(:,1:2:end),2);
IFE_CE=mean(IFE(:,2:2:end),2);


%save
save([locationsave31],'IFE_OE');
save([locationsave32],'IFE_CE');

end

