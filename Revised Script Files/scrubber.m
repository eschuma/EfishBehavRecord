function sr=scrubber(input,output,sr)



%
%	scrubber.m
%	changes polarity, normalizes, aligns, and removes abnormal waveforms from TDT exported files
%       for use in Spike Sorting using wave_clus
%		returns a .mat file
%		can be used with cluster
%		creates the Files.txt file necessary for wave_clus
%
%	USAGE
%		sr=scrubber(input,output)
%			OR
%		sr=scrubber(input,output,sr)
%
%	VARIABLES
%		sr = sampling rate; default = 48828.125
%		input = name of the raw data input file; input is an .xlsx file where spikes are in columns H:AT and times are in column D; MUST be sorted by channel!
%		output = desired basename of the output .mat file; output is used as the input into wave_clus and file will be saved as basename_spikes.mat
%		spikes = normalized EOD waveforms; required for wave_clus
%		index = times the EODs were emitted; required for wave_clus
%



if (nargin==2),
	sr=48828.125;
end

spikes=xlsread(input,'H:AT');

disp('Done Reading In')


A=size(spikes,1);

%% to determine direction
[max_num, max_idx]=max(spikes,[],2);
[min_num, min_idx]=min(spikes,[],2);

for i=1:size(spikes,1)
	dir=min_idx(i,:)-max_idx(i,:);
	if dir<0
		spikes(i,:)=-(spikes(i,:));
	else
		spikes(i,:)=spikes(i,:);
	end
	spikes(i,:)=spikes(i,:)/(max_num(i,:)-min_num(i,:));
end

disp('Done Flipping and Normalizing')


%% to align spikes by max
spikes = [zeros(A,15),spikes]; %% adds to front of data
spikes = [spikes,zeros(A,15)]; %% adds to back of data

[max_num, max_idx2]=max(spikes,[],2);
r=max(max_idx2);

for i=1:size(spikes,1)
%	spikes(i,:)=circshift(spikes(i,:),[1 b(i,:)]);
	spikes(i,:)=circshift(spikes(i,:),[1 (r-max_idx2(i,:))]);
end


%% to remove all zero columns
spikes(:, find(sum(abs(spikes)) == 0)) = [];

disp('Done Aligning')


%% to get the best representative waveform for each timepoint
%% amplitudes for each channel
[x x2]=size(amp(:,1));
y=floor(x/8);
mamp(1:y,1)=amp(1:y,1);
mamp(1:y,2)=amp((y+1):(2*y),1);
mamp(1:y,3)=amp(((2*y)+1):(3*y),1);
mamp(1:y,4)=amp(((3*y)+1):(4*y),1);
mamp(1:y,5)=amp(((4*y)+1):(5*y),1);
mamp(1:y,6)=amp(((5*y)+1):(6*y),1);
mamp(1:y,7)=amp(((6*y)+1):(7*y),1);
mamp(1:y,8)=amp(((7*y)+1):(8*y),1);

[max_num, max_idx]=max(mamp,[],2);

for i=1:size(max_idx,1)
	max_spikes(i,:)=spikes((((max_idx(i,1)-1)*y)+i),:);
end

disp('Done Finding Respresentative EOD')

%% to only look at the max_spikes aka one for each time stamp
%clear spikes
%spikes=max_spikes;

%% if you need to modify the index variable to only include the time points once instead of for each channel
%times=xlsread(input,'D:D');
%index(1:y,:)=times(1:y,:)



%% to remove strange peaks via standard deviations; chosen by calculation
%% finding the average EOD and st dev at each time 
avgEOD = nanmean(spikes); 
sdEOD = std(spikes, 'omitnan'); 

%% generating upper and lower limits at each time 
upperEOD = avgEOD + 4*sdEOD; 
lowerEOD = avgEOD - 4*sdEOD; 

%% including curves that remain within st dev at all points
[m,n] = size(spikes); %m = number of rows, n = number of columns


for i = 1:m %% loops through rows
    for j = 1:n %% loops through columns
        if spikes(i,j) > upperEOD(1,j)
            spikes(i,:) = zeros;
            break
        elseif spikes(i,j) < lowerEOD(1,j)
        	spikes(i,:) = zeros; 
            break 
        end 
    end
end

disp('Done Removing Weird Waveforms')



%% to save for wave_clus
%%import data by chan in data matrix called spikes
%%import time into a column vector called index


index=xlsread(input,'D:D');
save(strcat(output,'_spikes.mat'),'index','spikes','sr','max_spikes')

%% creates the Files.txt file necessary to run wave_clus
fid=fopen('Files.txt','w');
fprintf(fid,strcat(output,'_spikes.mat'));
fclose(fid);