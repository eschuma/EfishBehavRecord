function [Fish1,Fish2,Ipi1,Ipi2,NoFish]=fish_sort_nfish(input)



%
%	fish_sort_nfish.m
%	calls wave_clus for spike sorting; takes results of wave_clus and separates them by fish and returns
%       for use in Spike Sorting using wave_clus
%		can be used with cluster
%
%	USAGE
%		[Fish1...FishN,Ipi1...IpiN,NoFish]=fish_sort(input)
%
%	VARIABLES
%		input = basename of the scrubber output_spikes.mat file; i.e. what comes before the _spikes.mat; output file will be saved as basename_sorted.mat
%		FishN = timestamps of EODs produced by fish1N
%		IpiN = interpulse intervals of fishN
%		NoFish = timestamps of signals captured but not attributed to either fish; i.e. noise
%


%% batch clustering using wave_clus
%Do_clustering('Files.txt')

%% combining channels, if sort together
load(strcat('times_',input,'.mat'))

%% replace rows of all 0s with NaN
for i=1:size(spikes,1)
	if sum(spikes(i,:))==0
		cluster_class2(i,1)=NaN;
	else
		cluster_class2(i,1)=cluster_class(i,1);
	end
end
cluster_class2(:,2)=cluster_class(:,2);

%% start commenting out here if you only want to look at max_spikes, see line 90 for more details

%% makes a table of the channels
[x x2]=size(cluster_class(:,1));
y=floor(x/8);
mc(1:y,1)=cluster_class2(1:y,2);
mc(1:y,2)=cluster_class2(1:y,1);
mc(1:y,3)=cluster_class2((y+1):(2*y),1);
mc(1:y,4)=cluster_class2(((2*y)+1):(3*y),1);
mc(1:y,5)=cluster_class2(((3*y)+1):(4*y),1);
mc(1:y,6)=cluster_class2(((4*y)+1):(5*y),1);
mc(1:y,7)=cluster_class2(((5*y)+1):(6*y),1);
mc(1:y,8)=cluster_class2(((6*y)+1):(7*y),1);
mc(1:y,9)=cluster_class2(((7*y)+1):(8*y),1);

%% determing sort across all 8 channels
%% to get sort results consolidated between all channels, sort is resultant vector
sort(:,1)=mc(:,1);
sort(:,2)=mode(mc,2);


%% makes timestamp NaN if there is not a strong majority across all eight channels or if most channels are removed from analysis due to breaking standard deviation
for i=1:size(sort,1)
	if sort(i,1)==sort(i,2)
		sort(i,2)=0;
	end	
	for j=1:max(sort(:,2))
		if sort(i,2)==j
			Z(i,:)=sum(mc(i,:)==j);
			Y(i,:)=8-Z(i,:)-sum(mc(i,:)==0);
			if Z(i,:)<3 | Y(i,:)>=(Z(i,:)/2)-1				
				sort(i,2)=0;
			end
			break
		else
			continue
		end
	end
end

%% histogram of proportion of channels that agree with sort code (mode) from raw sort data; i.e. no correction for NaN's or rows with <3 channels not NaN or where it was half (or nearly half) 1s and half 2s
%for i=1:size(Z,1)
%	W(i,:)=8-(Z(i,:)+Y(i,:));
%	X(i,:)=(abs(Z(i,:)-Y(i,:)))/(8-W(i,:));
%end
%figure
%histogram(X) %% proportion of channels that agree with the majority
%print('mode_consensus_hist','-dpng')
%figure
%histogram(W) %% number of channels removed from analysis
%print('chan_remove_hist','-dpng')

%% if you want to only use the max spikes: uncomment this make sure that spikes=max_spikes in the end scrubber.m before the file is saved
%sort=cluster_class2;

%% to consider NaN's
%% to separate fish into separate tables, including NaN's
for j=1:max(sort(:,2))
	for i=1:size(sort,1)
		if sort(i,2)==j | sort(i,2)==0
			fishNaN(i,:)=sort(i,:);
		end
	end
	
	fishNaN(all(fishNaN==0,2),:)=[];
	fishNaN=sortrows(fishNaN);
	
	%%%%%%%%%%%%%%%%%%%% Comp Based NaN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	%% comparing derivatives using with NaN and without NaN, whichever is smaller should be sort value
	v=1; %% this is the threshold value used in determining whether the questionable EOD should be included or not

	for i=3:size(fishNaN,1)-2
		if fishNaN(i,2)~=j
			while fishNaN(i+1,1)~=j | fishNaN(i+2,1)~=j
			
			
				E(i,:)=fishNaN(i+1,1)-fishNaN(i,1); %% these are ipis, need to calculate ipi'
				F(i,:)=fishNaN(i,1)-fishNaN(i-1,1);
				G(i,:)=fishNaN(i+1,1)-fishNaN(i-1,1); %% this is without considering NaN
				H(i,:)=fishNaN(i+2,1)-fishNaN(i+1,1); 
				I(i,:)=fishNaN(i-1,1)-fishNaN(i-2,1);
				J(i,1)=abs(I(i,1)-F(i,1));
				J(i,2)=abs(F(i,1)-E(i,1));
				J(i,3)=abs(E(i,1)-H(i,1));
				K(i,1)=abs(I(i,1)-G(i,1));
				K(i,2)=abs(G(i,1)-H(i,1));
				%[Mi,max_idx]=max(J,[],2);
				%[Mr,m_idx]=max(K,[],2);
				%Mi=median(J,2);
				%Mr=median(K,2);
				%Mi=var(J,0,2);
				%Mr=var(K,0,2);
				%Mi=min(J,[],2);
				%Mr=min(K,[],2);
				Mi=mean(J,2);
				Mr=mean(K,2);
				%maj=max(J,[],2);
				%mij=min(J,[],2);
				%mak=max(K,[],2);
				%mik=min(K,[],2);
				%Mi=maj/mij;
				%Mr=mak/mik;
				%Mi=max(J,[],2)-min(J,[],2);
				%Mr=max(K,[],2)-min(K,[],2);
				if Mr(i,:) > v*Mi(i,:)
					fishNaN(i,2)=j;
				else
					fishNaN(i,2)=0;
				end
			end
		end
	end

	clear E F G H I Mi Mr

	
	%% rewrite sort from combining fish1NaN and fish2NaN
	if j==1
		sort2=[fishNaN];
		for i=1:size(fishNaN,1)
			if fishNaN(i,2)~=j
				notf(i,:)=fishNaN(i,:);
			else
				Fish(i,:)=fishNaN(i,:);
			end
		end
		notf(all(notf==0,2),:)=[];
		NoFish=notf;
		Fish(all(Fish==0,2),:)=[];
		Fish=sortrows(Fish);
	else
		sort2=[sort2; fishNaN];
		for i=1:size(fishNaN,1)
			if fishNaN(i,2)~=j
				notf(i,:)=fishNaN(i,:);
			else
				Fish(i,:)=fishNaN(i,:);
			end
		end
		notf(all(notf==0,2),:)=[];
		NoFish=intersect(NoFish,notf,'rows');
		Fish(all(Fish==0,2),:)=[];
		Fish=sortrows(Fish);	
	end


	%% finding interpulse intervals
	for i=1:(size(Fish,1))-1
		Ipi(i,:)=Fish(i+1,1)-Fish(i,1);
	end

	%% plotting
	figure
	plot(Fish(1:length(Fish)-1,1),Ipi,'.')
	print(strcat('Fish',int2str(j)),'-dpng')
	
	%% saving for Fish j and Ipi j
	eval(strcat('Fish',int2str(j),'=Fish;'));
	eval(strcat('Ipi',int2str(j),'=Ipi;'));
	
	%% to save the file
	if j==1
		outfile=strcat('Fish',int2str(j),',Ipi',int2str(j));
	else
		outfile=strcat(outfile,',Fish',int2str(j),',Ipi',int2str(j));
	end
	
	clear fishNaN Fish Ipi notf
end


output=strcat(input,'_sorted.mat')
save(output,outfile,'NoFish','max_spikes')
