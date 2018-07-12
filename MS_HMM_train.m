clear;
addpath(genpath('G:\HCPdata\Roma\HMM-MAR-master'))
subjects={'105923','106521','108323','109123','113922','116726','133019'};
data=[];
for itsubj=1:7

subject=subjects{itsubj};
folder=['G:\HCPdata\Roma\output\',subject];
listFile=dir(folder);

runs={'3-Restin','4-Restin','5-Restin'};

wl=10;
wstep=0.5;
firsttp=wl/2;

states=12;


for itfile=1:size(listFile,1);

    buffName=listFile(itfile).name;
    [startIndex,endIndex] = regexp(buffName,'MEG_.*_icablpdyn_alpha_windowlength');
    
    if(~isempty(startIndex))
        [startIndex2,endIndex2] = regexp(buffName,'._icablpdyn_alpha_windowlength');
        RUN=buffName(startIndex+4:startIndex2);
        
        [startIndex3,endIndex3] = regexp(buffName,'timepoint.*sconn.mat');
        timepoint=str2num(buffName(startIndex3+9:endIndex3-9));

        for itrun=1:size(runs,2)
            if(strcmp(RUN,runs{itrun}))
                load([folder,'\',buffName]);
                bufcon=conn.parcelled;
                bufcon=bufcon(find(triu(bufcon)));
                tpindex=round((timepoint-firsttp)/wstep)+1;
                data{itrun+(itsubj-1)*3}(:,tpindex)=bufcon;
                
                
            end
        end
    
    end
    
end
end
K=12;
ndim=size(data{1},1);
N=3;
Fs=1/wstep;
T=[];
for it1=1:(size(data,2))
	T{it1}=size(data{it1},2);
end

options = struct();
options.K = K; 
options.Fs = Fs; 
options.covtype = 'full';
options.order = 0;
options.DirichletDiag = 2; 
options.zeromean = 0;
options.verbose = 1;
for it1=1:(size(data,2))
	data{it1}=data{it1}';
end

[hmm, Gamma, Xi, vpath] = hmmmar(data,T,options);
