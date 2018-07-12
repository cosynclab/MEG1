function [ vcomp ] = MS_plotComp( subj,run, ic )
load(['F:\Roma\output\',num2str(subj),'\',num2str(subj),'_MEG_',run,'_icaclass_vs.mat']);
sinfo=comp_class.sampleinfo;
trial=comp_class.trial;
lastpoint=sinfo(end,2);
vcomp=zeros(1,lastpoint);
sampleaxis=1:lastpoint;
for it1=1:size(sinfo,1)
   
    vcomp(sinfo(it1,1):sinfo(it1,2))=trial{it1}(ic,:);
    
end
plot (vcomp);
end

