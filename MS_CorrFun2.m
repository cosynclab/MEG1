function [ conn ] = MS_CorrFun2( connect_stat, CorrOpt )

SIndex=CorrOpt.SIndex;
EIndex=CorrOpt.EIndex;
Savepath=CorrOpt.Savepath;
FileSuff=CorrOpt.FileSuff;
SelectedVertices=CorrOpt.SelectedVertices;
EDist=CorrOpt.EDist;
load(CorrOpt.YeoNetworks);
load(CorrOpt.SourceModelPath);
Wind=CorrOpt.Wind;




conn.complete=connect_stat;

nanind=find(~isnan(connect_stat(1,:)));
withoutNAN=connect_stat(nanind,nanind);


conn.withoutNAN=withoutNAN;



for it1=1:size(Yeo17NetworksLabels,1)
    
    if(strcmp(Yeo17NetworksLabels{it1,1},'MD'))
        Yeo17NetworksLabels{it1,3}=0;
    else
        Yeo17NetworksLabels{it1,3}=str2num(Yeo17NetworksLabels{it1,1}(2:end));
    end
    Yeo17NetworksLabels{it1,2}=str2num(Yeo17NetworksLabels{it1,2});
end
Yeo17NetworksLabelsIND=cell2mat(Yeo17NetworksLabels(:,2:3));

v=Yeo17NetworksVertices(SelectedVertices,:);

ordered=zeros(size(SelectedVertices,1));
ord=zeros(size(SelectedVertices,1),4);
nextind=1;
incind=1;
for it1=0:17
   indb=find(Yeo17NetworksLabelsIND(:,2)==it1);
   redLab=Yeo17NetworksLabelsIND(indb,:);
   for it2=1:max(redLab(:,1))
       indbb=find(redLab(:,1)==it2);
       if(~isempty(indbb))
           finind=indb(indbb);
           vindex=find(v(:,2)==finind);
           ll=size(vindex,1);
           ord(nextind:nextind+ll-1,1)=v(vindex,1);

           ord(nextind:nextind+ll-1,2)=it1;
           ord(nextind:nextind+ll-1,3)=it2;

           ord(nextind:nextind+ll-1,4)=incind;

           nextind=nextind+ll;
           incind=incind+1;
       end
   end
   
end

ordered(:,:)=conn.complete(ord(:,1),ord(:,1));

conn.ordered=ordered;
conn.ord=ord;

eudist=squareform(pdist(sourcemodel2d.pos));
mask_c=find(eudist<EDist);
connect_stat(mask_c)=NaN;



patches=max(ord(:,4));

patched=zeros(patches);
for it1=1:patches
    for it2=1:patches
        p1ind=find(ord(:,4)==it1);
        p2ind=find(ord(:,4)==it2);
        squareV=connect_stat(ord(p1ind,1),ord(p2ind,1));
        
        patched(it1,it2)=nanmean(nanmean(squareV));
    
    end
end
conn.patched=patched;


parcels=max(ord(:,2))+1;
parcelled=zeros(parcels);
for it1=0:parcels-1
    for it2=0:parcels-1
        p1ind=find(ord(:,2)==it1);
        p2ind=find(ord(:,2)==it2);
        squareV=conn.complete(ord(p1ind,1),ord(p2ind,1));
        
        parcelled(it1+1,it2+1)=nanmean(nanmean(squareV));
    
    end
end

conn.parcelled=parcelled;


conn.NetLabels={ 'MedialWall'
        '17Networks1VIS1'
        '17Networks2VIS2'
        '17Networks3MOT1'
        '17Networks4MOT2'
        '17Networks5DAN2'
        '17Networks6DAN1'
        '17Networks7VAN1'
        '17Networks8FP1'
        '17Networks9LIM1'
        '17Networks10LIM2'
        '17Networks11FP2'
        '17Networks12FP3'
        '17Networks13FP4'
        '17Networks14MOT3'
        '17Networks15DMN3'
        '17Networks16DMN1'
        '17Networks17DMN2'};
    
    



save([Savepath,FileSuff,'conn.mat'],'conn');

end

