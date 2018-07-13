function [ conn ] = MS_CorrFun( BLPvect, CorrOpt )

SIndex=CorrOpt.SIndex;
EIndex=CorrOpt.EIndex;
Savepath=CorrOpt.Savepath;
FileSuff=CorrOpt.FileSuff;
SelectedVertices=CorrOpt.SelectedVertices;
EDist=CorrOpt.EDist;
Yeo17NetworksLabels=CorrOpt.YeoNetworks.Yeo17NetworksLabels;
Yeo17NetworksVertices=CorrOpt.YeoNetworks.Yeo17NetworksVertices;
sourcemodel2d=CorrOpt.SourceModel;
Wind=CorrOpt.Wind;
returnDense=CorrOpt.returnDense;
returnPatched=CorrOpt.returnPatched;
returnParcelled=CorrOpt.returnParcelled;


NWind=floor((EIndex-SIndex+1)/Wind);
connect_stat=zeros(size(BLPvect,1));

for it1=1:NWind
    vect=[((it1-1)*Wind)+1:Wind*it1];
    connect_stat=corr(BLPvect(:,vect)')+connect_stat;
end
connect_stat = connect_stat/NWind;

eudist=squareform(pdist(sourcemodel2d.pos));
mask_c=find(eudist<EDist);
connect_stat(mask_c)=NaN;
conn.mask=mask_c;
conn.complete=connect_stat;




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
conn.patchLab=unique(ord(:,2:4),'rows');




patches=max(ord(:,4));

patched=zeros(patches);
for it1=1:patches
    for it2=1:patches
        p1ind=find(ord(:,4)==it1);
        p2ind=find(ord(:,4)==it2);
        squareV=conn.complete(ord(p1ind,1),ord(p2ind,1));
        
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

ordering_n=[17 18 16 4 5 15 9 12 13 14 7 6 8 10 11 2 3 1];
conn.NetLabels={ 'MedialWall'
        '1:VIS1'
        '2:VIS2'
        '3:MOT1'
        '4:MOT2'
        '5:DAN2'
        '6:DAN1'
        '7:VAN1'
        '8:FP1'
        '9:LIM1'
        '10:LIM2'
        '11:FP2'
        '12:FP3'
        '13:FP4'
        '14:MOT3'
        '15:DMN3'
        '16:DMN1'
        '17:DMN2'};
    
    
conn.parcelled2=conn.parcelled(ordering_n,ordering_n);

if(returnParcelled)
conn2.parcelled=conn.parcelled2;
conn2.parcelledOrdering=conn.NetLabels;

end

if(returnPatched)
conn2.patched=conn.patched;
conn2.patchedOrdering=conn.patchLab;
end

if(returnDense)
conn2.dense=conn.ordered;
conn2.denseOrdering=ord;
end


conn=conn2;
conn.NWind=NWind;

save([Savepath,FileSuff,'conn.mat'],'conn');

end

