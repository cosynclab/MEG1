outpath='G:\HCPdata\Roma\output\';
subj='105923';
list=dir([outpath,subj]);
files=[];
filesind=1;
for it=1:size(list,1)
   name=list(it).name;
   [s,e]=regexp(name,'icablpdyn');
   if(~isempty(s))
      files{filesind,1}=name;
      filesind=filesind+1;
   end
end

quit=0;
while(~quit)
   disp(files);
   cmd=lower(input('Command:','s'));
   
   switch cmd
   case 'quit'
      close all;
      quit=1;
   case 'stop'
      quit=1;
   case 'open'
      cmd=input('Filename:','s');
      ind=find(strcmp(files,cmd));
      if(~isempty(ind))
         close all;
         n=cmd(1:end-4);
         S=['Subject: ',subj, ' - '];
         [startIndex,endIndex] = regexp(n,'_MEG_.*_icablpdyn');
         R=['Run: ', n(startIndex+5:endIndex-10), ' - '];
         [startIndex,endIndex] = regexp(n,'dyn_.*_window');
         F=['Band: ' , n(startIndex+4:endIndex-7), ' - '];
         [startIndex,endIndex] = regexp(n,'length.*_timep');
         W=['WL: ' n(startIndex+6:endIndex-6), ' - '];
         [startIndex,endIndex] = regexp(n,'timepoint.*sconn');
         T=['TP: ', n(startIndex+9:endIndex-5)];
         TI=[];
         TT=[];
         [startIndex,endIndex] = regexp(n,'Motort');
         if(~isempty(startIndex))
                [startIndex,endIndex] = regexp(n,'TrialIndex_.*_TrialType');
                TI=[' - TI: ', n(startIndex+11:endIndex-10)];
                [startIndex,endIndex] = regexp(n,'TrialType_.*_windowlength');
                TT=[' - TI: ', n(startIndex+10:endIndex-13)];
                [startIndex,endIndex] = regexp(F,'.*_TrialIndex');
                F=[F(startIndex:endIndex-11),' - '];
         end
         tit=[S,R,F,W,T,TI,TT];
         load([outpath,subj,'\',cmd]);
         figure('units','normalized','outerposition',[0 0 1 1]);
         imagesc(conn.patched);
         colorbar;
         title([tit, ' - Patched']);
         figure('units','normalized','outerposition',[0 0 1 1]);
         imagesc(conn.parcelled);
         set(gca, 'XTick',[1:size(conn.parcelledOrdering,1)]);
         set(gca, 'XTickLabel',(conn.parcelledOrdering));
         set(gca, 'YTick',[1:size(conn.parcelledOrdering,1)]);
         set(gca, 'YTickLabel',(conn.parcelledOrdering));
         set(gca, 'XTickLabelRotation',90);
         colorbar;
         title([tit, ' - Parcelled']);
      end
    case 'change'
      cmd=input('Subject:','s');     
      subj=cmd;
      list=dir([outpath,subj]);
      files=[];
      filesind=1;
      for it=1:size(list,1)
         name=list(it).name;
         [s,e]=regexp(name,'icablpdyn');
         if(~isempty(s))
            files{filesind,1}=name;
            filesind=filesind+1;
         end
      end 
     case 'saveall'
       close all;
       for it=1:size(files,1)
           n=files{it};
           S=['Subject: ',subj, ' - '];
           [startIndex,endIndex] = regexp(n,'_MEG_.*_icablpdyn');
           R=['Run: ', n(startIndex+5:endIndex-10), ' - '];
           [startIndex,endIndex] = regexp(n,'dyn_.*_window');
           F=['Band: ' , n(startIndex+4:endIndex-7), ' - '];
           [startIndex,endIndex] = regexp(n,'length.*_timep');
           W=['WL: ' n(startIndex+6:endIndex-6), ' - '];
           [startIndex,endIndex] = regexp(n,'timepoint.*sconn');
           T=['TP: ', n(startIndex+9:endIndex-5)];
           TI=[];
           TT=[];
           [startIndex,endIndex] = regexp(n,'Motort');
           if(~isempty(startIndex))
                  [startIndex,endIndex] = regexp(n,'TrialIndex_.*_TrialType');
                  TI=[' - TI: ', n(startIndex+11:endIndex-10)];
                  [startIndex,endIndex] = regexp(n,'TrialType_.*_windowlength');
                  TT=[' - TI: ', n(startIndex+10:endIndex-13)];
                  [startIndex,endIndex] = regexp(F,'.*_TrialIndex');
                  F=[F(startIndex:endIndex-11),' - '];
           end
           tit=[S,R,F,W,T,TI,TT]; 
           load([outpath,subj,'\',n]);
           figure('units','normalized','outerposition',[0 0 1 1]);
           imagesc(conn.patched);
           colorbar;
           title([tit, ' - Patched']);
           savefig(gcf,[outpath,strrep(strrep([subj,'\',tit, ' - Patched.fig'],':',[]),' ', [])]);
           saveas(gcf,[outpath,strrep(strrep([subj,'\',tit, ' - Patched.png'],':',[]),' ', [])]);
           figure('units','normalized','outerposition',[0 0 1 1]);
           imagesc(conn.parcelled);
           set(gca, 'XTick',[1:size(conn.parcelledOrdering,1)]);
           set(gca, 'XTickLabel',(conn.parcelledOrdering));
           set(gca, 'YTick',[1:size(conn.parcelledOrdering,1)]);
           set(gca, 'YTickLabel',(conn.parcelledOrdering));
           set(gca, 'XTickLabelRotation',90);
           colorbar;
           title([tit, ' - Parcelled']);
           savefig(gcf,[outpath,strrep(strrep([subj,'\',tit, ' - Parcelled.fig'],':',[]),' ', [])]);
           saveas(gcf,[outpath,strrep(strrep([subj,'\',tit, ' - Parcelled.png'],':',[]),' ', [])]);
           close all;
       end
   end
    
    
end
