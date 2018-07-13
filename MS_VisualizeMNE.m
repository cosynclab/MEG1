clear
clc
load('lista.mat');
MS_SubjectsListF=MS_SubjectsList;

MS_STARTSUBJ='105923';
MS_ENDSUBJ='';
MS_INCLUDESUBJ=[''];
MS_EXCLUDESUBJ=[''];

MS_STARTSUBJ=max([1,find(strcmp(MS_SubjectsList,MS_STARTSUBJ))]);
MS_ENDSUBJ=min([size(MS_SubjectsList,1),find(strcmp(MS_SubjectsList,MS_ENDSUBJ))]);
MS_SubjectsList=MS_SubjectsList(MS_STARTSUBJ:MS_ENDSUBJ);

for it1=1:size(MS_EXCLUDESUBJ,1)
    MS_indx=find(strcmp(MS_SubjectsList,MS_EXCLUDESUBJ(it1,:)));
    if (~isempty(MS_indx))
        MS_SubjectsList(MS_indx)=[];
    end
end

for it1=1:size(MS_INCLUDESUBJ,1)
    MS_indx=find(strcmp(MS_SubjectsList,MS_INCLUDESUBJ(it1,:)));
    if (isempty(MS_indx))
        MS_indx=find(strcmp(MS_SubjectsListF,MS_INCLUDESUBJ(it1,:)));
        if (~isempty(MS_indx))
            MS_SubjectsList=[MS_SubjectsList;{MS_INCLUDESUBJ(it1,:)}];
        end
    end
end

MS_DRIVE=cd;
MS_DRIVE=MS_DRIVE(1);
MS_Root=[MS_DRIVE,':/HCPdata/Roma'];
load('OPTIONS.mat');
load('FOLDERS.mat');
MS_DataFolder=[MS_Root,'/',MS_DataFolder];
MS_FTPath=[MS_Root,'/',MS_FTPath];
MS_HCPPath=[MS_Root,'/',MS_HCPPath];
MS_OutputFolder=[MS_Root,'/',MS_OutputFolder];
MS_ScriptFolder=[MS_Root,'/',MS_ScriptFolder];
x='1';
it1=1;
while(~strcmp(x,'q'))
    MS_SBJNM=num2str(MS_SubjectsList{it1});
    MS_ExperimentOutputPath=[MS_OutputFolder,'/',MS_SBJNM];
    MS_DIR=dir(MS_ExperimentOutputPath);
    MS_EXPID=[];
    for it2=1:size(MS_DIR,1)
        MS_NAME=MS_DIR(it2).name;
        [startIndex,endIndex] = regexp(MS_NAME,'MEG_\d[^_]*_icaclass_vs.mat');
        if(~isempty(startIndex)&&(size(MS_NAME,2)==endIndex))
            MS_EXPID=[MS_EXPID; {MS_NAME(startIndex+4:endIndex-16)}];
        end
    end
    it2=1;
    x='1';
    while((~strcmp(x(1),'g'))&&(~strcmp(x,'ns'))&&(~strcmp(x,'ps'))&&(~strcmp(x,'q')))
        disp([MS_SBJNM,' - ',MS_EXPID{it2}]);
        load([MS_ExperimentOutputPath,'/',MS_SBJNM,'_MEG_',MS_EXPID{it2},'_icaclass_vs.mat']);
        disp('TOTAL IC:');
        disp(comp_class.class.total_ic_number);
        disp('ECG/EOG IC:');
        disp(comp_class.class.ecg_eog_ic);
        artifact=[];
        for it3=1:comp_class.class.total_ic_number
            if((isempty(find(comp_class.class.ecg_eog_ic==it3)))&&(isempty(find(comp_class.class.brain_ic_vs==it3))))
                artifact=[artifact,it3];
            end
        end
        disp('ARTIFACT IC:');
        disp(artifact);
        disp('BRAIN IC:');
        disp(comp_class.class.brain_ic_vs);
        
        MS_good=[];
        for it3=1:comp_class.class.total_ic_number
            if((isempty(find(comp_class.class.brain_ic==it3)))&&(~isempty(find(comp_class.class.brain_ic_vs==it3))))
                MS_good=[MS_good,it3];
            end
        end
        
        disp('MANUAL BRAIN IC:');
        disp(MS_good);
        
        MS_bad=[];
        for it3=1:comp_class.class.total_ic_number
            if((~isempty(find(comp_class.class.brain_ic==it3)))&&(isempty(find(comp_class.class.brain_ic_vs==it3))))
                MS_bad=[MS_bad,it3];
            end
        end
        disp('MANUAL NOT BRAIN IC:');
        disp(MS_bad);
        
        it3=1;
        x='1';
        while((~strcmp(x(1),'g'))&&(~strcmp(x(1),'J'))&&(~strcmp(x,'n'))&&(~strcmp(x,'p'))&&(~strcmp(x,'ns'))&&(~strcmp(x,'ps'))&&(~strcmp(x,'q')))
%         for it3=1:comp_class.class.total_ic_number
            disp(['Component: ',num2str(it3),'/',num2str(comp_class.class.total_ic_number)]);
            
            MSCLASS='Brain IC';
            if(~isempty(find(comp_class.class.ecg_eog_ic==it3)))
                MSCLASS='ECG/EOG IC';
            elseif(~isempty(find(artifact==it3)))
                MSCLASS='ARTIFACT IC';
            end
            
            disp(['Classification: ',MSCLASS]);
            
            MS_F=openfig([MS_ExperimentOutputPath,'/',MS_SBJNM,'_MEG_',MS_EXPID{it2},'_icamne_',num2str(it3),'fig.fig']);
            set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
            MS_caxis=caxis;
            Deltaaxis=(MS_caxis(2)-MS_caxis(1))/100;
            pause(3);
            MS_IM1=imread([MS_ExperimentOutputPath,'/',MS_SBJNM,'_MEG_',MS_EXPID{it2},'_icamne_',num2str(it3),'.png']);
            MS_IM2=imread([MS_ExperimentOutputPath,'/',MS_SBJNM,'_MEG_',MS_EXPID{it2},'_icaclass_vs_',num2str(it3),'.png']);
            figure(2);
            imshow(MS_IM1);
            set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
            figure(3);
            imshow(MS_IM2);
            set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
            commandwindow;
            x = input('Command\n','s');
            
            while((~strcmp(x(1),'g'))&&(~strcmp(x(1),'j'))&&(~strcmp(x(1),'J'))&&(~strcmp(x,'a'))&&(~strcmp(x,'d'))&&(~strcmp(x,'n'))&&(~strcmp(x,'p'))&&(~strcmp(x,'ns'))&&(~strcmp(x,'ps'))&&(~strcmp(x,'q')))
                if(strcmp(x,'1')||strcmp(x,'2')||strcmp(x,'3'))
                    figure(str2num(x));
                elseif((length(x)>2)&&(strcmp(x(2),'+')))
                    
                    figure(1);
                    MS_caxis=caxis
                    Deltaaxis=(MS_caxis(2)-MS_caxis(1))/100
                    if(strcmp(x(1),'u'))
                        MS_axis=[MS_caxis(1),MS_caxis(2)+Deltaaxis*str2num(x(3:end))];
                    elseif(strcmp(x(1),'l'))
                        MS_axis=[MS_caxis(1)+Deltaaxis*str2num(x(3:end)),MS_caxis(2)];
                    end
                    MS_caxis1=min([MS_caxis(1),MS_caxis(2)]);
                    MS_caxis2=max([MS_caxis(1),MS_caxis(2)]);
                    MS_caxis=[max([MS_caxis1,0]),MS_caxis2];
                    caxis
                    caxis(MS_caxis);
                    
                elseif((length(x)>2)&&(strcmp(x(2),'-')))
                    
                    figure(1);
                    MS_caxis=caxis
                    Deltaaxis=(MS_caxis(2)-MS_caxis(1))/100
                    if(strcmp(x(1),'u'))
                        MS_axis=[MS_caxis(1),MS_caxis(2)-Deltaaxis*str2num(x(3:end))];
                    elseif(strcmp(x(1),'l'))
                        MS_axis=[MS_caxis(1)-Deltaaxis*str2num(x(3:end)),MS_caxis(2)];
                    end
                    MS_caxis1=min([MS_caxis(1),MS_caxis(2)]);
                    MS_caxis2=max([MS_caxis(1),MS_caxis(2)]);
                    MS_caxis=[max([MS_caxis1,0]),MS_caxis2];
                    caxis
                    caxis(MS_caxis);
                   
                end
                commandwindow;
                x = input('Command\n','s');
            end
            
            if(strcmp(x,'a'))
                it3=max([1,it3-1]);
            elseif((length(x)>1)&&(strcmp(x(1),'j')))
                it3=min([comp_class.class.total_ic_number,str2num(x(2:end))]);  
                it3=max([1,it3]);
            elseif((length(x)>1)&&(strcmp(x(1),'J')))
                it2=min([size(MS_EXPID,1),str2num(x(2:end))]);  
                it2=max([1,it2]);
            elseif(strcmp(x,'d'))
                it3=min([comp_class.class.total_ic_number,it3+1]);
            elseif(strcmp(x,'n'))
                it2=min([size(MS_EXPID,1),it2+1]);
            elseif(strcmp(x,'p'))
                it2=max([1,it2-1]);
            elseif(strcmp(x,'ns'))
                it1=min([size(MS_SubjectsList,1),it1+1]);
            elseif(strcmp(x,'ps'))
                it1=max([1,it1-1]);
            elseif((length(x)>2)&&(strcmp(x(1:2),'go')))
                [startIndexS,endIndexS] = regexp(x,'go_\d*');
                MS_subj=x(startIndexS+3:endIndexS);
                MS_subj=find(strcmp(MS_SubjectsList,MS_subj));
                if(~isempty(MS_subj))
                    it1=MS_subj;
                end
                
            end
            
            close(MS_F);
            
            
        end
    end
    
    
    
end
close all;    