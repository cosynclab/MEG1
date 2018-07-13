clear, clc
load('lista.mat');
MS_SpecSubjRun=1;

if (MS_SpecSubjRun)
    load('MS_SubjRunStruct.mat');
    MS_STARTSUBJ=MS_SubjRunStruct{1};
    MS_ENDSUBJ=MS_SubjRunStruct{1};
    MS_INCLUDESUBJ=[''];
    MS_EXCLUDESUBJ=[''];
    MS_EXPIDS_manual=MS_SubjRunStruct(2);

else

    

    MS_STARTSUBJ='';
    MS_ENDSUBJ='';
    MS_INCLUDESUBJ=[''];
    MS_EXCLUDESUBJ=[''];

    MS_EXPIDS_manual={''};
end

    MS_SubjectsListF=MS_SubjectsList;
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



MS_DIR=dir;
MS_DRIVE=MS_DIR(1).folder;
MS_DRIVE=MS_DRIVE(1);
MS_Root=[MS_DRIVE,':/Roma'];

load('OPTIONS.mat');
load('YeoNetworks.mat');
MS_VERTICES_SELECTION=[];
for it1=1:size(Yeo17NetworksVertices,1)
    if(Yeo17NetworksVertices(it1,2)~=0)
        MS_net=Yeo17NetworksLabels{Yeo17NetworksVertices(it1,2),1};
        MS_subnet=str2num(Yeo17NetworksLabels{Yeo17NetworksVertices(it1,2),2});
        MS_indNet=[];
        for it2=1:size(MS_NETWORKS,1)
            if(strcmp(MS_NETWORKS{it2,1},MS_net))
                MS_indNet=it2;
            end
        end
        if (~isempty(MS_indNet))
            if(~isempty(MS_NETWORKS{MS_indNet,2}))
                if((strcmp(MS_NETWORKS{MS_indNet,2},'all'))||(~isempty(find(MS_NETWORKS{MS_indNet,2}==MS_subnet))))
                    MS_VERTICES_SELECTION=[MS_VERTICES_SELECTION;Yeo17NetworksVertices(it1,1)];
                end
            end
        end
    end
end

load('FOLDERS.mat');
MS_DataFolder=[MS_Root,'/',MS_DataFolder];
MS_FTPath=[MS_Root,'/',MS_FTPath];
MS_HCPPath=[MS_Root,'/',MS_HCPPath];
MS_OutputFolder=[MS_Root,'/',MS_OutputFolder];
MS_ScriptFolder=[MS_Root,'/',MS_ScriptFolder];
cd(MS_FTPath);
ft_defaults;

addpath(genpath(MS_HCPPath));
addpath(MS_DataFolder);
addpath(MS_OutputFolder);
addpath(MS_ScriptFolder);
cd(MS_ScriptFolder);

for it1=1:size(MS_SubjectsList,1)
    clear subjectid
    clear experimentid
    clear scanid
    MS_FILE=fopen([MS_ScriptFolder,'/MS_log.txt'],'a');
    MS_SBJNM=num2str(MS_SubjectsList{it1});
    fprintf(MS_FILE,'%s\t%s\n',['Inizio Analisi Soggetto: ',MS_SBJNM],datestr(now,'yyyy-mm-dd-HH-MM-SS'));
    MS_SBJDataFolder=[MS_DataFolder,'/',MS_SBJNM,'/MEG'];
    MS_ExperimentOutputPath=[MS_SBJNM];
    if(~exist([MS_OutputFolder,'/',MS_ExperimentOutputPath],'dir'))
        mkdir([MS_OutputFolder,'/',MS_ExperimentOutputPath]);
    end
    MS_ExperimentOutputPath=[MS_OutputFolder,'/',MS_ExperimentOutputPath];
    pipelinedatadir=MS_ExperimentOutputPath;
    MS_DIR1=dir(MS_SBJDataFolder);
    MS_EXPIDS={};
    MS_itEXPIDS=1;
    for it2=3:size(MS_DIR1,1)
        if((MS_DIR1(it2).isdir)&&(~strcmp(MS_DIR1(it2).name,'anatomy')))
            MS_DIR2=dir([MS_DIR1(it2).folder,'/',MS_DIR1(it2).name]);
            for it3=3:size(MS_DIR2,1)
                if(MS_DIR2(it3).isdir)
                    MS_DIR3=dir([MS_DIR2(it3).folder,'/',MS_DIR2(it3).name]);
                    for it4=3:size(MS_DIR3,1)
                        if(~MS_DIR3(it4).isdir)
                            if(~exist([MS_ExperimentOutputPath,'/',MS_DIR3(it4).name]))
                                copyfile([MS_DIR3(it4).folder,'/',MS_DIR3(it4).name], MS_ExperimentOutputPath);
                            end
                            [startIndex,endIndex] = regexp(MS_DIR3(it4).name,'MEG_\d[^_]*_');
                            if ((isempty(MS_EXPIDS))&&(~isempty(MS_DIR3(it4).name(startIndex+4:endIndex-1))))
                                MS_EXPIDS{MS_itEXPIDS,1}=MS_DIR3(it4).name(startIndex+4:endIndex-1);
                                MS_itEXPIDS=MS_itEXPIDS+1;
                            else
                                MS_FOUND=0;
                                for it5=1:size(MS_EXPIDS,1)
                                    if(strcmp(MS_DIR3(it4).name(startIndex+4:endIndex-1),MS_EXPIDS{it5}))
                                       MS_FOUND=1;
                                    end
                                end
                                if ((~MS_FOUND)&&(~isempty(MS_DIR3(it4).name(startIndex+4:endIndex-1))))
                                    MS_EXPIDS{MS_itEXPIDS,1}=MS_DIR3(it4).name(startIndex+4:endIndex-1);
                                    MS_itEXPIDS=MS_itEXPIDS+1;
                                end
                            end
                        else
                            if(strcmp(MS_DIR3(it4).name,'figures'))
                                MS_DIR4=dir([MS_DIR3(it4).folder,'/',MS_DIR3(it4).name]);
                                for it5=3:size(MS_DIR4,1)
                                    if(~MS_DIR4(it5).isdir)
                                        if(~exist([MS_ExperimentOutputPath,'/',MS_DIR4(it5).name]))
                                            copyfile([MS_DIR4(it5).folder,'/',MS_DIR4(it5).name], MS_ExperimentOutputPath);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        if((MS_DIR1(it2).isdir)&&(strcmp(MS_DIR1(it2).name,'anatomy')))
            MS_DIR2=dir([MS_DIR1(it2).folder,'/',MS_DIR1(it2).name]);
            for it3=3:size(MS_DIR2,1)
                if(~MS_DIR2(it3).isdir)
                    if(~exist([MS_ExperimentOutputPath,'/',MS_DIR2(it3).name]))
                        copyfile([MS_DIR2(it3).folder,'/',MS_DIR2(it3).name], MS_ExperimentOutputPath);
                    else
                        if(strcmp(MS_DIR2(it3).name,'figures'))
                            MS_DIR3=dir([MS_DIR2(it3).folder,'/',MS_DIR2(it3).name]);
                            for it4=3:size(MS_DIR3,1)
                                if(~MS_DIR3(it4).isdir)
                                    if(~exist([MS_ExperimentOutputPath,'/',MS_DIR3(it4).name]))
                                        copyfile([MS_DIR3(it4).folder,'/',MS_DIR3(it4).name], MS_ExperimentOutputPath);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    fclose(MS_FILE);
    if(~isempty(MS_EXPIDS_manual))
        disp('Run manuali:');
        disp(MS_EXPIDS_manual);
        MS_EXPIDS=MS_EXPIDS_manual;
    end
    for it2=1:size( MS_EXPIDS,1)
        MS_FILE=fopen([MS_ScriptFolder,'/MS_log.txt'],'a');
        fprintf(MS_FILE,'\t%s\t%s\n',['Inizio esperimento:',MS_EXPIDS{it2}],datestr(now,'yyyy-mm-dd-HH-MM-SS'));
        filename=[MS_DataFolder,'/',MS_SBJNM,'/Experiments/',MS_SBJNM,'_MEG/Scans/',MS_EXPIDS{it2},'/Resources/4D/c,rfDC'];
        clear subjectid
        clear experimentid
        clear scanid
        
        if (MS_SpecSubjRun)
              buffunc=MS_SubjRunStruct{3};
              funstr='';
              for it3=1:size(buffunc,1)
                 funstr=[funstr,buffunc{it3},'; ']; 
              end
              
              fprintf(MS_FILE,'\t\t%s\t%s\n',['Script da eseguire:',funstr],datestr(now,'yyyy-mm-dd-HH-MM-SS'));
              for it3=1:size(buffunc,1)
                    fprintf(MS_FILE,'\t\t%s\t%s\n',['Inizio Script:',buffunc{it3}],datestr(now,'yyyy-mm-dd-HH-MM-SS'));
                    eval([buffunc{it3},';']);
                    fprintf(MS_FILE,'\t\t%s\t%s\n',['Script: ',buffunc{it3},' terminato'],datestr(now,'yyyy-mm-dd-HH-MM-SS'));      
              end
                
        else

% 
%             hcp_baddata;
%             hcp_icaclass;
%             hcp_icaclass_qc;
%             hcp_tmegpreproc;
%             hcp_icamne;
% 
%             MS_hcp_icablpenv;
%             hcp_icablpenv;

            
        end

        fprintf(MS_FILE,'\t%s\t%s\n',['Esperimento: ',MS_EXPIDS{it2},' terminato'],datestr(now,'yyyy-mm-dd-HH-MM-SS'));
        fclose(MS_FILE);
    end
    
end
