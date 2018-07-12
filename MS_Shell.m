clear;
close all;
clc;






%140117
cd('G:\HCPdata\Roma\script miei\Analisi Dati');
MS_SubjRunStruct={'140117','10-Motort',{'hcp_baddata'}};
save('MS_SubjRunStruct.mat','MS_SubjRunStruct');
MS_DataAnalysis_withouthcopy; clear; clc; close all;

cd('G:\HCPdata\Roma\script miei\Analisi Dati');
MS_SubjRunStruct={'140117','10-Motort',{'hcp_icaclass'}};
save('MS_SubjRunStruct.mat','MS_SubjRunStruct');
MS_DataAnalysis_withouthcopy; clear; clc; close all;

cd('G:\HCPdata\Roma\script miei\Analisi Dati');
MS_SubjRunStruct={'140117','11-Motort',{'hcp_baddata'}};
save('MS_SubjRunStruct.mat','MS_SubjRunStruct');
MS_DataAnalysis_withouthcopy; clear; clc; close all;

cd('G:\HCPdata\Roma\script miei\Analisi Dati');
MS_SubjRunStruct={'140117','11-Motort',{'hcp_icaclass'}};
save('MS_SubjRunStruct.mat','MS_SubjRunStruct');
MS_DataAnalysis_withouthcopy; clear; clc; close all;