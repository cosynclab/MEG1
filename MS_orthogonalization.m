function [ SIGort ] = orthogonalization( MIS, SIG )


[L1,L2]=size(MIS);
TRENDS=SIG;

for it1=1:L1
    
    SIG(it1,:)=detrend(SIG(it1,:));
    
    TRENDS(it1,:)=TRENDS(it1,:)-SIG(it1,:);
    
    
end



MISnorm=sqrt(sum((MIS).^2,2));
MISunity=MIS./repmat(MISnorm,1,L2);
SIGproj=sum(SIG.*MISunity ,2);
SIGort=SIG-MISunity.*repmat(SIGproj,1,L2);
SIGort=SIGort+TRENDS;

end

