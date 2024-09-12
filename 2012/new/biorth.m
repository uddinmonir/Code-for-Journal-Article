function [V,W] = biorth(V,W)
[n,r]=size(V);

for i=1:r
    
    v=V(:,i); w=W(:,i); 
    p1=1; p2=1;
    
    for k=1:i
        
     p1=p1*(eye(n)-V(:,k)*W(:,k)');   
     p2=p2*(eye(n)-W(:,k)*V(:,k)');   
        
    end
    v=p1*v; w=p2*w; 
    v=v/norm(v); w=w/norm(w); 
    v=v/(w'*v);
    V(:,i)=v; W(:,i)=w;
    
end

end

