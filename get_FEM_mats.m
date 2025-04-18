function [S,M] = get_FEM_mats(V,F)
Nv = size(V,1);
Nf = size(F,1);
S = sparse(Nv,Nv);
M = sparse(Nv,Nv);

for face_idx = 1:Nf
    % face = F(face_idx,:);
    k1 = F(face_idx,1);
    k2 = F(face_idx,2);
    k3 = F(face_idx,3);
    e1=V(k2,:)-V(k3,:);
    e2=V(k3,:)-V(k1,:);
    e3=V(k1,:)-V(k2,:);
    % e1=V(k3,:)-V(k2,:);
    % e2=V(k1,:)-V(k3,:);
    % e3=V(k2,:)-V(k1,:);
    area = (0.5)*norm(cross(e1,e2));

    %update global stiffness
    S(k1,k1) = S(k1,k1) + (.25/area)*dot(e1,e1);
    S(k1,k2) = S(k1,k2) + (.25/area)*dot(e1,e2);
    S(k1,k3) = S(k1,k3) + (.25/area)*dot(e1,e3);
    
    S(k2,k1) = S(k2,k1) + (.25/area)*dot(e2,e1);
    S(k2,k2) = S(k2,k2) + (.25/area)*dot(e2,e2);
    S(k2,k3) = S(k2,k3) + (.25/area)*dot(e2,e3);

    S(k3,k1) = S(k3,k1) + (.25/area)*dot(e3,e1);
    S(k3,k2) = S(k3,k2) + (.25/area)*dot(e3,e2);
    S(k3,k3) = S(k3,k3) + (.25/area)*dot(e3,e3);

    %update global mass
    M(k1,k1)=M(k1,k1)+area/6;
    M(k1,k2)=M(k1,k2)+area/12;
    M(k1,k3)=M(k1,k3)+area/12;

    M(k2,k1)=M(k2,k1)+area/12;
    M(k2,k2)=M(k2,k2)+area/6;
    M(k2,k3)=M(k2,k3)+area/12;

    M(k3,k1)=M(k3,k1)+area/12;
    M(k3,k2)=M(k3,k2)+area/12;
    M(k3,k3)=M(k3,k3)+area/6;
end
end