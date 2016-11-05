function acc=leggi_acc(nomefile,nlinee,ncol)
% legge un accelerogramma con dati equispaziati, disposti su ncol colonne e con nlinee inutili 

A=textread(nomefile,'','headerlines',nlinee);
siz=size(A)
k=0;
for i=1:siz(1),
    for j=1:ncol,
        acc(k+j)=A(i,j);
    end
    k=k+ncol;
end