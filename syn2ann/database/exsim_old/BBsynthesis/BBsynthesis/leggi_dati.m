function dati=leggi_dati(nomedati)
D=textread(nomedati,'%s');
s=size(D);
for i=1:s(1),
    dati(i)=D(i);
end