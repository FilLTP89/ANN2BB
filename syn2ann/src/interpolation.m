function [Y_m] = interpolation(X,Y,xi)


if length(find(X==xi))==0
   r_low=sum(X<xi);
   m=(Y(r_low+1)-Y(r_low))/(X(r_low+1)-X(r_low));
   q=Y(r_low)-m*X(r_low);
   Y_m=m*xi+q;
else
    i=find(X==xi)
    Y_m=Y(i);
end

end

