function res = bound(v,p,max)

  n=size(v,1);
  for i=1:n,
    temp = p(i);
    if ( v(i) >= 0 && temp < 0 ),
      temp = -v(i)/temp;
      if( temp < max ),
        max = temp;
      end
    end
  end
  
  res= max;

